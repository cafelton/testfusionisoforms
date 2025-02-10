from comparefusiondetection_Jan2025 import *
import matplotlib.pyplot as plt
import subprocess, pysam, os

class IsoAln(object):
    def __init__(self, name=None, q=None, p=None, cigar=None, tlen=None, als=None):
        self.name = name
        self.q = q
        self.startpos = p
        self.cigar = cigar
        self.tlen = tlen
        self.alignscore = als


def getannotinfo(isoformfile):
    ###rewrite for fusion isoforms
    transcripttoexons = {}
    for line in open(isoformfile):
        line = line.rstrip().split('\t')
        name, left, right, chrom, strand = line[3], int(line[1]), int(line[2]), line[0], line[5]
        name = name.split('__')[0]
        blocksizes = [int(n) for n in line[10].rstrip(',').split(',')]
        if strand == '-': blocksizes = blocksizes[::-1]
        if name not in transcripttoexons: transcripttoexons[name] = []
        transcripttoexons[name].extend(blocksizes)
    return transcripttoexons


def process_cigar(cigarblocks, startpos):
    matchpos = 0
    coveredpos = [0] * (startpos - 1)
    queryclipping = []
    tendpos = startpos - 1
    blockstarts, blocksizes = [], []
    for btype, blen in cigarblocks:
        if btype in {4, 5}:  # btype == 'S' or btype == 'H':
            queryclipping.append(blen)
        elif btype == 0:  # btype == 'M' and (args.stringent or args.check_splice):
            # coveredpos.extend([1] * blen)
            coveredpos.extend([1] * blen)  ##still not checking matchvals
            blockstarts.append(tendpos)
            blocksizes.append(blen)
            matchpos += blen
            tendpos += blen
        elif btype in {2, 3}:  # btype == 'D' or btype == 'N':
            coveredpos.extend([0] * blen)
            tendpos += blen
            if blen > large_indel_tolerance: return True, None, None, None, None, None  #
        elif btype == 1:  # == 'I':
            coveredpos[-1] += blen
            if blen > large_indel_tolerance: return True, None, None, None, None, None
    return False, coveredpos, queryclipping, blockstarts, blocksizes, tendpos


def check_singleexon(read_start, read_end, tlen):
    if read_start < 25 and read_end > tlen - 25:
        return True
    else:
        return False


def check_exonenddist(blocksize, disttoend, trust_ends, disttoblock):
    if blocksize < 25:  # if the terminal exon is sufficiently short, relax by 5 bp bc drna misses bases on 5' end
        return disttoend < 5
    elif trust_ends:
        return disttoend <= trust_ends_window
    else:
        return disttoblock > 25


def check_firstlastexon(first_blocksize, last_blocksize, read_start, read_end, tlen, trust_ends):
    left_coverage = check_exonenddist(first_blocksize, read_start, trust_ends, first_blocksize - read_start)
    right_coverage = check_exonenddist(last_blocksize, tlen - read_end, trust_ends, read_end - (tlen - last_blocksize))
    # print(first_blocksize-read_start, read_end - (tlen - last_blocksize), read_start, read_end)
    return right_coverage and left_coverage


def check_stringent(coveredpos, exonpos, tlen, blockstarts, blocksizes, trust_ends):
    matchpos = len([x for x in coveredpos if x == 1])
    ###I think that the 80% of the transcript rule is less important than the exists in both first and last exons. Could add a toggle for this though.
    # if matchpos/tlen < 0.8:
    # 	return False
    # else:
    read_start, read_end = blockstarts[0], blockstarts[-1] + blocksizes[-1]
    first_blocksize, last_blocksize = exonpos[0], exonpos[-1]
    # covers enough bases into the first and last exons
    if len(blocksizes) == 1:  # single exon transcript
        return check_singleexon(read_start, read_end, tlen)
    else:
        return check_firstlastexon(first_blocksize, last_blocksize, read_start, read_end, tlen, trust_ends)


def check_splicesites(coveredpos, exonpos, tstart, tend, tname):
    currpos = 0
    for i in range(len(exonpos) - 1):
        elen = exonpos[i]
        currpos += elen
        if tstart < currpos < tend:
            ssvals = coveredpos[currpos - 3:currpos + 3]
            totinsert = sum([x - 1 for x in ssvals if x > 1])  ###value is match = 1 + insertsize
            totmatch = sum([1 for x in ssvals if x >= 1])  ##insert at pos still counts as match
            if totmatch - totinsert <= num_match_in_ss_window:
                return False
    return True


def checktranscriptinannot(exondict, tname):
    try:
        exoninfo = exondict[tname]
    except KeyError:
        raise Exception(
            "The transcript names in the annotation fasta do not appear to match the ones in the isoforms file. You may be able to fix this by using gtf_to_bed and bed_to_sequence on your annotation gtf and using the resulting file as your annotation fasta input to this program")
    except Exception as ex:
        raise Exception("** check_splice FAILED for %s" % (tname)) from ex
    return exoninfo


def getbesttranscript(tinfo, transcripttoexons):
    ###parse CIGAR + MD tag to ID transcript pos covered by alignment
    ###get start + end of transcript on read, alignment block positions
    ###also save soft/hard clipping at ends of read
    ###not positions of insertions larger than min_insertion_len, apply those to check_splice
    ##filter out reads with long indels
    ###generate list of 0s and 1s - transcript pos with match to query, val > 1 = insertion
    passingtranscripts = []
    # print(tinfo)
    for tname in tinfo:
        thist = tinfo[tname]
        ###process MD tag here to query positions with mismatches
        ##for MD tag, keep track of position of mismatch in all match positions
        indel_detected, coveredpos, queryclipping, blockstarts, blocksizes, tendpos = process_cigar(thist.cigar,
                                                                                                    thist.startpos)
        if not indel_detected:
            if check_stringentandsplice(transcripttoexons, thist.name, coveredpos, thist.tlen, blockstarts,
                                        blocksizes, thist.startpos, tendpos, trust_ends=True):
                passingtranscripts.append(
                    [-1 * thist.alignscore, -1 * sum([1 for x in coveredpos if x == 1]), sum(queryclipping), thist.tlen,
                     tname])
    ###order passing transcripts by alignment score
    ###then order by amount of query covered
    ##then order by amount of transcript covered
    if len(passingtranscripts) > 0:
        passingtranscripts.sort()
        return passingtranscripts[0][-1]
    else:
        return None


def parsesam(sam, transcripttoexons):
    lastread = None
    curr_transcripts = {}
    transcripttoreads = {}
    samfile = pysam.AlignmentFile(sam, 'r')
    for read in samfile:
        if read.is_mapped:
            readname = read.query_name
            transcript = read.reference_name
            pos = read.reference_start
            alignscore = read.get_tag('AS')
            cigar = read.cigartuples
            # mdtag = read.get_tag('MD')
            ###not using MD tag, not caring about mismatches
            tlen = samfile.get_reference_length(transcript)
            quality= read.mapping_quality
            if lastread and readname != lastread:
                assignedt = getbesttranscript(curr_transcripts, transcripttoexons)
                if assignedt:
                    if assignedt not in transcripttoreads: transcripttoreads[assignedt] = []
                    transcripttoreads[assignedt].append(lastread)

                curr_transcripts = {}
            curr_transcripts[transcript] = IsoAln(transcript, quality, pos, cigar, tlen, alignscore)
            lastread = readname
    if lastread:
        assignedt = getbesttranscript(curr_transcripts, transcripttoexons)
        if assignedt:
            if assignedt not in transcripttoreads: transcripttoreads[assignedt] = []
            transcripttoreads[assignedt].append(lastread)

    return transcripttoreads


def check_stringentandsplice(transcripttoexons, tname, coveredpos, tlen, blockstarts, blocksizes, tstart, tend,
                             trust_ends):
    exoninfo = checktranscriptinannot(transcripttoexons, tname)
    passesstringent = check_stringent(coveredpos, exoninfo, tlen, blockstarts, blocksizes, trust_ends)
    passessplice = check_splicesites(coveredpos, exoninfo, tstart, tend, tname)
    # passessplice = True
    return passesstringent and passessplice

min_insertion_len = 3
# ss_window = 3 unused
num_match_in_ss_window = 5 #/6
trust_ends_window = 50
large_indel_tolerance = 25

xvals = [5, 10, 50, 100]
toolnames = ['jaffal', 'ctat', 'flair']
scorelabels = ['recall', 'precision', 'f1score']

for test in ['constbp', 'multbp']:#
    truthfile = '/private/groups/brookslab/cafelton/simulatedfusionsjaffal/fusionisos/truth/fusioniso' + test + '-Feb25.fa'
    truthbedfile = '/private/groups/brookslab/cafelton/simulatedfusionsjaffal/fusionisos/truth/fusioniso' + test + '-Feb25.bed'
    transcripttoexons = getannotinfo(truthbedfile)
    toolyvals = [[[], [], []] for x in range(len(toolnames))]
    for cov in ['5', '10', '50', '100']:
        jaffalfile = '/private/groups/brookslab/cafelton/simulatedfusionsjaffal/fusionisos/jaffal/jaffal_results/fusioniso' + test + '.0225.badread' + cov + 'x_jaffa_results.fasta'
        ctatfile = '/private/groups/brookslab/cafelton/simulatedfusionsjaffal/fusionisos/ctat/fusioniso' + test + '.0225.badread' + cov + 'x.ctat.fusion_transcripts.fa'
        flairfile = '/private/groups/brookslab/cafelton/simulatedfusionsjaffal/fusionisos/flairfusion/fusioniso' + test + '.0225.badread' + cov + 'x.fusions.isoforms.fa'
        toolfiles = [jaffalfile, ctatfile, flairfile]
        for t in range(len(toolfiles)):
            print(test, cov, toolnames[t])
            mm2_cmd = ['minimap2', '-a', '-N', '4', '-t', '12', truthfile, toolfiles[t]]
            samname = '/private/groups/brookslab/cafelton/simulatedfusionsjaffal/fusionisos/fusionisoeval/' + '.'.join([test, cov, toolnames[t], 'alignedtotruth', 'sam'])
            # if toolnames[t] == 'flair':
            #     subprocess.check_call(mm2_cmd, stdout=open(samname, 'w'))
            transcripttoreads = parsesam(samname, transcripttoexons)
            TP = len(set(transcripttoreads.keys()))
            FN = len(set(transcripttoexons.keys()) - set(transcripttoreads.keys()))
            TP2 = sum(len(x) for x in transcripttoreads.values())
            allfound = 0
            for line in open(toolfiles[t]):
                if line[0] == '>': allfound += 1
            FP = allfound - TP
            FP2 = allfound - TP2
            recall = round(TP/(TP+FN), 3) if TP+FN > 0 else 0
            precision = round(TP/(TP+FP), 3) if TP+FP > 0 else 0
            f1score = round((2*TP)/((2*TP)+FN+FP), 3) if TP+FN+FP > 0 else 0
            recall2 = round(TP2 / (TP2 + FN), 3) if TP2 + FN > 0 else 0
            precision2 = round(TP2 / (TP2 + FP2), 3) if TP2 + FP2 > 0 else 0
            f1score2 = round((2 * TP2) / ((2 * TP2) + FN + FP2), 3) if TP2 + FN + FP2 > 0 else 0

            scores = [recall, precision, f1score]
            # scores = [recall2, precision2, f1score2]
            for i in range(3):
                toolyvals[t][i].append(scores[i])

    fig, axs = plt.subplots(3)
    fig.set_figwidth(6)
    fig.set_figheight(18)
    for t in range(len(toolnames)):
        for i in range(3):
            yvals = toolyvals[t][i]
            axs[i].plot(xvals, yvals, label=toolnames[t])

    for i in range(3):
        axs[i].set_ylabel(scorelabels[i])
        axs[i].set_xscale('log')
        axs[i].set_ylim(0,1)
        axs[i].legend()
    # plt.savefig(test + '_isoforms_prf1-nochecksplice-generousprecision.png', dpi=600)
    # plt.savefig(test + '_isoforms_prf1-generousprecision.png', dpi=600)
    plt.savefig(test + '_isoforms_prf1-friafternoon.png', dpi=600)
    # plt.savefig(test + '_isoforms_prf1-friafternoon-nochecksplice-generousprecision.png', dpi=600)

