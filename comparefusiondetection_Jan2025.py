import os
import matplotlib.pyplot as plt

def getTrueFusions():
    fusioninfo = {'Pac':{'75':{}, '80':{}, '85':{}, '90':{}, '95':{}}, 'ONT':{'75':{}, '80':{}, '85':{}, '90':{}, '95':{}}}
    potentialcov = [1,5,10,50,100]
    for line in open('simfusionstocovandhg38bp_withchr_Jan2025.tsv'):
        line = line.rstrip().split('\t')
        fusion, paccovinfo, ontcovinfo, g1info, g2info = line
        paccovinfo = paccovinfo.split('|')
        ontcovinfo = ontcovinfo.split('|')
        g1info = g1info.split('|')
        g2info = g2info.split('|')
        fusion = fusion.split('--')
        fusiongenes = (fusion[0].split('|')[1].split('.')[0], fusion[1].split('|')[1].split('.')[0])
        paccov = int(paccovinfo[2])
        if paccov > 0:
            paccov = potentialcov[min(range(len(potentialcov)), key = lambda i: abs(potentialcov[i]-paccov))]
            if paccov > 1:
                fusioninfo[paccovinfo[0]][paccovinfo[1]][((g1info[2], int(g1info[3])), (g2info[2], int(g2info[3])))] = paccov
        ontcov = int(ontcovinfo[2])
        if ontcov > 0:
            ontcov = potentialcov[min(range(len(potentialcov)), key = lambda i: abs(potentialcov[i]-ontcov))]
            if ontcov > 1:
                fusioninfo[ontcovinfo[0]][ontcovinfo[1]][((g1info[2], int(g1info[3])), (g2info[2], int(g2info[3])))] = ontcov
    return fusioninfo

potentialcov = [5,10,50,100]

def getLonggfFusions(filename, only1bp=False):
    allfusions = {}
    for line in open(filename):
        if line[:5] == 'SumGF':
            line = line.rstrip().split()
            fname = line[1]
            readsup = int(line[2])
            chr1, pos1 = line[3].split(':')
            chr2, pos2 = line[4].split(':')
            if readsup >= 2:
                # if fname not in allfusions: allfusions[fname] = set()
                # allfusions[fname].add(((chr1, int(pos1)), (chr2, int(pos2))))
                info = ((chr1, int(pos1)), (chr2, int(pos2)))
                nearestcov = potentialcov[min(range(len(potentialcov)), key=lambda i: abs(potentialcov[i] - readsup))]
                if fname not in allfusions:
                    allfusions[fname] = {}
                    if only1bp:
                        if info not in allfusions[fname]:
                            allfusions[fname][info] = nearestcov
                        else:
                            allfusions[fname][info] += nearestcov
                if not only1bp:
                    if info not in allfusions[fname]:
                        allfusions[fname][info] = nearestcov
                    else:
                        allfusions[fname][info] += nearestcov
    return allfusions

def getJaffalFusions(filename, only1bp=False):
    allfusions = {}
    for line in open(filename):
        line = line.rstrip().split(',')
        if line[0] != 'sample':
            #sample,fusion genes,chrom1,base1,strand1,chrom2,base2,strand2,gap (kb),spanning pairs,spanning reads,inframe,aligns,rearrangement,contig,contig break,classification,known
            fname = line[1]
            readsup = int(line[10])
            chr1, pos1 = line[2], line[3]
            chr2, pos2 = line[5], line[6]
            if readsup >= 2:
                # if fname not in allfusions: allfusions[fname] = set()
                # allfusions[fname].add(((chr1, int(pos1)), (chr2, int(pos2))))
                info = ((chr1, int(pos1)), (chr2, int(pos2)))
                nearestcov = potentialcov[min(range(len(potentialcov)), key=lambda i: abs(potentialcov[i] - readsup))]
                if fname not in allfusions:
                    allfusions[fname] = {}
                    if only1bp:
                        if info not in allfusions[fname]:
                            allfusions[fname][info] = nearestcov
                        else:
                            allfusions[fname][info] += nearestcov
                if not only1bp:
                    if info not in allfusions[fname]:
                        allfusions[fname][info] = nearestcov
                    else:
                        allfusions[fname][info] += nearestcov
    return allfusions

def getCtatFusions(filename, only1bp=False):
    allfusions = {}
    for line in open(filename):
        if line[0] != '#':
            line = line.rstrip().split('\t')
            ##FusionName	num_LR	LeftGene	LeftLocalBreakpoint	LeftBreakpoint	RightGene	RightLocalBreakpoint	RightBreakpoint	SpliceType	LR_accessions	LR_FFPM	annots	max_LR_FFPM	frac_dom_iso	above_frac_dom_iso
            fname = line[0]
            readsup = int(line[1])
            chr1, pos1, strand1 = line[4].split(':')
            chr2, pos2, strand2 = line[7].split(':')
            if readsup >= 2:
                # if fname not in allfusions: allfusions[fname] = set()
                # allfusions[fname].add(((chr1, int(pos1)), (chr2, int(pos2))))
                info = ((chr1, int(pos1)), (chr2, int(pos2)))
                nearestcov = potentialcov[min(range(len(potentialcov)), key = lambda i: abs(potentialcov[i]-readsup))]
                if fname not in allfusions:
                    allfusions[fname] = {}
                    if only1bp:
                        if info not in allfusions[fname]: allfusions[fname][info] = nearestcov
                        else: allfusions[fname][info] += nearestcov
                if not only1bp:
                    if info not in allfusions[fname]: allfusions[fname][info] = nearestcov
                    else: allfusions[fname][info] += nearestcov
    return allfusions

def getFlairCoverage(prefix):
    c = 0
    ognametonewname = {}
    nametocount = {}
    for line in open(prefix + '.syntheticAligned.flair.isoforms.bed'):
        c += 1
        line = line.split('\t')
        ogname = line[3]
        newname = 'fusioniso' + str(c) + '_' + '_'.join(ogname.split('_')[1:]).split(':')[0]
        ognametonewname[ogname] = newname
    for line in open(prefix + '.syntheticAligned.flair.firstpass.q.counts'):
        name, count = line.rstrip().split('\t')
        if name in ognametonewname:
            nametocount[ognametonewname[name]] = int(count)
    return nametocount

def getFlairFusions(filename, getcov=False, only1bp=False):
    allfusions = {}
    lastname, lastinfo, lastfname, lastcov = None, [], None, None
    prefix = filename.split('.fusions')[0]
    foundcounts = getFlairCoverage(prefix) if getcov else None
    # print(foundcounts)
    for line in open(filename):
        line = line.rstrip().split('\t')
        chr, p1, p2 = line[0], int(line[1]), int(line[2]) #float(line[1]), float(line[2])#
        strand = line[5]
        tname = '_'.join(line[3].split('_')[1:])
        fname = line[3].split('_')[-1]
        # if 'chr' in fname: continue
        glabel = line[3].split('_')[0]
        readsup = foundcounts['_'.join(line[3].split('_')[1:])] if getcov else -1
        # nearestcov = potentialcov[min(range(len(potentialcov)), key=lambda i: abs(potentialcov[i] - readsup))]
        ###only parsing all as two gene fusions, not worrying about internal gene breakpoints
        # if int(glabel[-1]) > 2: continue
        bppos = None
        if glabel == 'gene1' and strand == '+': bppos = p2
        elif glabel == 'gene1' and strand == '-': bppos = p1
        elif glabel == 'gene2' and strand == '+': bppos = p1
        elif glabel == 'gene2' and strand == '-': bppos = p2
        if tname != lastname:
            if lastname and len(lastinfo) == 2:
                # if fname not in allfusions: allfusions[fname] = set()
                # allfusions[fname].add(tuple(lastinfo))
                lastinfo = tuple(lastinfo)
                ###ADD reference counts file to get counts for fusions!!
                if lastfname not in allfusions:
                    allfusions[lastfname] = {}
                    if only1bp:
                        if lastinfo not in allfusions[lastfname]: allfusions[lastfname][lastinfo] = lastcov
                        else: allfusions[lastfname][lastinfo] += lastcov
                if not only1bp:
                    if lastinfo not in allfusions[lastfname]:
                        allfusions[lastfname][lastinfo] = lastcov
                    else:
                        allfusions[lastfname][lastinfo] += lastcov
            lastname, lastfname, lastcov = tname, fname, readsup
            lastinfo = [(chr, bppos)]
        else: lastinfo.append((chr, bppos))
    if lastname and lastinfo == 2:
        lastinfo = tuple(lastinfo)
        ###ADD reference counts file to get counts for fusions!!
        # readsup = foundcounts['_'.join(line[3].split('_')[1:])]
        # nearestcov = potentialcov[min(range(len(potentialcov)), key=lambda i: abs(potentialcov[i] - readsup))]
        if lastfname not in allfusions:
            allfusions[lastfname] = {}
        # else:
            # print(lastfname, lastname, lastinfo, lastcov)
            if only1bp:
                if lastinfo not in allfusions[lastfname]:
                    allfusions[lastfname][lastinfo] = lastcov
                else:
                    allfusions[lastfname][lastinfo] += lastcov
        if not only1bp:
            if lastinfo not in allfusions[lastfname]:
                allfusions[lastfname][lastinfo] = lastcov
            else:
                allfusions[lastfname][lastinfo] += lastcov
    for f in allfusions:
        for i in allfusions[f]:
            readsup = allfusions[f][i]
            allfusions[f][i] = potentialcov[min(range(len(potentialcov)), key=lambda i: abs(potentialcov[i] - readsup))]

    return allfusions


def calculateRecall(truefusions, foundfusions, buffersize, returnmatches=False, returnrecall=True):
    tottrue, totmatch = 0, 0
    matchf = set()
    for g1, g2 in truefusions:
        g1chr, g1pos = g1
        g2chr, g2pos = g2
        fullfusionfound = False
        for f in foundfusions:
            for bpset in foundfusions[f]:
                f1chr, f1pos = bpset[0]
                f2chr, f2pos = bpset[1]
                if g1chr == f1chr and g1pos-buffersize <= f1pos <= g1pos+buffersize and \
                   g2chr == f2chr and g2pos - buffersize <= f2pos <= g2pos + buffersize:
                    fullfusionfound = True
                    break
            if fullfusionfound == True: break
        tottrue += 1
        if fullfusionfound:
            totmatch += 1
            matchf.add((g1, g2))
    if returnmatches: return matchf
    else:
        if returnrecall: return round(totmatch/tottrue, 3)
        else: return totmatch, tottrue-totmatch

def calculatePrecision(truefusions, foundfusions, buffersize, returnprecision=True):
    totfound, totmatch = 0, 0
    for f in foundfusions:
        for bpset in foundfusions[f]:
            f1chr, f1pos = bpset[0]
            f2chr, f2pos = bpset[1]
            totfound += 1
            for g1, g2 in truefusions:
                g1chr, g1pos = g1
                g2chr, g2pos = g2
                if g1chr == f1chr and g1pos - buffersize <= f1pos <= g1pos + buffersize and \
                        g2chr == f2chr and g2pos - buffersize <= f2pos <= g2pos + buffersize:
                    totmatch += 1
                    break
    if returnprecision: return round(totmatch / totfound, 3)
    else: return totfound-totmatch

def calculatePrecisionByCov(truefusions, foundfusions, buffersize):
    totfound, totmatch = {5:0,10:0,50:0,100:0}, {5:0,10:0,50:0,100:0}
    for f in foundfusions:
        for bpset in foundfusions[f]:
            isfound = False
            f1chr, f1pos = bpset[0]
            f2chr, f2pos = bpset[1]
            foundcov = foundfusions[f][bpset]
            for g1, g2 in truefusions:
                g1chr, g1pos = g1
                g2chr, g2pos = g2
                if g1chr == f1chr and g1pos - buffersize <= f1pos <= g1pos + buffersize and \
                        g2chr == f2chr and g2pos - buffersize <= f2pos <= g2pos + buffersize:
                    truecov = truefusions[(g1, g2)]
                    totfound[truecov] += 1
                    totmatch[truecov] += 1
                    isfound = True
                    break
            if not isfound: totfound[foundcov] += 1
    return {c:totfound[c]-totmatch[c] for c in [5,10,50,100]}
    # return {c:round(totmatch[c]/totfound[c], 3) if totfound[c] > 0 else 0 for c in [5,10,50,100]}

def calculateRecallByCov(truefusions, foundfusions, buffersize, returnmatches=False):
    tottrue, totmatch = {5:0,10:0,50:0,100:0}, {5:0,10:0,50:0,100:0}
    matchf = set()
    for g1, g2 in truefusions:
        g1chr, g1pos = g1
        g2chr, g2pos = g2
        fullfusionfound = False
        truecov = truefusions[(g1, g2)]
        for f in foundfusions:
            for bpset in foundfusions[f]:
                f1chr, f1pos = bpset[0]
                f2chr, f2pos = bpset[1]
                if g1chr == f1chr and g1pos-buffersize <= f1pos <= g1pos+buffersize and \
                   g2chr == f2chr and g2pos - buffersize <= f2pos <= g2pos + buffersize:
                    fullfusionfound = True
                    break
            if fullfusionfound == True: break
        tottrue[truecov] += 1
        if fullfusionfound:
            totmatch[truecov] += 1
            matchf.add((g1, g2))
    if returnmatches: return matchf
    else:
        return {c:(totmatch[c], tottrue[c]-totmatch[c]) for c in [5,10,50,100]}
        # return {c:round(totmatch[c]/tottrue[c], 3) for c in [5,10,50,100]}


def calcPR(truefusions, foundfusions, buffersize):
    return calculatePrecision(truefusions, foundfusions, buffersize), calculateRecall(truefusions, foundfusions, buffersize)

# for buffersize in [0, 1, 2, 10, 100, 1000, 10000]:
#     totsim, totfound = 0, 0
#     for g1, g2 in fusioninfo['Pac']['95']:
#         if fusioninfo['Pac']['95'][g1, g2] == 1: continue
#         g1chr, g1pos = g1
#         g2chr, g2pos = g2
#         fullfusionfound = False
#         for f in allfusions:
#             # print(f, allfusions[f])
#             f1info, f2info = allfusions[f]['gene1'], allfusions[f]['gene2']
#             # if (g1chr == f1info[0] and f1info[1]-1000 <= g1pos <= f1info[2]+1000 and g2chr == f2info[0] and f2info[1]-1000 <= g2pos <= f2info[2]+1000) or \
#             #         (g2chr == f1info[0] and f1info[1]-1000 <= g2pos <= f1info[2]+1000 and g1chr == f2info[0] and f2info[1]-1000 <= g1pos <= f2info[2]+1000):
#             #     fullfusionfound = True
#             #     break
#             if (g1chr == f1info[0] and f1info[1] - buffersize <= g1pos <= f1info[2] + buffersize and g2chr == f2info[0] and f2info[1] - buffersize <= g2pos <= f2info[2] + buffersize):
#                 if (f1info[1]-buffersize <= g1pos <= f1info[1] + buffersize or f1info[2]-buffersize <= g1pos <= f1info[2] + buffersize) and (f2info[1]-buffersize <= g2pos <= f2info[1] + buffersize or f2info[2]-buffersize <= g2pos <= f2info[2] + buffersize):
#                     fullfusionfound = True
#                     break
#             elif (g2chr == f1info[0] and f1info[1]-buffersize <= g2pos <= f1info[2]+buffersize and g1chr == f2info[0] and f2info[1]-buffersize <= g1pos <= f2info[2]+buffersize):
#                 if (f1info[1]-buffersize <= g2pos <= f1info[1] + buffersize or f1info[2]-buffersize <= g2pos <= f1info[2] + buffersize) and (f2info[1]-buffersize <= g1pos <= f2info[1] + buffersize or f2info[2]-buffersize <= g1pos <= f2info[2] + buffersize):
#                     fullfusionfound = True
#                     break
#         totsim += 1
#         if fullfusionfound: totfound += 1
#         # else: print(g1, g2, fusioninfo['Pac']['95'][g1, g2])
#     print(buffersize, 'recall', totsim, totfound, totfound/totsim)
#     print(buffersize, 'precision', len(allfusions), totfound/len(allfusions))




# ###NOTE when checking based on coverage, only do it for 95err ont and pac for each
# truefusions = getTrueFusions()
# precisiondata, recalldata = [], []
# buffersize = 2
# # tnames = ['Pac-longgf', 'Pac-jaffal', 'Pac-ctat', 'ONT-longgf', 'ONT-jaffal', 'ONT-ctat']#, 'flair']
# # tnames = ['Pac-longgf', 'Pac-jaffal', 'Pac-ctat', 'Pac-FLAIR']#, 'ONT-longgf', 'ONT-jaffal', 'ONT-ctat', 'ONT-FLAIR']#, 'flair']
# # tnames = ['longgf-75', 'jaffal-75', 'ctat-75', 'flair-75', 'longgf-80', 'jaffal-80', 'ctat-80', 'flair-80', 'longgf-85', 'jaffal-85', 'ctat-85', 'flair-85', 'longgf-90', 'jaffal-90', 'ctat-90', 'flair-90', 'longgf-95', 'jaffal-95', 'ctat-95', 'flair-95']
# names = ['longgf', 'jaffal', 'ctat', 'flair']
# # error = '95'
# orecalls, oprecision, op1bp, of1scores, of11bp, onames = [], [], [], [], [], []
# for seqtype in ['Pac', 'ONT']:
#     for error in ['75', '80', '85', '90', '95']:
#         longgf = getLonggfFusions('1224.longgf/' + seqtype + '_' + error + '_NA12878chim10pct_lrgasp_fusions.txt')
#         jaffal = getJaffalFusions(
#             '1224.jaffal/jaffal_results/' + seqtype + '_' + error + 'err_NA12878chim10pct_jaffa_results.csv')
#         ctat = getCtatFusions(
#             '1224.ctat/' + seqtype + '_' + error + 'err_NA12878chim10pct/ctat-LR-fusion.fusion_predictions.tsv')
#         flair = getFlairFusions(seqtype + '_fus_sim_' + error + 'err_NA12878chim10pct.fusions.isoforms.bed', getcov=True)
#         longgf2 = getLonggfFusions('1224.longgf/' + seqtype + '_' + error + '_NA12878chim10pct_lrgasp_fusions.txt', only1bp=True)
#         jaffal2 = getJaffalFusions(
#             '1224.jaffal/jaffal_results/' + seqtype + '_' + error + 'err_NA12878chim10pct_jaffa_results.csv', only1bp=True)
#         ctat2 = getCtatFusions(
#             '1224.ctat/' + seqtype + '_' + error + 'err_NA12878chim10pct/ctat-LR-fusion.fusion_predictions.tsv', only1bp=True)
#         flair2 = getFlairFusions(seqtype + '_fus_sim_' + error + 'err_NA12878chim10pct.fusions.isoforms.bed', only1bp=True,
#                                 getcov=True)
#         # for i in range(len(toolres)):
#         #     m = toolres[i]
#         #     precisiondata.append(calculatePrecisionByCov(truefusions[seqtype][error], m, buffersize))
#         #     recalldata.append(calculateRecallByCov(truefusions[seqtype][error], m, buffersize))
#         toolinfo = [[longgf, longgf2], [jaffal, jaffal2], [ctat, ctat2], [flair,flair2]]
#         for i in range(len(toolinfo)):
#             t1, t2 = toolinfo[i]
#             recalls, precisions, p1bp, f1scores, f11bp = [], [], [], [], []
#             for cov in potentialcov:
#                 TP, FN = calculateRecallByCov(truefusions[seqtype][error], t1, buffersize)[cov]
#                 FP = calculatePrecisionByCov(truefusions[seqtype][error], t1, buffersize)[cov]
#                 recalls.append(round(TP / (TP + FN), 3) if TP + FN > 0 else 0)
#                 precisions.append(round(TP / (TP + FP), 3) if TP + FP > 0 else 0)
#                 f1scores.append(round((2 * TP) / ((2 * TP) + FN + FP), 3) if ((2 * TP) + FN + FP) > 0 else 0)
#                 TP, FN = calculateRecallByCov(truefusions[seqtype][error], t2, buffersize)[cov]
#                 FP = calculatePrecisionByCov(truefusions[seqtype][error], t2, buffersize)[cov]
#                 p1bp.append(round(TP / (TP + FP), 3) if TP + FP > 0 else 0)
#                 f11bp.append(round((2 * TP) / ((2 * TP) + FN + FP), 3) if ((2 * TP) + FN + FP) > 0 else 0)
#             orecalls.append(recalls)
#             oprecision.append(precisions)
#             op1bp.append(p1bp)
#             of1scores.append(f1scores)
#             of11bp.append(f11bp)
#             onames.append(seqtype + '-' + error + '-' + names[i])
#
#         # for cov in potentialcov:
#             # vals = [str(cov) + '-' + error]
#             # for i in range(len(toolres)):
#             #     m = toolres[i]
#             #     # vals.append(calculateRecallByCov(truefusions[seqtype][error], m, buffersize)[cov])
#             #     # vals.append(calculatePrecisionByCov(truefusions[seqtype][error], m, buffersize)[cov])
#             #     TP, FN = calculateRecallByCov(truefusions[seqtype][error], m, buffersize)[cov]
#             #     FP = calculatePrecisionByCov(truefusions[seqtype][error], m, buffersize)[cov]
#             #     vals.append(round((2*TP)/((2*TP)+FN+FP), 3) if ((2*TP)+FN+FP) > 0 else 0)
#             # print(' '.join([str(x) for x in vals]))
#
# metricnames = ['recall', 'precision', 'precision-1bp', 'f1scores', 'f1scores1bp']
# metrics = [orecalls, oprecision, op1bp, of1scores, of11bp]
# xvals = potentialcov#[75, 80, 85, 90, 95]
# tooltocolor = {'longgf':'blue', 'jaffal':'red', 'ctat':'orange', 'flair':'green'}
# seqtolinestyle = {'Pac':'solid', 'ONT':'dashed'}
# # errortoalpha = {'75':0.2, '80':0.4, '85':0.6, '90':0.8, '95':1}
# errortoindex = {'75':0, '80':1, '85':2, '90':3, '95':4}
# indextoerror = {0:'75', 1:'80', 2:'85', 3:'90', 4:'95'}
#
# fig, axs = plt.subplots(5, 5)
# fig.set_figwidth(20)
# fig.set_figheight(20)
# for i in range(5):
#     for j in range(len(onames)):
#         plotindex = errortoindex[onames[j].split('-')[1]]
#         axs[plotindex, i].plot(xvals, metrics[i][j], label=onames[j].split('-')[-1] + '-' + onames[j].split('-')[0], color=tooltocolor[onames[j].split('-')[-1]], linestyle=seqtolinestyle[onames[j].split('-')[0]])#, alpha=errortoalpha[onames[j].split('-')[1]])
#     for k in range(5):
#         axs[k, i].set_ylabel(indextoerror[k] + ':' + metricnames[i])
#         axs[k, i].set_xlabel('coverage')
#         axs[k, i].set_ylim(0, 1)
#         axs[k, i].set_xscale('log')
#     # axs[i].set_xscale('log')
#     # axs[i].legend()
# axs[-1, -1].legend()
# plt.savefig('bycov_precisionrecallf1plots_020625.png', dpi=600)
#     # print('\n' + metricnames[i] + '\n')
#     # print('error', ' '.join(tnames))
#     # for j in metrics[i]:
#     #     print(' '.join([str(x) for x in j]))


# print('precision\n')
# print('coverage', ' '.join(tnames))
# for cov in potentialcov:
#     vals = [cov]
#     for m in precisiondata:
#         vals.append(m[cov])
#     print(' '.join([str(x) for x in vals]))
# print('\nrecall\n')
# print('coverage', ' '.join(tnames))
# for cov in potentialcov:
#     vals = [cov]
#     for m in recalldata:
#         vals.append(m[cov])
#     print(' '.join([str(x) for x in vals]))




# truefusions = getTrueFusions()
# buffersize = 2
# # tnames = ['longgf', 'jaffal', 'ctat']#, 'flair']
# # tnames = ['Pac-longgf', 'Pac-jaffal', 'Pac-ctat', 'ONT-longgf', 'ONT-jaffal', 'ONT-ctat']#, 'flair']
# tnames = ['Pac-longgf', 'Pac-jaffal', 'Pac-ctat', 'Pac-flair', 'ONT-longgf', 'ONT-jaffal', 'ONT-ctat', 'ONT-flair']#, 'flair']
#
# # print('F1 score\n')
# # print('error', ' '.join(tnames))
# orecalls, oprecision, op1bp, of1scores, of11bp = [], [], [], [], []
# for error in ['75', '80', '85', '90', '95']:
#     # vals = [int(error)]
#     # recalls, precisions, p1bp, f1scores, f11bp = [int(error)], [int(error)], [int(error)], [int(error)], [int(error)]
#     recalls, precisions, p1bp, f1scores, f11bp = [], [], [], [], []
#     for seqtype in ['Pac', 'ONT']:
#         longgf = getLonggfFusions('1224.longgf/' + seqtype + '_' + error + '_NA12878chim10pct_lrgasp_fusions.txt')
#         jaffal = getJaffalFusions('1224.jaffal/jaffal_results/' + seqtype + '_' + error + 'err_NA12878chim10pct_jaffa_results.csv')
#         ctat = getCtatFusions('1224.ctat/' + seqtype + '_' + error + 'err_NA12878chim10pct/ctat-LR-fusion.fusion_predictions.tsv')
#         flair = getFlairFusions(seqtype + '_fus_sim_' + error + 'err_NA12878chim10pct.fusions.isoforms.bed')
#         longgf2 = getLonggfFusions('1224.longgf/' + seqtype + '_' + error + '_NA12878chim10pct_lrgasp_fusions.txt', only1bp=True)
#         jaffal2 = getJaffalFusions(
#             '1224.jaffal/jaffal_results/' + seqtype + '_' + error + 'err_NA12878chim10pct_jaffa_results.csv', only1bp=True)
#         ctat2 = getCtatFusions(
#             '1224.ctat/' + seqtype + '_' + error + 'err_NA12878chim10pct/ctat-LR-fusion.fusion_predictions.tsv', only1bp=True)
#         flair2 = getFlairFusions(seqtype + '_fus_sim_' + error + 'err_NA12878chim10pct.fusions.isoforms.bed', only1bp=True)
#         for m in [longgf, jaffal, ctat, flair]:
#             TP, FN = calculateRecall(truefusions[seqtype][error], m, buffersize, returnrecall=False)
#             FP = calculatePrecision(truefusions[seqtype][error], m, buffersize, returnprecision=False)
#             recalls.append(round(TP/(TP+FN), 3) if TP+FN > 0 else 0)
#             precisions.append(round(TP / (TP + FP), 3) if TP + FP > 0 else 0)
#             f1scores.append(round((2*TP)/((2*TP)+FN+FP), 3) if ((2*TP)+FN+FP) > 0 else 0)
#         for m in [longgf2, jaffal2, ctat2, flair2]:
#             TP, FN = calculateRecall(truefusions[seqtype][error], m, buffersize, returnrecall=False)
#             FP = calculatePrecision(truefusions[seqtype][error], m, buffersize, returnprecision=False)
#             p1bp.append(round(TP / (TP + FP), 3) if TP + FP > 0 else 0)
#             f11bp.append(round((2*TP)/((2*TP)+FN+FP), 3) if ((2*TP)+FN+FP) > 0 else 0)
#     orecalls.append(recalls)
#     oprecision.append(precisions)
#     op1bp.append(p1bp)
#     of1scores.append(f1scores)
#     of11bp.append(f11bp)
#
# metricnames = ['recall', 'precision', 'precision-1bp', 'f1scores', 'f1scores1bp']
# metrics = [orecalls, oprecision, op1bp, of1scores, of11bp]
# xvals = [75, 80, 85, 90, 95]
# tooltocolor = {'longgf':'blue', 'jaffal':'red', 'ctat':'orange', 'flair':'green'}
# seqtolinestyle = {'Pac':'solid', 'ONT':'dashed'}
#
# fig, axs = plt.subplots(1, 5)
# fig.set_figwidth(25)
# fig.set_figheight(6)
# for i in range(5):
#     for j in range(len(tnames)):
#         axs[i].plot(xvals, [x[j] for x in metrics[i]], label=tnames[j], color=tooltocolor[tnames[j][4:]], linestyle=seqtolinestyle[tnames[j][:3]])
#     axs[i].set_ylabel(metricnames[i])
#     # axs[i].set_xscale('log')
#     axs[i].set_ylim(0, 1)
#     # axs[i].legend()
# axs[-1].legend()
# plt.savefig('precisionrecallf1plots_020625.png', dpi=600)
#     # print('\n' + metricnames[i] + '\n')
#     # print('error', ' '.join(tnames))
#     # for j in metrics[i]:
#     #     print(' '.join([str(x) for x in j]))




# print('precision\n')
# print('error', ' '.join(tnames))
# for error in ['75', '80', '85', '90', '95']:
#     vals = [int(error)]
#     for seqtype in ['Pac', 'ONT']:
#         longgf = getLonggfFusions('1224.longgf/' + seqtype + '_' + error + '_NA12878chim10pct_lrgasp_fusions.txt')
#         jaffal = getJaffalFusions('1224.jaffal/jaffal_results/' + seqtype + '_' + error + 'err_NA12878chim10pct_jaffa_results.csv')
#         ctat = getCtatFusions('1224.ctat/' + seqtype + '_' + error + 'err_NA12878chim10pct/ctat-LR-fusion.fusion_predictions.tsv')
#         flair = getFlairFusions(seqtype + '_fus_sim_' + error + 'err_NA12878chim10pct.fusions.isoforms.bed')
#         for m in [longgf, jaffal, ctat, flair]:
#             vals.append(calculatePrecision(truefusions[seqtype][error], m, buffersize))
#     print(' '.join([str(x) for x in vals]))
#
# print('\nrecall\n')
# print('error', ' '.join(tnames))
# for error in ['75', '80', '85', '90', '95']:
#     vals = [int(error)]
#     for seqtype in ['Pac', 'ONT']:
#         longgf = getLonggfFusions('1224.longgf/' + seqtype + '_' + error + '_NA12878chim10pct_lrgasp_fusions.txt')
#         jaffal = getJaffalFusions('1224.jaffal/jaffal_results/' + seqtype + '_' + error + 'err_NA12878chim10pct_jaffa_results.csv')
#         ctat = getCtatFusions('1224.ctat/' + seqtype + '_' + error + 'err_NA12878chim10pct/ctat-LR-fusion.fusion_predictions.tsv')
#         flair = getFlairFusions(seqtype + '_fus_sim_' + error + 'err_NA12878chim10pct.fusions.isoforms.bed')
#         for m in [longgf, jaffal, ctat, flair]:
#             vals.append(calculateRecall(truefusions[seqtype][error], m, buffersize))
#     print(' '.join([str(x) for x in vals]))








# # longgf = getLonggfFusions('1224.longgf/Pac_95_NA12878chim10pct_lrgasp_fusions.txt')
# # jaffal = getJaffalFusions('1224.jaffal/jaffal_results/Pac_95err_NA12878chim10pct_jaffa_results.csv')
# # ctat = getCtatFusions('1224.ctat/Pac_95err_NA12878chim10pct/ctat-LR-fusion.fusion_predictions.tsv')
# # flair1 = getFlairFusions('Pac_fus_sim_95err_NA12878chim10pct.fusions.isoforms.standard.bed')
# # flair2 = getFlairFusions('Pac_fus_sim_95err_NA12878chim10pct.fusions.isoforms.multipleparalogs.bed')
# # flair3 = getFlairFusions('Pac_fus_sim_95err_NA12878chim10pct.fusions.isoforms.nocsnostringentmultpara.bed')
# # flair4 = getFlairFusions('Pac_fus_sim_95err_NA12878chim10pct.fusions.isoforms.onlycheckbpmultpara.bed')
# # flair5 = getFlairFusions('Pac_fus_sim_95err_NA12878chim10pct.fusions.isoforms.onlycheckbpmultparasscorr8.bed')
# flair6 = getFlairFusions('Pac_fus_sim_95err_NA12878chim10pct.fusions.isoforms.sunday020225.bed')
# # flair7 = getFlairFusions('Pac_fus_sim_95err_NA12878chim10pct.fusions.isoforms.bed')#, only1bp=True)
#
# # for f in flair6:
# #     print(f, flair6[f])
#
# # otherflair = getFlairFusions('Pac_fus_sim_95err_NA12878chim10pct.fusions.isoforms.multipleparalogs.bed')
# truefusions = getTrueFusions()
# tnames = ['flair']
# tfusions = [flair6]
# # tnames = ['longgf', 'jaffal', 'ctat', 'flair', 'flair-multiple-paralogs', 'flair-less-stringent', 'flair-only-check-bp', 'flair-only-check-bp-sscorr8', 'flair-removehighqdist', 'flair-ssaroundbponly']
# # tfusions = [longgf, jaffal, ctat, flair1, flair2, flair3, flair4, flair5, flair6, flair7]



# goodotherflair = calculateRecall(truefusions['Pac']['95'], flair3, 100, True)
# goodjaffal = calculateRecall(truefusions['Pac']['95'], jaffal, 100, True)
# goodflair = calculateRecall(truefusions['Pac']['95'], flair6, 100, True)
# goodctat = calculateRecall(truefusions['Pac']['95'], ctat, 100, True)
# # print('flair only')
# # for f in goodflair-goodjaffal:
# #     print(' '.join([':'.join([str(y) for y in x]) for x in f] ))
# # print('jaffal only')
# # for f in (goodjaffal|goodctat)-goodflair:
# #     print(' '.join([':'.join([str(y) for y in x]) for x in f] ))
# for f in goodotherflair-goodflair:
#     print(' '.join([':'.join([str(y) for y in x]) for x in f] ))
# # print('ctat only')
# # for f in (goodctat -goodjaffal)-goodflair:
# #     print(' '.join([':'.join([str(y) for y in x]) for x in f] ))
# # for f in goodotherflair-goodflair:
# #     print(' '.join([':'.join([str(y) for y in x]) for x in f] ))
# # print('f')
# # for f in (set(truefusions['Pac']['95'].keys())-goodflair):
# #     print(' '.join([':'.join([str(y) for y in x]) for x in f] ))
# # print('jaffal')
# # for f in (set(truefusions['Pac']['95'].keys())-goodjaffal):
# #     print(' '.join([':'.join([str(y) for y in x]) for x in f] ))



# print('precision\n')
# print('buffersize', ' '.join(tnames))
# for buffersize in [2]:#[0, 1, 2, 10, 100, 1000, 10000]:
#     vals = [buffersize]
#     for i in range(len(tfusions)):
#         vals.append(calculatePrecision(truefusions['Pac']['95'], tfusions[i], buffersize))
#     print(' '.join([str(x) for x in vals]))
# print('\nrecall\n')
# print('buffersize', ' '.join(tnames))
# for buffersize in [2]:#[0, 1, 2, 10, 100, 1000, 10000]:
#     vals = [buffersize]
#     for i in range(len(tfusions)):
#         vals.append(calculateRecall(truefusions['Pac']['95'], tfusions[i], buffersize))
#     print(' '.join([str(x) for x in vals]))


# # flairfiles = ['Pac_fus_sim_95err_NA12878chim10pct.fusions.isoforms.sunday020225.bed', 'Pac_fus_sim_95err_NA12878chim10pct.fusions.isoforms.monday020325.bed', 'Pac_fus_sim_95err_NA12878chim10pct.fusions.isoforms.bed']
# # names = ['flair-sunday', 'flair-monday', 'flair-now']
# # buffersize = 2
# # truefusions = getTrueFusions()
# # r, p1, p2, f1, f12 = [], [], [], [], []
# # for i in range(len(flaifiles)):
#     info = getFlairFusions(flairfiles[i])
#     info1bp = getFlairFusions(flairfiles[i], only1bp=True)
#     r.append(calculateRecall(truefusions['Pac']['95'], info, 2))
#     p1.append(calculatePrecision(truefusions['Pac']['95'], info, 2))
#     p2.append(calculatePrecision(truefusions['Pac']['95'], info1bp, 2))
#     TP, FN = calculateRecall(truefusions['Pac']['95'], info1bp, buffersize, returnrecall=False)
#     FP = calculatePrecision(truefusions['Pac']['95'], info1bp, buffersize, returnprecision=False)
#     f1.append(round((2 * TP) / ((2 * TP) + FN + FP), 3) if ((2 * TP) + FN + FP) > 0 else 0)
#     TP, FN = calculateRecall(truefusions['Pac']['95'], info, buffersize, returnrecall=False)
#     FP = calculatePrecision(truefusions['Pac']['95'], info, buffersize, returnprecision=False)
#     f12.append(round((2 * TP) / ((2 * TP) + FN + FP), 3) if ((2 * TP) + FN + FP) > 0 else 0)
# #
# # print('metric', ' '.join(names))
# # print('recall', ' '.join([str(x) for x in r]))
# # print('precision', ' '.join([str(x) for x in p1]))
# # print('precision-1breakpoint', ' '.join([str(x) for x in p2]))
# # print('f1score', ' '.join([str(x) for x in f12]))
# # print('f1score-1breakpoint', ' '.join([str(x) for x in f1]))


