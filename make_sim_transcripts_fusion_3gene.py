import sys, csv, os, random
import seaborn as sns
import random
import pysam


try:
	fa = sys.argv[1]#open(sys.argv[1])  # transcript sequences
	gtf = sys.argv[2]
	outfilename = sys.argv[3]
except:
	# sys.stderr.write('usage: script.py transcripts.fa outputbase [bed]\n')
	sys.stderr.write('usage: script.py genome.fa gtfannot outputbase \n')
	sys.exit()

# def convertToGenomePos(tname, tpos, tchrom, tstrand, texons, oldbase, newbase):
# 	for exon in texons:
# 		if exon[1] - exon[0] > tpos:
# 			if tstrand == '+':
# 				return [tchrom, str(exon[0] + tpos - 1), '.', oldbase, newbase, tname]
# 			else:
# 				return [tchrom, str(exon[1] - tpos + 1), '.', oldbase, newbase, tname]
# 		else:
# 			tpos -= (exon[1] - exon[0] + 1)


# read in acceptable gene names
gene_names = set()

genetoinfo = {}
for line in open(gtf):
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	chrom, ty, start, end, strand = line[0], line[2], int(line[3])-1, int(line[4]), line[6]
	if ty != 'exon': continue
	if 'tag "CCDS"' not in line[8]: continue
	transcript_id = line[8].split('transcript_id "')[1].split('"')[0]
	gene_id = line[8].split('gene_id "')[1].split('"')[0]
	if gene_id not in genetoinfo: genetoinfo[gene_id] = {'chrom':chrom, 'strand':strand}
	if transcript_id not in genetoinfo[gene_id]: genetoinfo[gene_id][transcript_id] = []
	genetoinfo[gene_id][transcript_id].append((start, end))

for gene in genetoinfo:
	exoncounts = []
	for t in list(genetoinfo[gene]):
		if t not in {'chrom', 'strand'}:
			###this removes short exons and long introns!!!!
			genetoinfo[gene][t] = sorted(genetoinfo[gene][t])
			if len(genetoinfo[gene][t]) < 4: #or any(x[1] - x[0] <= 40 for x in genetoinfo[gene][t]) or any(genetoinfo[gene][t][x + 1][0] - genetoinfo[gene][t][x][1] >= 200000 for x in range(len(genetoinfo[gene][t]) - 1)):
				genetoinfo[gene].pop(t)
			else:
				exoncounts.append(len(genetoinfo[gene][t]))

	genetoinfo[gene]['ecounts'] = exoncounts
print('loaded annotation')
print(len(genetoinfo.keys()))

b = ['A', 'T', 'G', 'C']  # bases

compbase = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
def revcomp(seq):
	newseq = []
	for base in seq: newseq.append(compbase[base])
	return ''.join(newseq[::-1])

genome=pysam.FastaFile(fa)
print('loaded genome')
# print(yeastgenome.fetch('chrX', 74615, 74630))


seqout = open(outfilename + '.fa', 'w')
# vcfout = open(outfilename + '.vcf', 'w')
infoout = open(outfilename + '.bed', 'w')

# totmutcount = 2 ##must be even number for fusions
numisos = 1
fusioncount = 100
# genecount = 500
numfgenes = 3

allgenes = list(genetoinfo.keys())
c = 0
numsingle = 0
while c < fusioncount:
	fgenes = []
	while len(fgenes) < numfgenes:
		g = random.choice(allgenes)
		allgenes.remove(g)
		if sum([1 if x > 3 else 0 for x in genetoinfo[g]['ecounts']]) >= 2:
			fgenes.append(g)
	isonames = [list(set(genetoinfo[fgenes[x]]) - {'chrom', 'strand', 'ecounts'}) for x in range(numfgenes)]
	chrom, strand = [genetoinfo[fgenes[x]]['chrom'] for x in range(numfgenes)], [genetoinfo[fgenes[x]]['strand'] for x in range(numfgenes)]
	for fusiso in range(numisos):
		isoseq, isoname, isoinfo, isomuts = [], [], [], []
		for gindex in range(numfgenes):
			thisiso = None
			while not thisiso:
				tempiso = random.choice(isonames[gindex])
				isonames[gindex].remove(tempiso)
				if len(genetoinfo[fgenes[gindex]][tempiso]) > 2:
				# if len(genetoinfo[fgenes[gindex]][tempiso]) > 2 and all(x[1] - x[0] > 40 for x in genetoinfo[fgenes[gindex]][tempiso]) and all(genetoinfo[fgenes[gindex]][tempiso][x + 1][0] - genetoinfo[fgenes[gindex]][tempiso][x][1] < 200000 for x in range(len(genetoinfo[fgenes[gindex]][tempiso]) - 1)):
					thisiso = tempiso
			theseexons = genetoinfo[fgenes[gindex]][thisiso]


			###reorder exons for negative strand genes
			e1pos = None
			if strand[gindex] == '-': theseexons.sort(reverse=True) #= theseexons[::-1]
			##select exons for fusion
			numexons = random.randrange(1, len(theseexons))
			if gindex == 0: croppedexons = theseexons[:numexons]
			elif gindex == numfgenes - 1: croppedexons = theseexons[len(theseexons)-numexons:]
			else:
				e1pos = random.randrange(1, len(theseexons)-1)
				numexons = random.randrange(e1pos+1, len(theseexons)) - e1pos
				croppedexons = theseexons[e1pos:e1pos+numexons]
			# print(e1pos, numexons, len(croppedexons))
			if numexons == 1: numsingle += 1
			##get sequence for all exons
			exonseq = []
			for i in range(numexons):
				exonseq.append(genome.fetch(chrom[gindex], croppedexons[i][0], croppedexons[i][1]))

			##revcomp seq if required
			if strand[gindex] == '-':
				for i in range(numexons):
					exonseq[i] = revcomp(exonseq[i])
			isoseq.extend(exonseq)
			isoname.append(thisiso + '_' + fgenes[gindex])

			###build bed line
			croppedexons.sort()
			startpos, endpos = croppedexons[0][0], croppedexons[-1][-1]
			estarts = [str(x[0] - startpos) for x in croppedexons]
			esizes = [str(x[1] - x[0]) for x in croppedexons]
			# chr1	958245	959256	ENST00000469563.1_ENSG00000188976.11	1000	-	958245	959256	0	2	836,42,	0,969,
			bedline = [chrom[gindex], str(startpos), str(endpos), 'temp__' + str(gindex), '1000', strand[gindex], str(startpos), str(endpos), '0',
					   str(len(croppedexons)), ','.join(esizes) + ',', ','.join(estarts) + ',']
			isoinfo.append(bedline)
			# isoinfo.append(','.join([chrom, strand] + ['.'.join([str(y) for y in x]) for x in croppedexons]))
			#chr1	917495	.	C	T
		isoseq = ''.join(isoseq)
		isoname = '--'.join(isoname)
		seqout.write('>' + isoname + '\n')
		seqout.write(isoseq + '\n')
		# for vline in isomuts:
		# 	vcfout.write('\t'.join(vline + [isoname]) + '\n')
		for bedline in isoinfo:
			bedline[3] = isoname + '__' + bedline[3].split('__')[-1]
			infoout.write('\t'.join(bedline) + '\n')
		# infoout.write('\t'.join([isoname] + isoinfo) + '\n')
	c += 1


print(numsingle)





