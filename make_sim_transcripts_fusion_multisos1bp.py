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
	leftpos, rightpos = [], []
	for t in list(genetoinfo[gene]):
		if t not in {'chrom', 'strand'}:
			###this removes short exons and long introns!!!!
			genetoinfo[gene][t] = sorted(genetoinfo[gene][t])
			if len(genetoinfo[gene][t]) < 3: #or any(x[1] - x[0] <= 40 for x in genetoinfo[gene][t]) or any(genetoinfo[gene][t][x + 1][0] - genetoinfo[gene][t][x][1] >= 200000 for x in range(len(genetoinfo[gene][t]) - 1)):
				genetoinfo[gene].pop(t)
			else:
				exoncounts.append(len(genetoinfo[gene][t]))
				leftpos.append(genetoinfo[gene][t][0][0])
				rightpos.append(genetoinfo[gene][t][-1][1])

	genetoinfo[gene]['ecounts'] = exoncounts
	if len(exoncounts) >= 2:
		genetoinfo[gene]['coords'] = (min(leftpos), max(rightpos))
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
# fullisobed = open(outfilename + '.fulliso.bed', 'w')

# totmutcount = 2 ##must be even number for fusions
numisos = 2
fusioncount = 100
# genecount = 500

allgenes = list(genetoinfo.keys())
c = 0
while c < fusioncount:
	fgenes = []
	isonames = []
	genebreakpoints = []
	isoexons = []
	while len(fgenes) < 2:
		g = random.choice(allgenes)
		allgenes.remove(g)
		if sum([1 if x > 2 else 0 for x in genetoinfo[g]['ecounts']]) >= 2:
			genecoords = genetoinfo[g]['coords']
			z = 0
			goodisos = []
			while len(goodisos) < numisos and z < 5:
				tempbp = random.randint(genecoords[0], genecoords[1])
				potisos = list(set(genetoinfo[g]) - {'chrom', 'strand', 'ecounts', 'coords'})
				goodisos = []
				goodexons = []
				for iso in potisos:
					theseexons = genetoinfo[g][iso]
					thisstrand = genetoinfo[g]['strand']
					if theseexons[0][1] < tempbp < theseexons[-1][0]:
						goodisos.append(iso)
						croppedexons = []
						for e in theseexons:
							if (len(fgenes) == 0 and thisstrand == '+') or (len(fgenes) == 1 and thisstrand == '-'):
								if e[1] < tempbp: croppedexons.append(e)
							else:
								if e[0] > tempbp: croppedexons.append(e)
						if thisstrand == '-': croppedexons = [x[::-1] for x in croppedexons][::-1] ###reverse internally as well #croppedexons.sort(reverse=True)  # = theseexons[::-1]
						goodexons.append(croppedexons)
			if len(goodisos) >= numisos:
				for p in range(len(goodisos) - 1):
					###get internal exon boundaries
					intbound = []
					for i in range(2):
						b = []
						for j in range(len(goodexons[p+i])):
							if j == 0 and j != len(goodexons[p+i]) - 1: b.append(goodexons[p+i][j][1])
							elif j != 0 and j == len(goodexons[p+i]) - 1: b.append(goodexons[p+i][j][0])
							elif j != 0 and j != len(goodexons[p+i]) - 1:
								b.append(goodexons[p + i][j][0])
								b.append(goodexons[p + i][j][1])
						intbound.append(b)
					if ((len(fgenes) == 0 and goodexons[p][-1][1] == goodexons[p+1][-1][1]) or (len(fgenes) == 1 and goodexons[p][0][0] == goodexons[p+1][0][0])) and intbound[0] != intbound[1] and intbound[0] != [] and intbound[1] != []: ##same breakpoint, diff internal exons
						fgenes.append(g)
						isonames.append(goodisos[p:p+2])
						genebreakpoints.append(tempbp)
						isoexons.append(goodexons[p:p+2])
						break

	chrom, strand = [genetoinfo[fgenes[x]]['chrom'] for x in range(2)], [genetoinfo[fgenes[x]]['strand'] for x in range(2)]

	for fusiso in range(numisos):
		isoseq, isoname, isoinfo, isomuts = [], [], [], []
		for gindex in range(2):
			thisiso = isonames[gindex][fusiso]
			# theseexons = genetoinfo[fgenes[gindex]][thisiso]

			###reorder exons for negative strand genes - already sorted in forward direction above
			# if strand[gindex] == '+': theseexons.sort()
			##select exons for fusion
			# croppedexons = []
			# for e in theseexons:
			# 	if (gindex == 0 and strand[gindex] == '+') or (gindex==1 and strand[gindex] == '-'):
			# 		if e[1] < genebreakpoints[gindex]: croppedexons.append(e)
			# 	else:
			# 		if e[0] > genebreakpoints[gindex]: croppedexons.append(e)
			# numexons = len(croppedexons)
			# if strand[gindex] == '-': croppedexons.sort(reverse=True) #= theseexons[::-1]

			croppedexons = isoexons[gindex][fusiso]
			numexons = len(croppedexons)
			if strand[gindex] == '-': croppedexons = [x[::-1] for x in croppedexons]
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










