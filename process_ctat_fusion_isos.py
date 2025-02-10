from statistics import median
import pysam


for test in ['constbp', 'multbp', 'nongenic']:
    for cov in ['5', '10', '50', '100']:
        synthgenomefile = '/private/groups/brookslab/cafelton/simulatedfusionsjaffal/fusionisos/ctat/fusioniso' + test + '.0225.badread' + cov + 'x/fusion_intermediates_dir/LR-FI_targets.fa'
        fusiontranscriptsfile = '/private/groups/brookslab/cafelton/simulatedfusionsjaffal/fusionisos/ctat/fusioniso' + test + '.0225.badread' + cov + 'x/fusion_intermediates_dir/LR-FI.mm2.fusion_transcripts.gff3'
        outfile = '/private/groups/brookslab/cafelton/simulatedfusionsjaffal/fusionisos/ctat/fusioniso' + test + '.0225.badread' + cov + 'x.ctat.fusion_transcripts.fa'

        fusiontoisotoexons = {}
        for line in open(fusiontranscriptsfile):
            line = line.rstrip().split('\t')
            fname = line[0]
            start, end = int(line[3]), int(line[4])
            strand = line[6]
            tname = line[8].split(';')[0]
            if fname not in fusiontoisotoexons: fusiontoisotoexons[fname] = {}
            if tname not in fusiontoisotoexons[fname]: fusiontoisotoexons[fname][tname] = []
            fusiontoisotoexons[fname][tname].append((start, end))
        # print(len(fusiontoisotoexons['ANKEF1--ETNK1']))

        fusiontosjchains = {}
        for f in fusiontoisotoexons:
            if f not in fusiontosjchains: fusiontosjchains[f] = {}
            for iso in fusiontoisotoexons[f]:
                exons = sorted(fusiontoisotoexons[f][iso])
                ends = (exons[0][0], exons[-1][1])
                introns = []
                for i in range(len(exons) - 1):
                    introns.append((exons[i][1], exons[i+1][0]))
                introns = tuple(introns)
                if introns not in fusiontosjchains[f]: fusiontosjchains[f][introns] = []
                fusiontosjchains[f][introns].append(ends)
        # print(len(fusiontosjchains['ANKEF1--ETNK1']))
        # x = 0
        # for i in fusiontosjchains['ANKEF1--ETNK1']:
        #     x += len(fusiontosjchains['ANKEF1--ETNK1'][i])
            # print(len(fusiontosjchains['ANKEF1--ETNK1'][i]))
        # print(x)

        finalexonchains = []
        for f in fusiontosjchains:
            for ichain in fusiontosjchains[f]:
                ends = sorted(fusiontosjchains[f][ichain])
                groups, g = [], [(-200, -200)]
                for e in ends:
                    if e[0]-g[-1][0] > 100 and e[1]-g[-1][1] > 100:
                        if g[0] != (-200, -200):
                            groups.append(g)
                        g = [e]
                    else: g.append(e)
                groups.append(g)
                # if f == 'ANKEF1--ETNK1':
                #     print(f, len(groups), [len(g) for g in groups], ichain)
                centralexons = []
                for i in range(len(ichain) - 1):
                    centralexons.append((ichain[i][1], ichain[i+1][0]))
                for g in groups:
                    startpos = int(median([x[0] for x in g]))
                    endpos = int(median([x[1] for x in g]))
                    if len(ichain) > 0:
                        allexons = [(startpos, ichain[0][0])] + centralexons + [(ichain[-1][1], endpos)]
                    else: allexons = [(startpos, endpos)]
                    finalexonchains.append((f, allexons, len(g)))

        # genome = {}
        # last = None
        # for line in open(synthgenomefile):
        #     if line[0] == '>':
        #         last = line[1:].rstrip()
        #         genome[last] = []
        #     else: genome[last].append(line.rstrip())
        # for g in genome:
        #     genome[g] = ''.join(genome[g])
        genome = pysam.FastaFile(synthgenomefile)

        out = open(outfile, 'w')
        c = 0
        for f, exons, count in finalexonchains:
            seq = []
            c += 1
            name = f + '_iso' + str(c) + '_' + str(count) + 'reads'
            for e in exons:
                # seq.append(genome[f][e[0]:e[1]])
                seq.append(genome.fetch(f, e[0]-1, e[1]))
            seq = ''.join(seq)
            out.write('>' + name + '\n')
            out.write(seq + '\n')
        out.close()



