import os
from comparefusiondetection_Jan2025 import *
import matplotlib.pyplot as plt

def getTrueBreakpoints(bedfile):
    fnametobp = {}
    for line in open(bedfile):
        line = line.rstrip().split('\t')
        fname, pos = line[3].split('__')
        chrom, startpos, endpos, strand = line[0], int(line[1]), int(line[2]), line[5]
        if fname not in fnametobp: fnametobp[fname] = []
        if pos == '0' and strand == '+' or pos == '1' and strand == '-':  bppos = endpos
        else: bppos = startpos
        fnametobp[fname].append((chrom, bppos))
    fusionbreakpoints = set()
    for f in fnametobp:
        fusionbreakpoints.add(tuple(fnametobp[f]))
    return fusionbreakpoints



for test in ['multbp', 'nongenic']:
    truthfile = '/private/groups/brookslab/cafelton/simulatedfusionsjaffal/fusionisos/truth/fusioniso' + test + '-Feb25.bed'
    truebp = getTrueBreakpoints(truthfile)
    toolnames = ['longgf', 'jaffal', 'ctat', 'flair']
    fig, axs = plt.subplots(3)
    fig.set_figwidth(6)
    fig.set_figheight(18)

    xvals = [5, 10, 50, 100]
    toolyvals = [[[],[],[]] for x in range(len(toolnames))]

    for cov in ['5', '10', '50', '100']:
        longgffile = '/private/groups/brookslab/cafelton/simulatedfusionsjaffal/fusionisos/longgf/fusioniso' + test + '.0225.badread' + cov + 'x_fusions.txt'
        jaffalfile = '/private/groups/brookslab/cafelton/simulatedfusionsjaffal/fusionisos/jaffal/jaffal_results/fusioniso' + test + '.0225.badread' + cov + 'x_jaffa_results.csv'
        ctatfile = '/private/groups/brookslab/cafelton/simulatedfusionsjaffal/fusionisos/ctat/fusioniso' + test + '.0225.badread' + cov + 'x/ctat-LR-fusion.fusion_predictions.abridged.tsv'
        flairfile = '/private/groups/brookslab/cafelton/simulatedfusionsjaffal/fusionisos/flairfusion/fusioniso' + test + '.0225.badread' + cov + 'x.fusions.isoforms.bed'
        toolinfo = [getLonggfFusions(longgffile), getJaffalFusions(jaffalfile), getCtatFusions(ctatfile), getFlairFusions(flairfile)]
        for t in range(len(toolinfo)):
            fusions = toolinfo[t]
            TP, FN = calculateRecall(truebp, fusions, 2, False, False)
            FP = calculatePrecision(truebp, fusions, 2, False)
            recall = round(TP/(TP+FN), 3) if TP+FN > 0 else 0
            precision = round(TP/(TP+FP), 3) if TP+FP > 0 else 0
            f1score = round((2*TP)/((2*TP)+FN+FP), 3) if TP+FN+FP > 0 else 0
            scores = [recall, precision, f1score]
            for i in range(3):
                toolyvals[t][i].append(scores[i])
    scorelabels = ['recall', 'precision', 'f1score']

    for t in range(len(toolnames)):
        for i in range(3):
            yvals = toolyvals[t][i]
            axs[i].plot(xvals, yvals, label=toolnames[t])

    for i in range(3):
        axs[i].set_ylabel(scorelabels[i])
        axs[i].set_xscale('log')
        axs[i].set_ylim(0,1)
        axs[i].legend()
    plt.savefig(test + '_breakpoints_prf1-friafternoon.png', dpi=600)


