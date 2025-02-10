
import subprocess, pipettor, os, pysam
reference_genome = '/private/groups/brookslab/reference_sequence/GRCh38.primary_assembly.genome.fa'
reference_annot = '/private/groups/brookslab/reference_annotations/gencode.v38.annotation.gtf'
flair_path = '/private/groups/brookslab/cafelton/git-flair/flair/'

#export PATH="/private/groups/brookslab/cafelton/git-flair/flair/bin:/private/groups/brookslab/cafelton/git-flair/flair/src/flair:$PATH"

c= 0
for simtype in ['nongenic', 'constbp', 'nongenic']:#['constbp', 'multbp', 'nongenic']:#['nongenic']:#['multbp']: #
    for cov in ['5', '10', '50', '100']:#['100']:#
        print(simtype, cov)
        fileprefix = 'fusioniso' + simtype + '.0225.badread' + cov + 'x'
        fastqfile = fileprefix + '.fastq'


        # ####NOTE need to modify flair align to not delete unfiltered file
        # aligncommand = 'python3 ' + flair_path + 'flair.py align -r ' + fastqfile + \
        #            ' -t 12 -g ' + reference_genome + ' --filtertype separate --minfragmentsize 40 --maxintronlen 350k -o ' + \
        #             fileprefix + '.hg38aligned'
        # subprocess.call(aligncommand, shell=True)
        #
        #
        # genomealignedfile = fileprefix + '.hg38aligned_unfiltered.bam'
        # if not os.path.isfile(fileprefix + '.hg38aligned_unfiltered.namesorted.bam'):
        #     subprocess.call(['samtools', 'sort', '-n', genomealignedfile], stdout=open(fileprefix + '.hg38aligned_unfiltered.namesorted.bam', 'w'))
        # print('done resorting alignment')
        # subprocess.call(['/private/groups/brookslab/cafelton/bin/LongGF/bin/LongGF', fileprefix + '.hg38aligned_unfiltered.namesorted.bam', reference_annot, '100', '50', '100'], stdout=open('longgf/' + fileprefix + '_fusions.txt', 'w'))
        # print('done longgf')
        #
        #
        # subprocess.call(['docker', 'run', '--rm', '-v', '/private/groups/brookslab/cafelton/simulatedfusionsjaffal/fusionisos:/data',
        #      '-v',
        #      '/private/groups/brookslab/cafelton/bin/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir:/ctat_genome_lib',
        #      'trinityctat/ctat_lr_fusion', 'ctat-LR-fusion', '-T', '/data/' + fastqfile, '--genome_lib_dir',
        #      '/ctat_genome_lib', '-o', '/data/ctat/' + fileprefix, '--vis'])
        # print('done ctat')
        #
        #
        # if not os.path.isfile(fastqfile + '.gz'):
        #     subprocess.call(['gzip', '-c', fastqfile], stdout = open(fastqfile + '.gz', 'w'))
        # print('done gzipping')
        # subprocess.call(['mv', fastqfile + '.gz', 'jaffal/'])
        # os.chdir('jaffal')
        # subprocess.call(['/private/groups/brookslab/cafelton/bin/JAFFA-version-2.3/tools/bin/bpipe', 'run', '/private/groups/brookslab/cafelton/bin/JAFFA-version-2.3/JAFFAL.groovy', fastqfile + '.gz'])
        # print('done jaffal')
        # subprocess.call(['mv', 'jaffa_results.csv', 'jaffal_results/' + fileprefix + '_jaffa_results.csv'])
        # subprocess.call(['mv', 'jaffa_results.fasta', 'jaffal_results/' + fileprefix + '_jaffa_results.fasta'])
        # os.chdir('..')

        ###filter bam reads to chimeras
        ##align to transcriptome
        ##filter to chimeras
        ###align to transcriptome with --secondary=no

        # os.chdir('flairfusion')

        # mm2_cmd = ('minimap2', '-a', '-s', '40', '-t', '12', '--secondary=no', '/private/groups/brookslab/cafelton/simulatedfusionsjaffal/gencode.v38.flair.fa', '../' + fastqfile)
        # # samtools; the dash at the end means STDIN
        # samtools_sort_cmd = ('samtools', 'sort', '-@', '12', '-o', fileprefix + '_unfilteredtranscriptome.bam', '-')
        # samtools_index_cmd = ('samtools', 'index', fileprefix + '_unfilteredtranscriptome.bam')
        # pipettor.run([mm2_cmd, samtools_sort_cmd])
        # pipettor.run([samtools_index_cmd])
        #
        # ##filter transcriptome alignment to chimeric only and remove the rest
        # ##run filtering
        # samfile = pysam.AlignmentFile(fileprefix + '_unfilteredtranscriptome.bam', "rb")
        # withsup = pysam.AlignmentFile(fileprefix + '_transcriptomechimeric.bam', "wb", template=samfile)
        # for read in samfile.fetch():
        #     if read.is_mapped and not read.is_secondary:
        #         if read.has_tag('SA'):
        #             withsup.write(read)
        # samfile.close()
        # withsup.close()
        # pysam.index(fileprefix + '_transcriptomechimeric.bam')
        #
        # pipettor.run([('rm', fileprefix + '_unfilteredtranscriptome.bam', fileprefix + '_unfilteredtranscriptome.bam.bai')])

        # transcriptomechimbam = fileprefix + '_transcriptomechimeric.bam'
        # genomechimbam = '../' + fileprefix + '.hg38aligned_chimeric.bam'
        # print(os.path.exists(genomechimbam))
        # fusioncmd = 'python3 /private/groups/brookslab/cafelton/git-flair/flair/src/flair/flair_detectfusions_copy.py -g ' \
        #             + reference_genome + ' -f ' + reference_annot + ' -r ' + '../' + fastqfile + ' -b ' + genomechimbam \
        #             + ' -o ' + fileprefix + ' --annotated_fa ../gencode.v38.flair.fa -t 12 -s 2'
        # subprocess.call(fusioncmd, shell=True)
        # os.chdir('..')

        # break
    # break
