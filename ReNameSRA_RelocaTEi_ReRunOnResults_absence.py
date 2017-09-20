#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
import time
from Bio import SeqIO

def usage():
    test="name"
    message='''
python ReNameSRA_RelocaTEi.py --input Japonica_fastq

Run RelocaTEi for rice strain in Japonica_fastq

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 80 --lines %s --interval 120 --resource nodes=1:ppn=1,walltime=100:00:00,mem=10G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r' ',line)
                data[unit[0]] = line
                #print unit[0], line
    return data

def rerun_shell(outfile, script, topdir):
    script_cmd = readtable(script)
    script_list= re.split(r' ', script_cmd['python'])
    script_list[1] = '/rhome/cjinfeng/BigData/00.RD/RelocaTE2_mPing/scripts/relocaTE_absenceFinder.py'
    ofile = open(outfile, 'w')
    for i in range(1,13):
        script_list[3] = 'Chr%s' %(i)
        cmd = ' '.join(script_list)
        cmd = re.sub(r'.RepeatMasker.out', r'.mPing.RepeatMasker.out', cmd)
        cmd = re.sub(r'/shared/wesslerlab/', r'/bigdata/wesslerlab/shared/', cmd)
        print >> ofile, cmd
    print >> ofile, 'cat %s/repeat/results/*.all_ref_insert.txt > %s/repeat/results/ALL.all_ref_insert.txt' %(topdir, topdir)
    print >> ofile, 'cat %s/repeat/results/*.all_ref_insert.gff > %s/repeat/results/ALL.all_ref_insert.gff' %(topdir, topdir)
    ofile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-g', '--genome')
    parser.add_argument('-r', '--repeat')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = '%s' %(os.path.abspath(args.input))

    if not args.genome:
        args.genome = '/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU_r7.fa'
  
    if not args.repeat:
        args.repeat = '/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/ping.fa'
        #args.repeat = '/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/Rice.TE.short.unique.fa'

    #-t ../input/mping_UNK.fa -g /rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa -d ../input/FC52_7 -e HEG4 -o mPing_HEG4_UNK -r 1 -p 1 -a 1   
    #RelocaTE = 'python /rhome/cjinfeng/software/tools/RelocaTE_1.0.3_i/RelocaTE/scripts/relocaTE.py'
    #RelocaTE = 'python /rhome/cjinfeng/BigData/00.RD/RelocaTE2/scripts/relocaTE.py'
    RelocaTE = 'python /rhome/cjinfeng/BigData/00.RD/RelocaTE2_mPing/scripts/relocaTE.py'
    Reference= os.path.abspath(args.genome)
    Repeat   = os.path.abspath(args.repeat)
    project = os.path.split(args.output)[1]
    cpu = 16
    if not os.path.exists(project):
        os.mkdir(project)
    print project
    read_dirs = glob.glob('%s/ERS*' %(os.path.abspath(args.input)))
    ofile = open('%s.run.sh' %(args.output), 'w')
    for read_dir in sorted(read_dirs):
        outdir = '%s/%s' %(os.path.abspath(args.output), os.path.split(read_dir)[1])
        existingTE  = '%s.mPing.RepeatMasker.out' %(Reference)
        # relocate will not run if there is result exists
        #if not os.path.exists(outdir):
        if 1:
            #relocaTE = '%s --te_fasta %s --genome_fasta %s --fq_dir %s --outdir %s --reference_ins %s' %(RelocaTE, Repeat, Reference, read_dir, outdir, existingTE)
            os.system('cp /rhome/cjinfeng/Rice/Rice_population_sequence/Rice_3000/CAAS/existingTE.bed %s/repeat/' %(outdir))
            rerun_shell('%s/run_these_jobs_rerun.sh' %(outdir), '%s/shellscripts/step_6/0.repeat.absence.sh' %(outdir), outdir)
            shell    = 'bash %s/run_these_jobs_rerun.sh > %s/run.log 2> %s/run.log2' %(outdir, outdir, outdir)
            #os.system(relocaTE)
            #print >> ofile, relocaTE
            print >> ofile, shell
    ofile.close()
    runjob('%s.run.sh' %(args.output), 1)
 
if __name__ == '__main__':
    main()

