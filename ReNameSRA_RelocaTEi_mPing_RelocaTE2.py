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
    #cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl -q js --maxjob 5 --lines %s --interval 120 --resource nodes=1:ppn=8,walltime=100:00:00,mem=20G --convert no %s' %(lines, script)
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 5 --lines %s --interval 120 --task 8 --mem 20G --time 100:00:00 --convert no %s' %(lines, script)
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
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data


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
        args.output = '%s_RelocaTEi' %(os.path.abspath(args.input))

    if not args.genome:
        args.genome = '/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU_r7.fa'
  
    if not args.repeat:
        args.repeat = '/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa'
        #args.repeat = '/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/Rice.TE.short.unique.fa'

    #-t ../input/mping_UNK.fa -g /rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa -d ../input/FC52_7 -e HEG4 -o mPing_HEG4_UNK -r 1 -p 1 -a 1   
    #RelocaTE = 'python /rhome/cjinfeng/software/tools/RelocaTE_1.0.3_i/RelocaTE/scripts/relocaTE.py'
    RelocaTE = 'python /rhome/cjinfeng/BigData/00.RD/RelocaTE2/scripts/relocaTE2.py'
    #RelocaTE = 'python /rhome/cjinfeng/BigData/00.RD/RelocaTE2_mPing/scripts/relocaTE.py'
    Reference= os.path.abspath(args.genome)
    Repeat   = os.path.abspath(args.repeat)
    project = os.path.split(args.output)[1]
    cpu = 8
    if not os.path.exists(project):
        os.mkdir(project)
    print project
    read_dirs = glob.glob('%s/RIL*' %(os.path.abspath(args.input)))
    ofile = open('%s.run.sh' %(args.output), 'w')
    for read_dir in sorted(read_dirs):
        outdir = '%s/%s_RelocaTEi' %(os.path.abspath(args.output), os.path.split(read_dir)[1])
        existingTE  = '%s.mPing.RepeatMasker.out' %(Reference)
        # relocate will not run if there is result exists
        if not os.path.exists(outdir):
        #if 1:
            #bug in fq and fq.gz format, RIL232_3_ACTTGA_FC381L1_p1.te_repeat.read_repeat_name.chr.txt
            relocaTE = '%s --mate_1_id _p1 --mate_2_id _p2 --split --te_fasta %s --genome_fasta %s --fq_dir %s --outdir %s --reference_ins %s --size 500 --step 1234567 --mismatch 2 --run --cpu %s --aligner blat --verbose 3' %(RelocaTE, Repeat, Reference, read_dir, outdir, existingTE, cpu)
            #relocaTE = '%s --mate_1_id _p1 --mate_2_id _p2 --te_fasta %s --genome_fasta %s --fq_dir %s --outdir %s --reference_ins %s' %(RelocaTE, Repeat, Reference, read_dir, outdir, existingTE)
            #shell    = 'bash %s/run_these_jobs.sh > %s/run.log 2> %s/run.log2' %(outdir, outdir, outdir)
            #os.system(relocaTE)
            print >> ofile, relocaTE
            #print >> ofile, shell
    ofile.close()
    runjob('%s.run.sh' %(args.output), 1)
 
if __name__ == '__main__':
    main()

