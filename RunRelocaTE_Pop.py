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
python RunRelocaTE_Pop.py --input /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_fastq

Run RelocaTEi on population level data.
--input: directory of RIL subdir: RILs_ALL_fastq, which have RILs_ALL_fastq/RIL10/RIL10_1.fq
--genome: reference genome sequence, need to index using bwa first
--repeat: fasta repeat sequence, do not have too similar element like mPing and ping, which are >99% identity
--project: output directory

    '''
    print message

def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource walltime=100:00:00,mem=5G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-g', '--genome')
    parser.add_argument('-r', '--repeat')
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.genome:
        args.genome = '/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa'
  
    if not args.repeat:
        args.repeat = '/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/TE/mping.fa'
    
    if not args.project:
        args.project = 'Rice_pop'
    
    #-t ../input/mping_UNK.fa -g /rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa -d ../input/FC52_7 -e HEG4 -o mPing_HEG4_UNK -r 1 -p 1 -a 1   
    RelocaTE    = 'perl /rhome/cjinfeng/software/tools/RelocaTE_1.0.3/RelocaTE/scripts/relocaTE.pl'
    RelocaTE_bin= '/rhome/cjinfeng/software/tools/RelocaTE_1.0.3_i/RelocaTE/scripts'
    samtools    = '/usr/local/bin/samtools'
    bamdir      = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam'
    Reference   = os.path.abspath(args.genome)
    Repeat      = os.path.abspath(args.repeat)
    existingTE  = '/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Nipponbare.mPing.txt'
    if not os.path.isfile(existingTE):
        print 'Will give many false positive on existing TEs if you do not give RepeatMasker annotation: Reference_name.RepeatMasker.out'
        
    project = os.path.abspath(args.project)
    if not os.path.exists(project):
        os.mkdir(project)
    ofile = open ('%s.RelocaTE.run.sh' %(project), 'w')
    #print project
    read_dirs = glob.glob('%s/RIL*' %(args.input))
    count = 0
    for read_dir in sorted(read_dirs):
        count += 1
        ril    = os.path.split(read_dir)[1]
        bam    = '%s/%s.bam' %(bamdir, ril) 
        outdir = os.path.abspath('RelocaTE_%s' %(ril))
        print outdir
        if not os.path.exists('./%s/%s' %(project, outdir)):
            relocaTE = '%s -1 _1 -2 _2 -t %s -g %s -d %s -e %s -o %s -r %s -p 1 -a 0' %(RelocaTE, Repeat, Reference, read_dir, ril, outdir, existingTE)
            shell    = 'time (bash %s/run_these_jobs.sh > %s/run.log 2>> %s/run.log) 2>> %s/run.log' %(outdir, outdir, outdir, outdir)
            nonref   = 'grep "Non-ref" %s/mping/results/%s.mping.all_inserts.gff > %s/mping/results/%s.mping.all_inserts.Non_Ref.gff' %(outdir, ril, outdir, ril)
            shared   = 'grep "Reference-only" %s/mping/results/%s.mping.all_inserts.gff > %s/mping/results/%s.mping.all_inserts.Ref_Only.gff' %(outdir, ril, outdir, ril)
            refonly  = 'grep "Shared" %s/mping/results/%s.mping.all_inserts.gff > %s/mping/results/%s.mping.all_inserts.Shared.gff' %(outdir, ril, outdir, ril)
            charac   = 'perl %s/characterizer.pl -s %s/mping/results/%s.mping.all_nonref.txt -b %s -g %s --samtools %s' %(RelocaTE_bin, outdir, ril, bam, Reference, samtools)
            mv       = 'mv %s %s' %(outdir, project)
            print >> ofile, relocaTE
            print >> ofile, shell
            print >> ofile, nonref
            print >> ofile, shared
            print >> ofile, refonly
            print >> ofile, charac
            print >> ofile, mv
    ofile.close()
    runjob('%s.RelocaTE.run.sh' %(project), 7)
 
if __name__ == '__main__':
    main()

