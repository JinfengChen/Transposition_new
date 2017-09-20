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
python RunRelocaTEi_bam_Pop.py --input /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam/

Run RelocaTEi on population level data.
--input: directory of bam files, GN160.bam
--genome: reference genome sequence, need to index using bwa first
--repeat: fasta repeat sequence, do not have too similar element like mPing and ping, which are >99% identity
--existingRun: result direrctory that ran on the same dataset. we can link subset.bam to new project, save time if we run mPing first and run other TE after.
--project: output directory

    '''
    print message

def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 80 --lines %s --interval 120 --resource nodes=1:ppn=1,walltime=100:00:00,mem=12G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

#readgroup:RIL1_0_CGTACG_FC153L5	platform:Illumina	map:../input/bam/GN1.bam	readlen:100.18	lib:RIL1_0_CGTACG_FC153L5	num:10001	lower:76.38	upper:303.56	mean:141.39	std:25.74	SWnormality:-71.20	exe:samtools view
def getsize(infile):
    size = [] 
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                size = [unit[8], unit[9], unit[6], unit[7]]
                for i in range(len(size)):
                    size[i] = re.sub(r'.*:', '', size[i]) 
    return size



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-g', '--genome')
    parser.add_argument('-r', '--repeat')
    parser.add_argument('-e', '--existingRun')
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.genome:
        args.genome = '/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU_r7.fa'
  
    if not args.repeat:
        args.repeat = '/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa'
    
    if not args.project:
        args.project = 'Rice_pop'
    
    #-t ../input/mping_UNK.fa -g /rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa -d ../input/FC52_7 -e HEG4 -o mPing_HEG4_UNK -r 1 -p 1 -a 1   
    RelocaTE = 'python /rhome/cjinfeng/BigData/00.RD/RelocaTE2_mPing/scripts/relocaTE.py'
    Reference   = os.path.abspath(args.genome)
    Repeat      = os.path.abspath(args.repeat)
    existingTE  = '%s.mPing.RepeatMasker.out' %(Reference)
    if not os.path.isfile(existingTE):
        print 'Will give many false positive on existing TEs if you do not give RepeatMasker annotation: Reference_name.RepeatMasker.out'
        
    project = '%s_RelocaTEi' %(os.path.abspath(args.project))
    if not os.path.exists(project):
        os.mkdir(project)
    ofile = open ('%s.run.sh' %(project), 'w')
    #print project
    bams  = glob.glob('%s/*.bam' %(args.input))
    count = 0
    #print bams
    for bam in sorted(bams):
        count += 1
        outdir = '%s/%s_RelocaTEi' %(project, os.path.split(os.path.splitext(bam)[0])[1])
        subdata= '%s/RelocaTEi_%s/repeat/fastq' %(os.path.abspath(args.existingRun), os.path.split(os.path.splitext(bam)[0])[1]) if args.existingRun else ''
        print outdir
        if not os.path.exists('%s/%s' %(project, outdir)):
            prefix   = os.path.split(os.path.splitext(bam)[0])[1]
            cfg      = '%s/%s.config' %(args.input, prefix) 
            size     = getsize(cfg)
            relocaTE = '%s --te_fasta %s --genome_fasta %s --bam %s --outdir %s --reference_ins %s --size %s' %(RelocaTE, Repeat, Reference, bam, outdir, existingTE, size[0])
            shell    = 'time (bash %s/run_these_jobs.sh > %s/run.log 2>> %s/run.log) 2>> %s/run.log' %(outdir, outdir, outdir, outdir)
            print >> ofile, relocaTE
            print >> ofile, shell
            #os.system(relocaTE)
            #start_time = time.time()
            #os.system(shell)
            #end_time   = time.time()
            #ofile = open('./%s/run.log' %(outdir), 'a')
            #print >> ofile, 'Run Time: %s seconds' %(end_time - start_time)
            #ofile.close()
            #os.system(mv)
            #print relocaTE
            #print shell
    ofile.close()
    runjob('%s.run.sh' %(project), 2)
 
if __name__ == '__main__':
    main()

