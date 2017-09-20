#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
import glob
import time
from Bio import SeqIO

def usage():
    test="name"
    message='''
python ReNameSRA_failed_call.py --input RILs_ALL_fastq_correct_merged_RelocaTEi_mPing
check if some chr is failed for finished job
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def parse_gff(infile):
    num = 0
    if not os.path.isfile(infile) or not os.path.getsize(infile):
        return 0 
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'):
                num += 1 
    return num

#1       B001    China   Temperate japonica      ERS470219       anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/ERR/ERR622/ERR622583/ERR622583.sra
def parse_inf(infile, sufix):
    data = defaultdict(lambda : str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'):
                unit = re.split(r'\t', line)
                index = '%s_%s' %(unit[4], sufix)
                if not data.has_key(index):
                    data[index]= '%s\t%s\t%s' %(unit[1], unit[2], unit[3])
    return data



def parse_overlap_reloacte(infile, ref_te, call_te):
    #Total number of call, true call, call with breakpoint near TSD, No call, False call
    data = defaultdict(lambda : int)
    dupli= defaultdict(lambda : int())
    true = 0
    tsd  = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                #print line
                if dupli.has_key('%s_%s' %(unit[12], unit[13])):
                    continue
                else:
                    dupli['%s_%s' %(unit[12], unit[13])] == 1
                #print 'pass'
                pos  = map(int, [unit[3], unit[4], unit[12], unit[13]])
                dist_min = min([abs(pos[2]-pos[0]), abs(pos[2]-pos[1]), abs(pos[3]-pos[0]), abs(pos[3]-pos[0])])
                dist_max = max([abs(pos[2]-pos[0]), abs(pos[2]-pos[1]), abs(pos[3]-pos[0]), abs(pos[3]-pos[0])])
                true += 1
                if dist_min <= 5 and dist_max <= 5:
                    tsd += 1
    call = len(call_te.keys())
    ref  = len(ref_te.keys())
    data = [call, true, tsd, ref-true, call-true]
    return data




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input') 
    parser.add_argument('-l', '--list')
    parser.add_argument('-c', '--call')
    parser.add_argument('-ck', '--check')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.check:
        args.check = 0

    if not args.call: 
        args.call = 'RelocaTEi'

    if not args.list:
        args.list = 'rice_line_CAAS_534.download.list'

    #ofile = open(sum_file, 'w')
    #print >> ofile, 'Accession\tPing\tmPing_Ref\tmPing_Non_Ref\tName\tOrigin\tClass'
    for call in sorted(os.listdir(args.input)):
        dirname = os.path.abspath('%s/%s' %(args.input, call))
        #print dirname
        if args.call == 'RelocaTE':
            pass
        elif args.call == 'RelocaTEi':
            #cmd_ping    = 'python PickPing.py --input %s' %(dirname)
            #ping_gff    = '%s/repeat/results/ALL.all_nonref_insert.Ping.gff' %(dirname)
            #if not os.path.exists(ping_gff):
            #    os.system(cmd_ping) 
            #non_ref_gff = '%s/repeat/results/ALL.all_nonref_insert.characTErized.gff' %(dirname)
            #ref_gff     = '%s/repeat/results/ALL.all_ref_insert.gff' %(dirname)
            #call_ping   = parse_gff(ping_gff)
            #call_nonref = parse_gff(non_ref_gff)
            #call_ref    = parse_gff(ref_gff)
            #call_inf    = inf[call] if inf.has_key(call) else 'NA\tNA\tNA'
            if os.path.exists('%s/repeat/results' %(dirname)):
                #print >> ofile, '%s\t%s\t%s\t%s\t%s\tNo output' %(call, call_ping, call_ref, call_nonref, call_inf)
                #os.system('')
                count = 0
                for i in range(1,13):
                    nonref = '%s/repeat/results/Chr%s.repeat.all_nonref_insert.gff' %(dirname, i)
                    ref    = '%s/repeat/results/Chr%s.repeat.all_ref_insert.gff' %(dirname, i)
                    if os.path.exists(nonref) and os.path.exists(ref):
                        count += 1
                if count < 12:
                    print '%s: %s'  %(dirname, count)
                    #tail -n 26 run_these_jobs.sh > run_these_jobs_rerun.sh
                    #bash run_these_jobs_rerun.sh > log 2>&1
                    #os.system('rm -R %s/repeat/results' %(dirname))
                    #os.system('tail -n 26 %s/run_these_jobs.sh > %s/run_these_jobs_rerun.sh' %(dirname, dirname))
                    #os.system('bash %s/run_these_jobs_rerun.sh > %s/run_these_jobs_rerun.log 2>&1' %(dirname, dirname))
            else:
                print 'results not exists: %s' %(dirname)
                #remove = 'rm %s/repeat/fastq/GN*' %(dirname)
                #remove = 'bash %s/clean_intermediate_files.sh' %(dirname)
                #os.system(remove)
            #print >> ofile, '%s\t%s\t%s\t%s\t%s' %(call, call_ping, call_ref, call_nonref, call_inf)
        elif args.call == 'TEMP':
            pass
    #ofile.close()  

if __name__ == '__main__':
    main()

