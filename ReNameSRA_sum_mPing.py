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
import fnmatch

def usage():
    test="name"
    message='''
python ReNameSRA_sum_Ping.py --input Other_fastq_RelocaTEi_Ping

Summary RelocatTE call in current direcory using simulation in input directory
--call: RelocaTE, TEMP or other
--check: check unfinished job without results directory or non_ref.gff
--input: dir of simulation, where we can find insertion simulated gff file "MSU7.Chr4.mPing.rep1.gff"

    '''
    print message

def readtable(infile):
    data = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                data.append(line)
    return data

def write_file(fh, lines):
    for line in lines:
        print >> fh, line


def write_gff_file(outfile, infiles):
    ofile = open(outfile, 'w')
    for infile in infiles:
        if os.path.exists(infile):
            lines = readtable(infile)
            write_file(ofile, lines)
    ofile.close()

#Right_junction_reads:8;Left_junction_reads
#
def parse_gff(infile, left_cut=0, right_cut=0):
    num = 0
    if not os.path.isfile(infile) or not os.path.getsize(infile) or int(os.path.getsize(infile)) == 0:
        return 0 
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'):
                #print line
                #num += 1
                unit  = re.split(r'\t',line)
                start = int(unit[3]) 
                end   = int(unit[4])
                chro  = unit[0]
                strand= unit[6]
                temp  = defaultdict(str)
                attrs = re.split(r';', unit[8])
                for attr in attrs:
                    #print attr
                    if not attr == '':
                        #print 'yes'
                        idx, value = re.split(r'\=|\:', attr)
                        temp[idx] = value
                #print left_cut, right_cut, temp['Left_junction_reads'], temp['Right_junction_reads']
                if int(temp['Right_junction_reads']) >= right_cut and int(temp['Left_junction_reads']) >= left_cut:
                    #print 'both end'
                    num += 1
                else:
                    #print 'one end'
                    continue
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
        args.list = '/bigdata/wesslerlab/shared/Rice/Rice_population_sequence/Rice_3000/GigaScience/rice_line_3000.download.list'

    #inf  = parse_inf(args.list, args.call)
    inf = defaultdict(lambda : str())
    data = defaultdict(lambda : defaultdict(lambda : list)) 
    project = os.path.split(args.input)[1]
    sum_file = '%s.summary' %(args.input)
    ofile = open(sum_file, 'w')
    print >> ofile, 'Accession\tPing\tPing_NonRef\tPing_Ref\tName\tOrigin\tClass'
    for call in sorted(os.listdir(args.input)):
        if not fnmatch.fnmatch(call, '*_RelocaTEi'):
            continue
        ping_dir = re.sub(r'_mPing$', r'_Ping', args.input)
        pong_dir = re.sub(r'_mPing$', r'_Pong', args.input)
        dirname  = os.path.abspath('%s/%s' %(args.input, call))
        ping_dir = os.path.abspath('%s/%s' %(ping_dir, call))
        pong_dir = os.path.abspath('%s/%s' %(pong_dir, call))
        #print 'mPing: %s' %(dirname)
        #print 'Ping: %s' %(ping_dir)
        #print 'Pong: %s' %(pong_dir)
        if args.call == 'RelocaTE':
            pass
        elif args.call == 'RelocaTEi':
            mping_nonref_raw_gff = '%s/repeat/results/ALL.all_nonref_insert.gff' %(dirname)
            mping_nonref_ft_gff  = '%s/repeat/results/ALL.all_nonref_insert.no_ping_pong.gff' %(dirname)
            mping_ref_raw_gff    = '%s/repeat/results/ALL.all_ref_insert.gff' %(dirname)
            mping_ref_ft_gff     = '%s/repeat/results/ALL.all_ref_insert.no_ping_pong.gff' %(dirname)
            ping_nonref_gff      = '%s/repeat/results/ALL.all_nonref_insert.Ping.gff' %(ping_dir)
            ping_ref_gff         = '%s/repeat/results/ALL.all_ref_insert.Ping.gff' %(ping_dir)
            pong_nonref_gff      = '%s/repeat/results/ALL.all_nonref_insert.Pong.gff' %(pong_dir)
            pong_ref_gff         = '%s/repeat/results/ALL.all_ref_insert.Pong.gff' %(pong_dir)
            nonref_ft_gff        = '%s/repeat/results/ALL.all_nonref_insert.ping_pong.gff' %(dirname)
            ref_ft_gff           = '%s/repeat/results/ALL.all_ref_insert.ping_pong.gff' %(dirname)
            call_inf    = inf[call] if inf.has_key(call) else 'NA\tNA\tNA'
            if os.path.exists('%s/repeat/results' %(dirname)):
                #print 'have results' 
                write_gff_file(nonref_ft_gff, [ping_nonref_gff, pong_nonref_gff])
                write_gff_file(ref_ft_gff, [ping_ref_gff, pong_ref_gff])
                cmd1 = '/opt/bedtools/2.17.0-25-g7b42b3b/bin/bedtools window -w 100 -a %s -b %s -v > %s' %(mping_nonref_raw_gff, nonref_ft_gff, mping_nonref_ft_gff)
                cmd2 = '/opt/bedtools/2.17.0-25-g7b42b3b/bin/bedtools window -w 100 -a %s -b %s -v > %s' %(mping_ref_raw_gff, ref_ft_gff, mping_ref_ft_gff)
                #print cmd1
                #print cmd2
                #nonref
                if parse_gff(nonref_ft_gff) > 0: 
                    os.system(cmd1)
                else:
                    os.system('cp %s %s' %(mping_nonref_raw_gff, mping_nonref_ft_gff))
                #ref
                if parse_gff(ref_ft_gff) > 0:
                    os.system(cmd2)
                else:
                    os.system('cp %s %s' %(mping_ref_raw_gff, mping_ref_ft_gff))
                call_nonref_mping  = parse_gff(mping_nonref_ft_gff, 1, 1)
                call_ref_mping     = parse_gff(mping_ref_ft_gff, 1, 1)
                call_mping         = call_nonref_mping + call_ref_mping
                print >> ofile, '%s\t%s\t%s\t%s\t%s' %(call, call_mping, call_nonref_mping, call_ref_mping, call_inf)
            else:
                print 'no results'
                print >> ofile, '%s\tNA\tNA\tNA\t%s' %(call, call_inf)
            #cmd_ping    = 'python PickPing.py --input %s' %(dirname)
            #ping_nonref_gff  = '%s/repeat/results/ALL.all_nonref_insert.Ping.gff' %(dirname)
            #ping_ref_gff     = '%s/repeat/results/ALL.all_ref_insert.Ping.gff' %(dirname)
            #os.system(cmd_ping) 
            #non_ref_gff = '%s/repeat/results/ALL.all_nonref_insert.gff' %(dirname)
            #ref_gff     = '%s/repeat/results/ALL.all_ref_insert.gff' %(dirname)
            #call_nonref_ping  = parse_gff(ping_nonref_gff)
            #call_ref_ping     = parse_gff(ping_ref_gff)
            #call_ping  = call_nonref_ping + call_ref_ping
            #call_inf    = inf[call] if inf.has_key(call) else 'NA\tNA\tNA' 
            #if int(args.check) == 1 and (not os.path.exists('%s/repeat/results' %(dirname)) or not os.path.getsize(non_ref_gff) > 0):
            #    print >> ofile, '%s\t%s\t%s\t%s\t%s\tNo output' %(call, call_ping, call_ref, call_nonref, call_inf)
            #print >> ofile, '%s\t%s\t%s\t%s\t%s' %(call, call_ping, call_nonref_ping, call_ref_ping, call_inf)
        elif args.call == 'TEMP':
            pass
    ofile.close()  

if __name__ == '__main__':
    main()

