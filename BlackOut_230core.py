#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser, createdir

def usage():
    test="name"
    message='''
python BlackOut_230Core.py --input Bam.Core.blacklist

Generate files from 275 to 230 core by removing problem RILs.
RIL275_RelocaTEi.CombinedGFF.ALL.gff
RIL275_RelocaTEi.CombinedGFF.characterized.gff
RIL275_RelocaTEi.summary.table
RIL275_RelocaTEi.summary_clean.table
RIL275_RelocaTE.sofia.ping_code.table
RIL275_RelocaTE.CombinedGFF.Ref_only.gff
RIL275_RelocaTE.CombinedGFF.Shared.gff

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#RIL104
def read_blacklist(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[0]] = 1
    return data


#blackout rils in gff and write a new file
#Chr1    RIL231_0        transposable_element_attribute  4228091 4228092
def gff_blackout(infile, blacklist):
    new_file = re.sub(r'275', r'230', infile)
    print 'new table: %s' %(new_file)
    ofile = open(new_file, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                if not line.startswith(r'#'):
                    unit = re.split(r'\t',line)
                    ril  = re.split(r'_', unit[1])[0]
                    if not blacklist.has_key(ril):
                        print >> ofile, line
                else:
                    print >> ofile, line
    ofile.close()

#blackout rils in table and write a new file
#RIL1    141.39  25.74   6.35    6.19    0.98    336     228     108     226     209     17      0
def table_blackout(infile, blacklist):
    new_file = re.sub(r'275', r'230', infile)
    print 'new table: %s' %(new_file)
    ofile = open(new_file, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                if line.startswith(r'RIL'):
                    unit = re.split(r'\t',line)
                    if not blacklist.has_key(unit[0]):
                        print >> ofile, line
                else:
                    print >> ofile, line
    ofile.close()

#blackout rils in table and write a new file
#Pings   Ping_Code       RIL
#NA      NA      RIL158
#NA      NA      RIL242
#0       NA      RIL39
def table_blackout_ping(infile, blacklist):
    new_file = re.sub(r'275', r'230', infile)
    print 'new table: %s' %(new_file)
    ofile = open(new_file, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                if not line.startswith(r'Pings'):
                    unit = re.split(r'\t',line)
                    if not blacklist.has_key(unit[2]):
                        print >> ofile, line
                else:
                    print >> ofile, line
    ofile.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    #RIL275_RelocaTEi.CombinedGFF.ALL.gff
    #RIL275_RelocaTEi.CombinedGFF.characterized.gff
    #RIL275_RelocaTEi.summary.table
    #RIL275_RelocaTEi.summary_clean.table
    #RIL275_RelocaTE.sofia.ping_code.table
    #RIL275_RelocaTE.CombinedGFF.Ref_only.gff
    #RIL275_RelocaTE.CombinedGFF.Shared.gff
    blacklist = read_blacklist(args.input)
    gff_blackout('RIL275_RelocaTEi.CombinedGFF.characterized.gff', blacklist)
    gff_blackout('RIL275_RelocaTEi.CombinedGFF.ALL.gff', blacklist)
    gff_blackout('RIL275_RelocaTE.CombinedGFF.Ref_only.gff', blacklist)
    gff_blackout('RIL275_RelocaTE.CombinedGFF.Shared.gff', blacklist)
    
    table_blackout('RIL275_RelocaTEi.summary.table', blacklist)
    table_blackout('RIL275_RelocaTEi.summary_clean.table', blacklist)
    table_blackout_ping('RIL275_RelocaTE.sofia.ping_code.table', blacklist)    

if __name__ == '__main__':
    main()

