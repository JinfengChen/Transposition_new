#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import glob

def usage():
    test="name"
    message='''
python PickPing.py --input Japonica_fastq_RelocaTEi_Pong/ERS470370_RelocaTEi


    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#repeat  Chr1:480846..480848     Junction_reads  ERR624499.2788059:start:5,ERR624507.3823168:start:5,ERR624505.1031664:end:5
#repeat  Chr1:480846..480848     Left_supporting_reads   
#repeat  Chr1:480846..480848     Right_supporting_reads  
def parse_support_reads(infile, support_inf):
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if len(unit) == 4:
                    reads = re.split(r',', unit[3])
                    support_inf[unit[1]][unit[2]] = reads
                else:
                    reads = []
                    support_inf[unit[1]][unit[2]] = reads

#Japonica_fastq_RelocaTEi_Pong/ERS470370_RelocaTEi/repeat/results/ALL.all_nonref_insert.gff
#Chr1    not.give        RelocaTE_i      480846  480848  .       -       .       ID=repeat_Chr1_480846_480848;TSD=TTA;Note=Non-reference, not found in reference;
#Right_junction_reads:2;Left_junction_reads:1;Right_support_reads:0;Left_support_reads:0;
def gff_parse(infile, flank_inf, support_inf, pong_NB, ref_pong_gff):
    pong_evidence = defaultdict(lambda : defaultdict(lambda : list()))
    data = defaultdict(lambda : list())
    r = re.compile(r'\=')
    ping_gff = re.sub(r'.gff', r'.Pong.gff', infile)
    ping_log = re.sub(r'.gff', r'.Pong.inf', infile)
    #print infile
    #print ping_gff
    ofile = open(ping_gff, 'w')
    ofile_log = open(ping_log, 'w')
    ping_special = 0
    with open (infile, 'r') as filehd:
        for line in filehd: 
            line = line.rstrip()
            if len(line) > 2 and not line.startswith('#'): 
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
                        if r.search(attr):
                            idx, value = re.split(r'\=', attr)
                            temp[idx] = value
                        else:
                            idx, value = re.split(r'\:', attr)
                            temp[idx] = value
                Rjun    = int(temp['Right_junction_reads'])
                Rsup    = int(temp['Right_support_reads'])
                Ljun    = int(temp['Left_junction_reads'])
                Lsup    = int(temp['Left_support_reads'])
#####################ref pong
                if pong_NB.has_key('%s:%s' %(chro, start)):
                    #print 'left: %s' %(line)
                    pong_evidence[pong_NB['%s:%s' %(chro, start)]['loci']][pong_NB['%s:%s' %(chro, start)]['direction']] = [Rjun, Ljun, Rsup, Lsup]
                    pong_evidence[pong_NB['%s:%s' %(chro, start)]['loci']]['strand'] = [strand]
                    pong_evidence[pong_NB['%s:%s' %(chro, start)]['loci']]['inf']    = pong_NB['%s:%s' %(chro, start)]['loci_i']
                    #print 'left', pong_evidence[pong_NB['%s:%s' %(chro, start)]['loci']]['left']
                    continue
                elif pong_NB.has_key('%s:%s' %(chro, start-1)):
                    #print 'right: %s' %(line)
                    pong_evidence[pong_NB['%s:%s' %(chro, start-1)]['loci']][pong_NB['%s:%s' %(chro, start-1)]['direction']] = [Rjun, Ljun, Rsup, Lsup]
                    pong_evidence[pong_NB['%s:%s' %(chro, start-1)]['loci']]['strand'] = strand
                    pong_evidence[pong_NB['%s:%s' %(chro, start-1)]['loci']]['inf']    = pong_NB['%s:%s' %(chro, start-1)]['loci_i']
                    #print 'right', pong_evidence[pong_NB['%s:%s' %(chro, start-1)]['loci']]['right']
                    continue
#####################ref pong
                #print '%s\t%s\t%s\t%s\t%s\t%s' %(repid, chro, start, end, repname, repfam)
                #print '%s\t%s\t%s' %(chro, start, end)
                repid   = '%s:%s..%s' %(chro, start, end)
                print >> ofile_log, 'Pong: %s' %(repid)
                left_supporting_picks  = check_supporting_reads(support_inf[repid]['Left_supporting_reads'], repid, flank_inf, ofile_log)
                right_supporting_picks = check_supporting_reads(support_inf[repid]['Right_supporting_reads'], repid, flank_inf, ofile_log)
                left_junction_picks  = check_supporting_reads(support_inf[repid]['Left_junction_reads'], repid, flank_inf, ofile_log)
                right_junction_picks = check_supporting_reads(support_inf[repid]['Right_junction_reads'], repid, flank_inf, ofile_log)
                #print left_supporting_picks[repid], right_supporting_picks[repid]
                #print left_junction_picks[repid], right_junction_picks[repid]
                #print line
                if len(left_supporting_picks[repid]) + len(left_junction_picks[repid]) > 0 and len(right_supporting_picks[repid]) + len(right_junction_picks[repid]) > 0:
                    #print 'Find Pong: %s' %(line)
                    if (chro == 'Chr1' and  start == 38092038) or (chro == 'Chr12' and start == 24587180):
                        unit[0] = 'Chr1'
                        unit[3] = '2640500'
                        unit[4] = '2640502'
                        line    = '\t'.join(unit)
                        if ping_special == 0:
                            print >> ofile, line
                            ping_special = 1
                    else:
                        print >> ofile, line 
                elif len(left_supporting_picks[repid]) + len(left_junction_picks[repid]) > 0 or len(right_supporting_picks[repid]) + len(right_junction_picks[repid]) > 0:
                    #print 'Candidate Pong: %s' %(line)
                    if (chro == 'Chr1' and  start == 38092038) or (chro == 'Chr12' and start == 24587180):
                        unit[0] = 'Chr1'
                        unit[3] = '2640500'
                        unit[4] = '2640502'
                        line    = '\t'.join(unit)
                        if ping_special == 0:
                            print >> ofile, line
                            ping_special = 1
                    else:
                        print >> ofile, line
                else:
                    pass
    ofile.close()
    ofile_log.close()
###############ref pong
    ofile = open(ref_pong_gff, 'a')
    for pong in sorted(pong_evidence.keys()):
        support = [0,0,0,0]
        #print pong
        if pong_evidence[pong].has_key('left'):
            support[0] += pong_evidence[pong]['left'][0]
            support[1] += pong_evidence[pong]['left'][1]
            support[2] += pong_evidence[pong]['left'][2]
            support[3] += pong_evidence[pong]['left'][3]
            #print 'left', support
        if pong_evidence[pong].has_key('right'):
            support[0] += pong_evidence[pong]['right'][0]
            support[1] += pong_evidence[pong]['right'][1]
            support[2] += pong_evidence[pong]['right'][2]          
            support[3] += pong_evidence[pong]['right'][3]
            #print 'right', support
        chro, start, end = pong_evidence[pong]['inf']
        strand = pong_evidence[pong]['strand']
        pong_inf = '%s\tnot.give\tRelocaTE_i\t%s\t%s\t.\t%s\t.\tID=repeat_%s_%s_%s;TSD=TSD;Note=Shared, in ref and reads;Right_junction_reads:%s;Left_junction_reads:%s;Right_support_reads:%s;Left_support_reads:%s;' %(chro, start, end, strand, chro, start, end, support[0], support[1], support[2], support[3]) 
        print >> ofile, pong_inf
    ofile.close()
###################ref pong
    return data
#Chr6    not.give        RelocaTE_i      12415640        12425974        .       -       .       ID=repeat_Chr6_12415640_12425974;TSD=TSD;Note=Shared, in ref and reads;Right_junction_reads:5;Left_junction_reads:1;Right_support_reads:11;



##start means that the TE was removed from the start of the read
##5 means the trimmed end mapps to the 5prime end of the TE
##3 means the trimmed end mapps to the 3prime end of the TE
def check_supporting_reads(reads_list, repid, flank_inf, ofile_log):
    data = defaultdict(lambda : list())
    #cover_int_ping = 0
    start_int      = 23    #253 for Ping
    end_int        = 5320  #5164 for Ping
    start_snp      = 17    #16 for Ping
    #start_int      = 253
    #end_int        = 5164
    #start_snp      = 16
    if len(reads_list) > 0:
        for read in reads_list:
            for i in [1, 2]:
                if not flank_inf[read][str(i)]['pos'] == '':
                    if flank_inf[read][str(i)]['pos'] == 'middle':
                        #if repid == 'Chr1:480846..480848':
                        #    print '%s\t%s\t%s\t%s' %(repid, read, flank_inf[read][str(i)]['start'], flank_inf[read][str(i)]['end'])
                        #if int(flank_inf[read][str(i)]['start']) >= start_int+10 and int(flank_inf[read][str(i)]['start']) <= end_int-10:
                        #    data[repid].append(read)
                        #elif int(flank_inf[read][str(i)]['end']) >= start_int+10 and int(flank_inf[read][str(i)]['end']) <= end_int-10:
                        #    data[repid].append(read)
                        #if int(flank_inf[read][str(i)]['start']) >= start_int+10 and int(flank_inf[read][str(i)]['start']) <= end_int-10:
                             #reads start located in ping specific region
                        #     data[repid].append(read)
                        #     print >> ofile_log, 'Midlle reads: %s' %(read)
                        #elif int(flank_inf[read][str(i)]['end']) >= start_int+10 and int(flank_inf[read][str(i)]['end']) <= end_int-10:
                             #reads end located in ping specific region
                        #     data[repid].append(read)
                        #     print >> ofile_log, 'Midlle reads: %s' %(read)
                        flank = 5
                        if int(flank_inf[read][str(i)]['start']) >= start_int+flank and int(flank_inf[read][str(i)]['start']) <= end_int-flank and int(flank_inf[read][str(i)]['end']) >= start_int+flank and int(flank_inf[read][str(i)]['end']) <= end_int-flank and int(flank_inf[read][str(i)]['mismatch']) <= 1:
                             data[repid].append(read)
                             print >> ofile_log, 'Midlle reads: %s' %(read)
                    #    elif int(flank_inf[read][str(i)]['start']) < start_snp - 10 and int(flank_inf[read][str(i)]['end']) > start_snp + 10 and int(flank_inf[read][str(i)]['mismatch']) == 0:
                    #         #reads cover snp
                    #         data[repid].append(read)
                    #         print >> ofile_log, 'SNP reads: %s' %(read)
                    #elif flank_inf[read][str(i)]['pos'] == 'start' or flank_inf[read][str(i)]['pos'] == 'end':
                        #junction reads
                    #    if int(flank_inf[read][str(i)]['end']) == 5: #read mapped to 5primer of TE
                    #        if int(flank_inf[read][str(i)]['start']) < start_snp - 10 and int(flank_inf[read][str(i)]['end']) > start_snp + 10 and int(flank_inf[read][str(i)]['mismatch']) == 0:
                    #            #junction reads cover snp 
                    #            data[repid].append(read)
                    #            print >> ofile_log, 'SNP reads: %s' %(read) 
                    #        elif int(flank_inf[read][str(i)]['end']) < start_snp - 10 and int(flank_inf[read][str(i)]['start']) > start_snp + 10 and int(flank_inf[read][str(i)]['mismatch']) == 0:
                    #            data[repid].append(read)
                    #            print >> ofile_log, 'SNP reads: %s' %(read)
                                #junction reads cover snp
        return data
    else:
        return data

#parse blatout file, if the reads were found in flankingreads we store match information also for further check
#Japonica_fastq_RelocaTEi_Pong/ERS470370_RelocaTEi/repeat/blat_output/ERR624495_1.te_repeat.blatout
def parse_blatout(blat_file, flank_inf):
    r        = re.compile(r'.*\/(\w+)\_(\d+)\.te\_repeat\.blatout')
    lib_n, lib_t = ['NA','0']
    if r.search(blat_file):
        lib_n = r.search(blat_file).groups(0)[0]
        lib_t = r.search(blat_file).groups(0)[1]
    #print 'blatout', lib_n, lib_t
    ##align_file
    count = 0
    with open (blat_file, 'r') as filehd:
        #for i in range(5):
        #    next(filehd)
        for line in filehd:
            count += 1
            if count <= 5:
                continue
            line = line.rstrip()
            unit = re.split(r'\t',line)
            match    = int(unit[0])
            mismatch = int(unit[1])
            strand   = unit[8]
            qName    = unit[9]
            qLen     = int(unit[10])
            qStart   = int(unit[11])
            qEnd     = int(unit[12]) - 1
            tName    = unit[13]
            tLen     = int(unit[14])
            tStart   = int(unit[15])
            #get all values into 1st base = 0 postion notation
            tEnd     = int(unit[16]) - 1
            if int(unit[17]) > 1:
                tStart = 0
                tEnd   = 0 
            boundary = 1 if int(qStart) == 0 or int(qEnd) + 1 == int(qLen) else 0
            #if qName == 'ERR624495.1563':
            #    print qName, lib_t, flank_inf[qName][str(lib_t)]['pos']
            if not flank_inf[qName][str(lib_t)]['pos'] == '':
                flank_inf[qName][str(lib_t)]['start'] = str(tStart)
                flank_inf[qName][str(lib_t)]['end']   = str(tEnd)
                flank_inf[qName][str(lib_t)]['mismatch'] = str(mismatch)
                #print qName, lib_t, tStart, tEnd, mismatch

#parse flankingReads.fq file to store these reads that matched to repeat confidently
#Japonica_fastq_RelocaTEi_Pong/ERS470370_RelocaTEi/repeat/flanking_seq/ERR624495_1.te_repeat.flankingReads.fq
#@ERR624495.2984685:middle
#@ERR624495.147060:start:5
def parse_flanking_fastq(fastq_file, flank_inf):
    r_middle = re.compile(r'(.*)\:middle')
    r_end    = re.compile(r'(.*)\:(start|end)\:(\d+)')
    r        = re.compile(r'.*\/(\w+)\_(\d+)\.te\_repeat\.flankingReads\.fq')
    lib_n, lib_t = ['NA','0']
    if r.search(fastq_file):
        lib_n = r.search(fastq_file).groups(0)[0]
        lib_t = r.search(fastq_file).groups(0)[1]
    count    = 3
    #print lib_n, lib_t
    with open (fastq_file, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            count += 1
            if count%4 == 0:
                read_name = str(line[1:])
                #print read_name
                if r_middle.search(read_name):
                    read_id = r_middle.search(read_name).groups(0)[0]
                    pos     = 'middle'
                    flank_inf[read_id][str(lib_t)]['pos'] = pos
                    #print read_id, lib_t, pos
                elif r_end.search(read_name):
                    read_id = r_end.search(read_name).groups(0)[0]
                    pos     = r_end.search(read_name).groups(0)[1]
                    end     = r_end.search(read_name).groups(0)[2]
                    flank_inf[read_id][str(lib_t)]['pos'] = pos
                    flank_inf[read_id][str(lib_t)]['end'] = end
                    #print read_id, lib_t, pos
                else:
                    pass

def setpong():
    pong_dict = defaultdict(lambda : defaultdict(lambda : str()))
    #pong['Chr11']['11436715']['direction'] = 'left'
    #pong['Chr11']['11441880']['direction'] = 'right'
    pongs = ['Chr11:11436715..11441880', 'Chr2:19904309..19909474', 'Chr6:12415640..12425974', 'Chr6:12415820..12420985', 'Chr6:21720706..21725871', 'Chr9:11302206..11307371']
    r = re.compile(r'(\w+):(\d+)\.\.(\d+)')
    for p in pongs:
        if r.search(p):
            chro  = r.search(p).groups(0)[0]
            start = r.search(p).groups(0)[1]
            end   = r.search(p).groups(0)[2]
            pong_dict['%s:%s' %(chro, start)]['direction'] = 'left'
            pong_dict['%s:%s' %(chro, end)]['direction']   = 'right'
            pong_dict['%s:%s' %(chro, start)]['loci']      = p
            pong_dict['%s:%s' %(chro, end)]['loci']        = p
            pong_dict['%s:%s' %(chro, start)]['loci_i']    = [chro, start, end]
            pong_dict['%s:%s' %(chro, end)]['loci_i']      = [chro, start, end]
    return pong_dict 

#parse ref or shared pong, only keep these known reference pong loci
#Japonica_fastq_RelocaTEi_Pong/ERS470370_RelocaTEi/repeat/results/ALL.all_nonref_insert.gff
#Chr1    not.give        RelocaTE_i      480846  480848  .       -       .       ID=repeat_Chr1_480846_480848;TSD=TTA;Note=Non-reference, not found in reference;
#Right_junction_reads:2;Left_junction_reads:1;Right_support_reads:0;Left_support_reads:0;
def ref_gff_parse(infile, pong_NB, pong_gff):
    data = defaultdict(lambda : str())
    r = re.compile(r'\=')
    #print infile
    #print ping_gff
    with open (infile, 'r') as filehd:
        for line in filehd: 
            line = line.rstrip()
            if len(line) > 2 and not line.startswith('#'): 
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
                        if r.search(attr):
                            idx, value = re.split(r'\=', attr)
                            temp[idx] = value
                        else:
                            idx, value = re.split(r'\:', attr)
                            temp[idx] = value
                Rjun    = temp['Right_junction_reads']
                Rsup    = temp['Right_support_reads']
                Ljun    = temp['Left_junction_reads']
                Lsup    = temp['Left_support_reads']
                if pong_NB.has_key('%s:%s' %(chro, start)) and pong_NB.has_key('%s:%s' %(chro, end)):
                    data[pong_NB['%s:%s' %(chro, start)]['loci']] = line
    ofile = open(pong_gff, 'w')
    for pong in sorted(data.keys()):
        print >> ofile, data[pong]
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
    
    flank_inf = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str())))

    #parse flankingReads.fq
    #Japonica_fastq_RelocaTEi_Pong/ERS470370_RelocaTEi/repeat/flanking_seq/ERR624495_1.te_repeat.flankingReads.fq
    flankingreads = glob.glob('%s/repeat/flanking_seq/*.te_repeat.flankingReads.fq' %(args.input))
    for fread in flankingreads:
        #print fread
        parse_flanking_fastq(fread, flank_inf)
 
    #print 'CK, ERR624495.1563, %s' %(flank_inf['ERR624495.1563'][str(1)]['pos'])   
    #for read in sorted(flank_inf.keys()):
    #    for lib in sorted(flank_inf[read].keys()):
    #        print '%s\t%s\t%s' %(read, lib, flank_inf[read][lib]['pos'])
 
    #parse blatout
    #Japonica_fastq_RelocaTEi_Pong/ERS470370_RelocaTEi/repeat/blat_output/ERR624495_1.te_repeat.blatout
    blatouts =glob.glob('%s/repeat/blat_output/*.te_repeat.blatout' %(args.input)) 
    for bout in blatouts:
        parse_blatout(bout, flank_inf)   

    #for read in sorted(flank_inf.keys()):
    #    for lib in sorted(flank_inf[read].keys()):
    #        print '%s\t%s\t%s' %(read, lib, flank_inf[read][lib]['pos'])

    #parse supporting_reads for each insertion
    support_inf = defaultdict(lambda : defaultdict(lambda : list()))
    #Japonica_fastq_RelocaTEi_Pong/ERS470370_RelocaTEi/repeat/results/Chr1.repeat.reads.list
    support_reads = glob.glob('%s/repeat/results/Chr*.repeat.reads.list' %(args.input)) 
    for sread in support_reads:
        parse_support_reads(sread, support_inf)

    pong_NB = setpong()
    #parse ref.gff and check with reference pong if exists
    ref_gff='%s/repeat/results/ALL.all_ref_insert.gff' %(args.input)
    ref_pong_gff = '%s/repeat/results/ALL.all_ref_insert.Pong.gff' %(args.input)
    if os.path.isfile(ref_gff):
        ref_gff_parse(ref_gff, pong_NB, ref_pong_gff) 

    #parse gff file and check junction reads and supporting reads if they cover mPing/Ping/Pong difference area
    all_gff='%s/repeat/results/ALL.all_nonref_insert.gff' %(args.input)
    if os.path.isfile(all_gff):
        gff_parse(all_gff, flank_inf, support_inf, pong_NB, ref_pong_gff)

if __name__ == '__main__':
    main()

