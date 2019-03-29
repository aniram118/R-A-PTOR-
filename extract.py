#!/usr/bin/python

import re
import sys
import getopt

import sys, getopt

from cigar import cigar_parse


def extract(argv):
    inputfile = argv
    outputfile = ''
    filename=inputfile[:-4]

    ## Create placeholder for input and output files

    file = open(inputfile,'r')
    file1 = open(inputfile, 'r')
    # file2 = open('5_tails.txt', 'a')
    file3 = open('%s_3_tails.txt'%filename, 'a')
    # f3 = open('5_tail_seqs.txt', 'a')
    f4 = open('%s_3_tail_seqs.txt'%filename, 'a')
    file5 = open('%s_3_UMR.fastq'%filename, 'a')
    file6 = open('%s_3_UMR.bed'%filename, 'a')

    sequences = []
    strings = []
    raw_seq = []
    cigar = []
    reads = []
    quality = []

    ### File and Cigar parsing

    for a, i in zip(file1, file):
        i1 = i.split()
        a1 = a.split()
        istring = i.split('MD:Z:')
        try:
            if not i.startswith('@'):
                sequences.append(a1[9])
                strings.append(istring[1])
                raw_seq.append(a1[9])
                cigar.append(a1[5])
                start_g = int(''.join(a1[3]))  #### added
                cigar_s = ''.join(a1[5])  #### added
                chr_no = ''.join(a1[2])

                # print(a1[9])
                def my_split(s):
                    return filter(None, re.split(r'(\d+)', s))

                s = list(my_split(a1[5]))
                strt = 0
                end = 0
                seql = ''.join(a1[9])
                seql.replace('T','U')
                qual = ''.join(a1[10])

                # print(s)
                ### Identify positive and negative strand and extract tails
                if int(i1[1]) == 16:
                    if s[1] == 'S':
                        #  print('Start'+' '+str(s[0]))
                        strt = (int(s[0]))
                        # print(strt)
                        # print(seql[0:][:strt])
                        if strt > 6:
                            reads.append(i1[0])
                            quality.append(qual[0:][:strt])
                            m = qual[0:][:strt]
                            # print(seql[0:][:strt] + ' ' + str(strt) + ' ' + 'start')
                            ##3' UMR info file
                            file3.write(i1[0] + ' ' + str(strt) + ' ' + seql[0:][:strt] + ' ' + chr_no + ' ' + str(start_g) + ' ' + str(cigar_parse(cigar_s, start_g)) + ' ' + '-' + '\n')
                            ##3' UMR Sequences
                            f4.write(seql[0:][:strt] + '\n')
                            ##3' UMR fastq file
                            file5.write('@' + i1[0] + ' ' + 'XXXX' + '\n' + seql[0:][:strt] + '\n' + '+' + '\n' + qual[0:][:strt] + '\n')
                            ##3' UMR Bed file
                            file6.write(chr_no + '\t' + str(start_g) + '\t' + str(cigar_parse(cigar_s,start_g)) + '\t' + i1[0] + '\t' + '-' + '\n')
                    # if s[len(s) - 1] == 'S':
                        #  print('End'+' '+str(s[len(s)-2])+' '+str(s[len(s)-1]))
                        # end = (int(s[len(s) - 2]))
                        # print(seql)
                        # print("end "+str(end))
                        # print(seql[-end:])
                        # if end > 6:
                            # print(seql[-end:] + ' ' + str(end) + ' ' + 'end')
                            # file2.write(i1[0] + ' ' + str(end) + ' ' + seql[-end:] + '\n')
                            # f3.write(seql[-end:] + '\n')
                elif int(i1[1]) == 0:
                    # if s[1] == 'S':
                        #  print('Start'+' '+str(s[0]))
                        # strt = (int(s[0]))
                        # print(strt)
                        # print(seql[0:][:strt])
                        # if strt > 6:
                            # print(seql[0:][:strt] + ' ' + str(strt) + ' ' + 'start')
                            # file2.write(i1[0] + ' ' + str(strt) + ' ' + seql[0:][:strt] + '\n')
                            # f3.write(seql[0:][:strt] + '\n')
                    if s[len(s) - 1] == 'S':
                        #  print('End'+' '+str(s[len(s)-2])+' '+str(s[len(s)-1]))
                        end = (int(s[len(s) - 2]))
                        # print(seql)
                        # print("end "+str(end))
                        # print(seql[-end:])
                        if end > 6:
                            # print(seql[-end:] + ' ' + str(end) + ' ' + 'end')
                            reads.append(i1[0])
                            quality.append(qual[-end:])
                            m = qual[-end:]
                            file3.write(i1[0] + ' ' + str(end) + ' ' + seql[-end:] + ' ' + chr_no + ' ' + str(start_g) + ' ' + str(cigar_parse(cigar_s, start_g)) + ' ' + '+' + '\n')
                            f4.write(seql[-end:] + '\n')
                            file5.write('@' + i1[0] + ' ' + 'runid=XXXX' + '\n' + seql[-end:] + '\n' + '+' + '\n' + qual[-end:] + '\n')
                            file6.write( chr_no + '\t' + str(start_g) + '\t' + str(cigar_parse(cigar_s, start_g)) + '\t' + i1[0] + '\t' + '+' + '\n')
        except IndexError:
            continue

    # print(''.join(strt))

# data = ("mapped_HepG2.sam")
# extract(data)
# extract(sys.argv[1:])

