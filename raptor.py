import sys, getopt
import os
import re
import argparse
import pandas as pd
from argparse import ArgumentParser
import pybedtools

import extract as fr
import base_composition as bc
import hexamer_count as hc

def main(argv):

    parser = argparse.ArgumentParser(usage="python3 raptor.py [-i <input>] [-b] [-hex] [-hum]")
    parser.add_argument('-i','--input',help="Input Sam file ", type=argparse.FileType('r'))
    parser.add_argument('-b', '--bases',help="Generate a csv file of nucleotide composition in UMR sequences",
                        action='store_true' )
    parser.add_argument('-hex', '--hexamers', help='Generate a csv file and plots of Hexamer frequencies in UMR sequences ',
                        action='store_true')
    parser.add_argument('-hum', '--human', help='Generate a csv file containing comprehensive information of UMR sequences for Human sample',
                        action='store_true')
    parser.add_argument('-mo', '--mouse', help='Generate a csv file of comprehensive information of UMR sequences for Mouse sample',
                        action='store_true')
    parser.add_argument('-y', '--yeast',
                        help='Generate a csv file of comprehensive information of UMR sequences for Saccharomyces Cerivisiae ',
                        action='store_true')
    # parser.add_argument('-f', '--fly', help='Generate a csv file of comprehensive information of UMR sequences for Fly sample ',
    #                     action='store_true')
    args = parser.parse_args()

    if args.input:
        file1 = args.input.name
        print("RAPTOR (v 1.0 ) - Comprehensive Analysis of 3' UnMappable Regions")
        print("Input file is", args.input.name)
        print("Process running...")
        fr.extract(file1)
        print("Process Completed.")
        print("3_tail_seqs.txt contains 3' UMR sequences.")
        print("3_tails.txt contains UMR sequences and associated information.")
        print("3_UMR.fastq contains UMR sequences in FASTQ format for analysis plots.")
    if args.bases:
        tail_data = pd.read_csv("3_tails.txt", sep=" ", error_bad_lines=False,
                                names=['ReadID', 'Length', 'Sequence', 'Chrom', 'Start', 'End', 'Strand'])
        bc.base_comp(tail_data)
        print("Base Composition analysis completed.")
        print("base_composition.csv contains composition of bases across UMR regions.")
    if args.hexamers:
        tail_data = pd.read_csv("3_tails.txt", sep=" ", error_bad_lines=False,
                                names=['ReadID', 'Length', 'Sequence', 'Chrom', 'Start', 'End', 'Strand'])
        hc.hexamerplot(tail_data)
        print("Hexamer composition analysis completed.")
        print("hexamer_composition.csv contains composition of Hexamers in UMR regions.")
        print("hexamer_distribution.jpg contains distribution of Hexamers in UMR regions.")
    if args.human:
        a=pybedtools.BedTool('3_UMR.bed')
        b=pybedtools.BedTool('human_genes.bed')
        c=a.intersect(b, wao=True, f=0.70)
        data1 = pd.read_table(c.fn, names = ['ReadChrom','ReadStart','ReadEnd','ReadID','Strand','TransChrom','TransStart','TransEnd','TID','GeneID','GeneName','TLength','Overlap'])
        data1 = (data1.sort_values(['ReadID','Overlap'], ascending=[True, False]).drop_duplicates(['ReadID']).reset_index(drop=True))
        data1 = data1.drop(['ReadChrom','ReadStart','ReadEnd','Strand'], axis=1)
        data2 = pd.read_table('3_tails.txt', names=['ReadID','UMRLength','UMRSequence','ReadChrom','ReadStart','ReadEnd','Strand'],sep=' ')
        data3 = data2.merge(data1, on='ReadID')
        # print(data3)
        data3.to_csv("3_Umr_Analysis.csv",mode='a')
        print("3_Umr_Analysis.csv contains Transcript and Gene information")
    if args.mouse:
        a=pybedtools.BedTool('3_UMR.bed')
        b=pybedtools.BedTool('mouse_genes.bed')
        c=a.intersect(b, wao=True, f=0.70)
        data1 = pd.read_table(c.fn, names = ['ReadChrom','ReadStart','ReadEnd','ReadID','Strand','TransChrom','TransStart','TransEnd','TID','GeneID','GeneName','TLength','Overlap'])
        data1 = (data1.sort_values(['ReadID','Overlap'], ascending=[True, False]).drop_duplicates(['ReadID']).reset_index(drop=True))
        data1 = data1.drop(['ReadChrom','ReadStart','ReadEnd','Strand'], axis=1)
        data2 = pd.read_table('3_tails.txt', names=['ReadID','UMRLength','UMRSequence','ReadChrom','ReadStart','ReadEnd','Strand'],sep=' ')
        data3 = data2.merge(data1, on='ReadID')
        # print(data3)
        data3.to_csv("3_Umr_Analysis.csv", mode='a')
        print("3_Umr_Analysis.csv contains Transcript and Gene information")
    if args.yeast:
        a=pybedtools.BedTool('3_UMR.bed')
        b=pybedtools.BedTool('yeast_genes.bed')
        c=a.intersect(b, wao=True, f=0.70)
        data1 = pd.read_table(c.fn, names = ['ReadChrom','ReadStart','ReadEnd','ReadID','Strand','TransChrom','TransStart','TransEnd','TID','GeneID','GeneName','TLength','Overlap'])
        data1 = (data1.sort_values(['ReadID','Overlap'], ascending=[True, False]).drop_duplicates(['ReadID']).reset_index(drop=True))
        data1 = data1.drop(['ReadChrom','ReadStart','ReadEnd','Strand'], axis=1)
        data2 = pd.read_table('3_tails.txt', names=['ReadID','UMRLength','UMRSequence','ReadChrom','ReadStart','ReadEnd','Strand'],sep=' ')
        data3 = data2.merge(data1, on='ReadID')
        data3.to_csv("3_Umr_Analysis.csv", mode='a')
        print("3_Umr_Analysis.csv contains Transcript and Gene information")

if __name__ == "__main__":
    main(sys.argv[1:])
