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
import new_umrplots as up

def main(argv):

    parser = argparse.ArgumentParser(usage="python3 raptor.py [-i <input>] [-b] [-hex] [-hum] [-u]")
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
    parser.add_argument('-u', '--umrplots',help='Generate analsysis plots for UMRs', action='store_true')
    # parser.add_argument('-f', '--fly', help='Generate a csv file of comprehensive information of UMR sequences for Fly sample ',
    #                     action='store_true')
    args = parser.parse_args()

    ##Take input sam file and generate UMR text file, UMR fastq file and UMR bed file
    if args.input:
        input_filename = args.input.name
        print("RAPTOR (v 1.0 ) - Comprehensive Analysis of 3' UnMappable Regions")
        print("Input file is", args.input.name)
        print("Process running...")
        fr.extract(input_filename)
        print("Process Completed.")
        print(input_filename[:-4] + "_3_tail_seqs.txt contains 3' UMR sequences.")
        print(input_filename[:-4] + "_3_tails.txt contains UMR sequences and associated information.")
        print(input_filename[:-4] + "_3_UMR.fastq contains UMR sequences in FASTQ format for analysis plots.")
        print(input_filename[:-4] + "_3_UMR.bed contains UMR sequences in bed format.")
    ##Take input UMR text file and generate a csv file of Nucleotide composition
    if args.bases:
        input_filename = args.input.name
        umr_textfile = input_filename[:-4] + '_3_tails.txt'
        bc.base_comp(umr_textfile)
        print("Base Composition analysis completed.")
        print(input_filename[:-4] + "_Nucleotidecomposition.csv contains composition of bases across UMR regions.")
    ##Take input UMR text file and generate a plot of top hexamers and hexamers percentage csv file
    if args.hexamers:
        input_filename = args.input.name
        umr_textfile = input_filename[:-4] + '_3_tails.txt'
        hc.plothexamers(umr_textfile)
        print("Hexamer composition analysis completed.")
        print(input_filename[:-4] + "_hexamer_composition.csv contains composition of Hexamers in UMR regions.")
        print(input_filename[:-4] + "_top_hexamers.pdf contains distribution of Hexamers in UMR regions.")
    ## Take input bed file and intersect with human genome and report overlap results as csv file
    if args.human:
        input_filename = args.input.name
        inputbedfile = input_filename[:-4] + '_3_Umr.bed'
        umr_textfile = input_filename[:-4] + '_3_tails.txt'
        umr_bedfile=pybedtools.BedTool(inputbedfile)
        genefile=pybedtools.BedTool('human_genes.bed')
        intersectfile=umr_bedfile.intersect(genefile, wao=True, f=0.70)
        intersect_df = pd.read_table(intersectfile.fn, names = ['ReadChrom','ReadStart','ReadEnd','ReadID','Strand','TransChrom','TransStart','TransEnd','TID','GeneID','GeneName','TLength','Overlap'])
        intersect_df = (intersect_df.sort_values(['ReadID','Overlap'], ascending=[True, False]).drop_duplicates(['ReadID']).reset_index(drop=True))
        intersect_df = intersect_df[intersect_df['TLength'] != '.']
        intersect_df['TLength'] = intersect_df['TLength'].astype(float)
        intersect_df['Identity'] = intersect_df['Overlap']/intersect_df['TLength']
        intersect_df = intersect_df[intersect_df['Identity'] >= 0.70]
        intersect_df = intersect_df.drop(['ReadChrom','ReadStart','ReadEnd','Strand','Identity'], axis=1)
        umr_df = pd.read_table(umr_textfile, names=['ReadID','UMRLength','UMRSequence','ReadChrom','ReadStart','ReadEnd','Strand'],sep=' ')
        result_df = umr_df.merge(intersect_df, on='ReadID')
        # print(result_df)
        result_df.to_csv(input_filename[:-4] + '_3\'Umr_analysis.csv',mode='a')
    ## Take input bed file and intersect with human genome and report overlap results as csv file
    # if args.mouse:
    #     file1 = args.input.name
    #     file2= file1[:-4] + '_3_UMR.bed'
    #     filename = file1[:-4]
    #     a=pybedtools.BedTool(file2)
    #     b=pybedtools.BedTool('mouse_genes.bed')
    #     c=a.intersect(b, wao=True, f=0.70)
    #     data1 = pd.read_table(c.fn, names = ['ReadChrom','ReadStart','ReadEnd','ReadID','Strand','TransChrom','TransStart','TransEnd','TID','GeneID','GeneName','TLength','Overlap'])
    #     data1 = (data1.sort_values(['ReadID','Overlap'], ascending=[True, False]).drop_duplicates(['ReadID']).reset_index(drop=True))
    #     data1 = data1[data1['TLength'] != '.']
    #     data1['TLength'] = data1['TLength'].astype(float)
    #     data1['Identity'] = data1['Overlap']/data1['TLength']
    #     data1 = data1[data1['Identity'] >= 0.70]
    #     data1 = data1.drop(['ReadChrom','ReadStart','ReadEnd','Strand','Identity'], axis=1)
    #     data2 = pd.read_table('3_tails.txt', names=['ReadID','UMRLength','UMRSequence','ReadChrom','ReadStart','ReadEnd','Strand'],sep=' ')
    #     data3 = data2.merge(data1, on='ReadID')
    #     print(data3)
    # if args.yeast:
    #     a=pybedtools.BedTool('3_UMR.bed')
    #     b=pybedtools.BedTool('yeast_genes.bed')
    #     c=a.intersect(b, wao=True, f=0.70)
    #     data1 = pd.read_table(c.fn, names = ['ReadChrom','ReadStart','ReadEnd','ReadID','Strand','TransChrom','TransStart','TransEnd','TID','GeneID','GeneName','TLength','Overlap'])
    #     data1 = (data1.sort_values(['ReadID','Overlap'], ascending=[True, False]).drop_duplicates(['ReadID']).reset_index(drop=True))
    #     data1 = data1.drop(['ReadChrom','ReadStart','ReadEnd','Strand'], axis=1)
    #     data2 = pd.read_table('3_tails.txt', names=['ReadID','UMRLength','UMRSequence','ReadChrom','ReadStart','ReadEnd','Strand'],sep=' ')
    #     data3 = data2.merge(data1, on='ReadID')
    if args.umrplots: 
        file1 = args.input.name
        file2= file1[:-4] + '_3_UMR.fastq'
        # print(file2)
        up.plotter(file2)
        print('3\'UMR Analysis plots generated')
        print(file2[:-6] +'.pdf contains 3\' UMR Analysis plots')

if __name__ == "__main__":
    main(sys.argv[1:])

