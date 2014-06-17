#!/usr/bin/env python
import Utils
from time import time
import argparse
from os import path
from os import makedirs
import shutil
import sys




###############
# Main menu
###############
main_parser = argparse.ArgumentParser(prog = 'MLST_Easy.py', description = 'Multi-function tool for NGS-MLST data analysis')
subparsers = main_parser.add_subparsers()

arg_all = subparsers.add_parser('all',help='run the whole pipeline for the data')
arg_all.add_argument('-c',dest='Config_File',help='config file for the pipeline, default: config.ini')
arg_all.set_defaults(func=all)

arg_align = subparsers.add_parser('align',help='search the barcodes and primers in the reads')
arg_align.set_defaults(func=align)
arg_align.add_argument('-c',dest='Config_File',help='config file for the pipeline, default: config.ini')

arg_cons = subparsers.add_parser('cons',help='generate the consensus sequences for the aligned reads')
arg_cons.add_argument('-c',dest='Config_File',help='config file for the pipeline, default: config.ini')
arg_cons.set_defaults(func=cons)

arg_unmapped = subparsers.add_parser('unmap',help="dump the reads can't be aligned to barcodes or primers")
arg_unmapped.add_argument('-c',dest='Config_File',help='config file for the pipeline, default: config.ini')
arg_unmapped.set_defaults(func=unmap)

arg_merge = subparsers.add_parser('merge',help='merge the aligned reads from different SMRT cells')
arg_merge.add_argument('-c',dest='Config_Files',required=True,nargs='+',metavar='config_file1 config_file2 ...', help='config files for the merge libraries')
arg_merge.add_argument('-o',dest='Output_Folder',required=True,help='output folder of the merged results')
arg_merge.set_defaults(func=merge)

arg_het = subparsers.add_parser()

arg_extract = subparsers.add_parser('exat',help='extract reads by strain and locus from aligned results')
arg_extract.add_argument('-f',dest='Input_Folder',required=True,help='folder contains aligned results')
arg_extract.add_argument('-s',dest='Strain',required=True,help='strain name')
arg_extract.add_argument('-l',dest='Locus',required=True,help='locus name')
arg_extract.add_argument('-m',dest='Out_Format',help='output format of the reads <fasta|fastq>,default:fasta')
arg_extract.add_argument('-t',dest='Out_Mode',type=int,help='output mode of the reads <0:trimed by primer|1:raw>,default:0')
arg_extract.set_defaults(func=exat)



arg_stats = subparsers.add_parser('stats',help='generate statistical information for the analysis')
arg_stats.add_argument('-c',dest='Config_File',help='config file for the pipeline, default: config.ini')
arg_stats.set_defaults(func=stats)

if(len(sys.argv)==1):
    main_parser.print_help()
    sys.exit(1)
args = main_parser.parse_args()
args.func(args)