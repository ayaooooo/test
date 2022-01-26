#!/usr/bin/env python3

import os
import re
import sys
import copy
import argparse

def main():
    BIN = os.path.dirname(os.path.abspath(__file__))

#---------------------------------------parameters start------------------------------------------
    parser = argparse.ArgumentParser(description = "collating the data from the pipeline raw results (include sample_4M_result, sample_mosaic_result, sample_qc_Karyotype_result and sample_hospital_result")
    parser.add_argument('-i', '--input', help = 'the path of raw results', required = True)
    parser.add_argument('-o', '--outdir', help = 'the out path of new results', required = True)
    args = parser.parse_args()
#---------------------------------------parameters end--------------------------------------------
    karyotype_file  = os.path.abspath(args.input) + '/sample_qc_Karyotype_result'  
    hospital_file   = os.path.abspath(args.input) + '/sample_hospital_result'
    qc_file = os.sep.join([os.path.abspath(args.outdir), 'sample_qc_halos_result.txt'])
    result_file1 = os.sep.join([os.path.abspath(args.outdir), 'sample_karyotype_4M_halos_result.txt'])
    result_file2 = os.sep.join([os.path.abspath(args.outdir), 'sample_karyotype_16M_halos_result.txt'])
    result_file3 = os.sep.join([os.path.abspath(args.outdir), 'sample_karyotype_5chr_halos_result.txt'])
