#!/usr/bin/env python3
# Author: Mengzhu Liu
# Date: 2018.5.28

import os, sys
import logging
import multiprocessing
import pandas as pd
import numpy as np
import argparse

def parse_args():
    
    parser = argparse.ArgumentParser(description='This is a program to call mutations in VDJ sequences (by Mengzhu)')
    parser.add_argument("-i", dest = "fafile", type = str, required = True,
                              help = "consensus fa file" )
    parser.add_argument("-o", dest = "outputdir", type = str, required = True,
                              help = "output directory" )
                              
    args = parser.parse_args()
    
    return args
    
def Consensus(args):
    
    filename = args.fafile
    workdir = args.outputdir
    basename = filename.rpartition('/')[2].rpartition('_')[0]
    fa_file = pd.read_table(filename, sep='\t')
    
    #igblast
    germline_db_V = "~/database/VDJ-database/imgt/mouse_IGHV_imgt"
    germline_db_J = "~/database/VDJ-database/imgt/mouse_IGHJ_imgt"
    germline_db_D = "~/database/VDJ-database/imgt/mouse_IGHD_imgt"
    germline_db_VL = "~/database/VDJ-database/imgt/mouse_IGLV_imgt"
    germline_db_JL = "~/database/VDJ-database/imgt/mouse_IGLJ_imgt"
    germline_db_VK = "~/database/VDJ-database/imgt/mouse_IGKV_imgt"
    germline_db_JK = "~/database/VDJ-database/imgt/mouse_IGKJ_imgt"
    
    organism = "mouse"
    aux_file = "~/software/ncbi-igblast-1.7.0/optional_file/mouse_gl.aux"

    print("\nStart runing igblastn for joined reads!")

    #H
    cmd = "igblastn -germline_db_V {0} \
           -germline_db_J {1} \
           -germline_db_D {2} \
           -organism {3} \
           -domain_system imgt \
           -auxiliary_data {4} \
           -ig_seqtype Ig \
           -outfmt '7 std qseq sseq btop' \
           -show_translation \
           -query {5} \
           -out {6}".format (\
           germline_db_V, germline_db_J, germline_db_D, organism, aux_file, workdir+"/"+filename, workdir+"/"+basename+".H.out")
    print(cmd)
    os.system(cmd)
    

    #makeDb
    germline_db = "~/database/VDJ-database/imgt/"
    cmd = "MakeDb.py igblast -i {0} -s {1} -r {2} --cdr3 --scores --partial".format(workdir+"/"+basename+".H.out", workdir+"/"+filename, germline_db)
    os.system(cmd)
    
    data = pd.read_table(workdir+"/"+basename+".H_db-pass.tab", sep='\t')
    
    data.to_csv(workdir+'/'+basename+'_consensus.tab',header=True,sep='\t',index=False,
                columns = ['SEQUENCE_ID','SEQUENCE_INPUT','FUNCTIONAL','IN_FRAME','V_CALL',
                'D_CALL','J_CALL','V_SEQ_LENGTH','D_SEQ_LENGTH','J_SEQ_LENGTH','V_SCORE',
                'V_IDENTITY','V_BTOP','J_BTOP','CDR3_IGBLAST_NT','CDR3_IGBLAST_AA'])

def main():
    args = parse_args()
    Consensus(args)
main()