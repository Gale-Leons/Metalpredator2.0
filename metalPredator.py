#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob

import subprocess
#import commands

import multiprocessing as mp

from lib.managers import Fragments_Manager, MFS_Manager, Pfam_Manager
from lib.utilities import read_fasta_as_dict, write_fasta_dict_to_file
from lib.utilities_m2 import parseCommandLineOptions
import time


def worker(method, userdir, metal, user_sequences, db_path):
    
    time0=time.time()

    outfile_summary=os.path.join(userdir, '%s_summary_annotation.csv' % method)
    outfile_detailed=os.path.join(userdir, '%s_detailed_annotation.csv' % method)
    
    if method=='fragments':
        rmg = Fragments_Manager(float(1e-3))
    elif method=='mfs':
        rmg = MFS_Manager(float(1e-3))
    elif method=='pfam':
        rmg = Pfam_Manager(float(1e-5))
    elif method=='pfam_no_structure':
        rmg = Pfam_Manager(float(1e-7))
    
    rmg.set_metal(metal)
    if db_path:
       rmg.set_db_path(db_path)
    black_exist = rmg.set_blacklist()
    rmg.set_database(method)
    
    rmg.load_input_sequences(user_sequences)
    file_not_empty = rmg.create_library_of_profiles()
    if file_not_empty:
       detected_sequences = rmg.get_detected_sequences(userdir) 

       rmg.filter_matched_domains(detected_sequences)
       if black_exist:
          rmg.filter_blacklisted_domains(userdir, detected_sequences)
    else:
       detected_sequences = []

    for det_seq in detected_sequences:
        for i in range(0, len(det_seq.predicted_ligands)):
            for j in range(0, len(det_seq.predicted_ligands[i])):
                det_seq.predicted_ligands[i][j] = sorted(det_seq.predicted_ligands[i][j])

    rmg.get_sequences_annotations(detected_sequences, outfile_summary, outfile_detailed)
    
    print('search', job, ':', time.time()-time0)
    
def worker_summary(userdir):
    
    summary = {}
    methods = ['pfam', 'mfs', 'fragments', 'pfam_no_structure']
   
    for method in methods:
        
        filepath=os.path.join(userdir, '%s_detailed_annotation.csv' % method)
        with open(filepath, 'r') as fl:
            results=fl.readlines()
        
        prima = True
        for line in results:
            
            if prima:
                prima = False
                continue
            
            fields=line.strip('\n').strip().split(',')

            sequence_id=fields[0]
            
            sequence = next((s for s in summary.keys() if s==sequence_id), None)
            
            if method=='fragments':
                e_value=float(fields[8])
                ligands=fields[9].split('-')
                if ligands == ['']:
                   if metal!='molybdenum':
                      continue
                   else:
                      ligands = ['No Ligands']
 
            
            elif method=='mfs':
                e_value=float(fields[7])
                ligands=fields[8].split('-')
                if ligands == ['']:
                   if metal!='molybdenum':
                      continue
                   else:
                      ligands = ['No Ligands']
                
            elif method=='pfam':
                e_value=float(fields[5])
                ligands=fields[6].split('-')
                
            elif method=='pfam_no_structure':
                e_value=float(fields[5])
                ligands = [fields[2]]
                
            if not sequence:
                summary[sequence_id]={method:[e_value, ligands]}
            
            elif method not in summary[sequence_id].keys():
                summary[sequence_id].update({method:[e_value, ligands]})
            
            else:
                min_e_value, all_ligands = summary[sequence_id][method]
                if e_value < min_e_value:
                    summary[sequence_id][method][0] = e_value
                    if method == 'pfam_no_structure' or ligands[0]=='No Ligands':
                        summary[sequence_id][method][1] = ligands
                        
                if method != 'pfam_no_structure' and ligands[0]!='No Ligands':
                    for lig in ligands:
                        if lig not in summary[sequence_id][method][1]:
                            summary[sequence_id][method][1].append(lig)

    summary_sorted_keys = sorted(summary.keys(), key=lambda t: len(summary[t].keys()), reverse=True)
    summary_sorted_ligands = sorted(summary_sorted_keys, key=lambda t: len(filter(None, map(lambda x: x[1], [v for k, v in summary[t].items() if k in summary[t]]))), reverse=True)
    
    summary_data = []
    for sequence_id in summary_sorted_ligands:
        line=sequence_id
        bad = False
        for method in methods:
            if method in summary[sequence_id].keys():
                ligands = summary[sequence_id][method][1]
                if ligands == ['']:
                   if metal!='molybdenum':
                      ligands_str = ''
                   else:
                      ligands_str = ligands[0]
                elif method == 'pfam_no_structure' or "".join(list(set(ligands)))=='No Ligands':
                    ligands_str = ligands[0]
                else:
                    if 'No Ligands' in ligands:
                        summary[sequence_id][method][1] = [x for x in ligands if x != 'No Ligands']                   
                    try:
                        ligands = sorted(summary[sequence_id][method][1], key=lambda x: int(x[1:]))
                        ligands_str = '-'.join(ligands)
                    except:
                        bad = True     
                        continue
                line += ',%s'%(ligands_str)
            else:
                line += ',-'
        summary_data.append(line+'\n')
        
    head_summ = "Sequence ID,Domain search with pattern filter,MFS search,Fragment search,Domain search without pattern filter\n"
    summary_file = os.path.join(userdir, 'methods_summary.csv')
    with open(summary_file, 'w') as fl:
        fl.write(head_summ)
        fl.writelines(summary_data)
        
        
def get_input_path(prog_name):
    
    resp = input("\nFull path to '%s' [leave it blank to exit]:"%prog_name)
    if resp == '':
        sys.exit()
    else:
        if os.path.exists(resp):
            return resp
        else:
            print("\nPath not correct.\n")
            count = 0
            while count<3:
                try_again = input("Try again?[y/n]:")
                if try_again == "y":
                    get_input_path(prog_name)
                    count == 4
                elif try_again == "n":
                    sys.exit()
                else:
                    count += 1
                    if count<3:
                        print("\nPress 'y' to give path or 'n' to terminate program.")
                    else:
                        sys.exit()
                        
                
                
def file_path_correct(path_hmm):
    
    if not path_hmm:
        this_file_path = os.path.realpath(__file__)
        local_path = "/".join(this_file_path.split("/")[0:-1])
        print("LOCAL PATH: %s"%local_path)
        full_path = os.path.join(local_path, "data/local_path_hmm.txt")
        if not os.path.exists(full_path):
            return False
        
        
        file_path = open(full_path, "r")
        line_path = file_path.readlines()
        file_path.close()
        
        path_hmm = line_path[0].strip()
        
        
    file_control_text = "%s -h" % path_hmm
    resp_control = subprocess.getstatusoutput(file_control_text)
    if resp_control[0] != 0:
        return False
    else:
        return True


        
        
def set_path_hmm():
        
        if os.name == "nt":
            where_is = 'where hmmscan'
        else:
            where_is = 'which hmmscan'

        path_hmm = subprocess.getstatusoutput(where_is)

        if (path_hmm[0] != 0) or (path_hmm[1]==''):
            print("\n---------WARNING----------")
            print("\nIt is not possible found 'hmmscan' software.\nMetalPredator needs 'hmmer' software installed.\nTo download it: https://www.ebi.ac.uk/Tools/hmmer/")
            print("\nIf you are sure that it is correctly installed you can type the path now.")
            path_hmm_tool = get_input_path('hmmscan')
            while not file_path_correct(path_hmm_tool):
                path_hmm_tool = get_input_path('hmmscan')
        else:
            path_hmm_tool = path_hmm[1]
            
        file_path = open("data/local_path_hmm.txt", "w")
        file_path.write(path_hmm_tool)
        file_path.close()
        
        return True
    

def getPathForTempFasta(user_path_seqs):
    
    fasta_name = user_path_seqs.split("/")[-1]
    this_file_path = os.path.realpath(__file__)
    seq_tmp_path = "%s/data/%s"%("/".join(this_file_path.split("/")[0:-1]), fasta_name)

    return seq_tmp_path 
    
def getPathForTempFasta2(user_path_seqs, now_out):

    fasta_name = user_path_seqs.split("/")[-1]
    seq_tmp_path_dir = os.path.join(now_out, "tmp_seq")
    os.mkdir(seq_tmp_path_dir)
    seq_tmp_path = os.path.join(seq_tmp_path_dir, fasta_name)

    return seq_tmp_path 



if __name__ == '__main__':
    
    try:
        argo = sys.argv[1]
    except:
        print("\nInput parameter for input file is mandatory.\nFor help use -h or --help\n")
        sys.exit()

    try:
        if sys.argv[1].strip()!='-h' and sys.argv[1].strip()!='--help':
           print('\nPre-processing the input data. Please, wait...')

        inputType, metal, inputList, pathToOutput, db_path = parseCommandLineOptions(sys.argv)
        print(file_path_correct)
        if not file_path_correct(False):
            set_path_hmm()

        print('Ready to go.\n')

        this_file_path = os.path.realpath(__file__)
        
        seq_list = []
        out_list = []
        if inputType == 'list':
            for fastafile in os.listdir(inputList):
                path_input = os.path.join(inputList, fastafile)
                pathToOutput = os.path.join(pathToOutput, fastafile)
                seq_list.append(path_input)
                out_list.append(out_list)
        else:
            seq_list.append(inputList)
            out_list.append(pathToOutput)
            
        
        for i in range(0, len(seq_list)):
            
            print('\nWorking %s ...\n'%seq_list[i])
            
            user_sequences = seq_list[i]
            pathToOutput = out_list[i]
            
            if not os.path.exists(pathToOutput):
                    os.mkdir(pathToOutput)
                
            root = pathToOutput
               
            lock = mp.Lock()
    
            try:
                
                lock.acquire()
                sequences = read_fasta_as_dict(user_sequences, True)
                seq_tmp_path = getPathForTempFasta(user_sequences)
                #seq_tmp_path = getPathForTempFasta2(user_sequences, pathToOutput) 
                write_fasta_dict_to_file(seq_tmp_path, sequences)
                lock.release()
        
                processes = []
                jobs = ['pfam', 'mfs', 'fragments', 'pfam_no_structure']
               
                for job in jobs:
        
                    p = mp.Process(target=worker, args=(job, pathToOutput, metal, seq_tmp_path, db_path))
                    processes.append(p)
                    p.start()
            
                for p in processes:
                    p.join()
                
                worker_summary(pathToOutput)
                            
            except Exception as e:
                error_message = os.path.join(pathToOutput, 'error')
                with open(error_message, 'w') as fl:
                    fl.write(str(e)+'\n')
                print(e)
                continue
        
            pathToOutput = root
                    
            print('Results are stored here: %s\n' % pathToOutput)
            os.remove(seq_tmp_path)
    
    except SystemExit:
        pass
