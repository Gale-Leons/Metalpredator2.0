#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

from lib.hmm_data_structures import Matched_Sequence, Domain

def read_fasta_as_dict(filepath, my_head=False):
    
    with open(filepath, 'r') as fl:
        filecontent = fl.readlines()
        
    sequences_dict = {}
    
    i=0
    j=0
    sequence = ''
    while i < len(filecontent):
        
        line = filecontent[i]
        
        if line.startswith(">"):
        
            header_tmp = line.strip().strip('\n').strip(">")
            if my_head:
               header = "%s@%s"%(str(j), header_tmp)
               j+=1
            else:
               header = header_tmp  
            sequence = ''
            i += 1
        else:
            while not line.startswith(">"):
                
                sequence += line.strip('\n')
                
                i += 1
                if i == len(filecontent):
                    break
                
                line = filecontent[i]
        
            sequences_dict[header] = sequence
    
    sequences_dict[header] = sequence
    
    return sequences_dict

def write_fasta_dict_to_file(filepath, fastadict):
    
    lines = []
    for header, sequence in fastadict.items():
        if not header.startswith(">"):
            header = ">"+header
        lines.append(header+'\n')
        lines.append(sequence+'\n')

    with open(filepath, "w") as fl:
        fl.writelines(lines)
        
def yield_interval(data):
    
    start = None
    end = None
    for entry in sorted(data):
        if start is None:
            start = entry[0]
        if end is None or entry[0] < end:
            end = entry[1]
        else:
            yield [start, end]
            start, end = entry
    yield [start, end]

def merge_ranges(ranges):
    
    """
    Merge overlapping and adjacent ranges and yield the merged ranges
    in order. The argument must be an iterable of pairs (start, stop).
    """
    
    ranges = iter(sorted(ranges))
    current_start, current_stop = next(ranges)
    for start, stop in ranges:
        if start > current_stop:
            # Gap between segments: output current segment and start a new one.
            yield current_start, current_stop
            current_start, current_stop = start, stop
        else:
            # Segments adjacent or overlapping: merge.
            current_stop = max(current_stop, stop)
    yield current_start, current_stop
    
def get_number_of_gaps(sequence, ind, gap_type):
    
    i=0
    j=0
    gaps=0
    while j <= ind:
        if sequence[i] == gap_type:
            gaps += 1
        else:
            j += 1
        i += 1
        
    return gaps

def parse_hmmscan_output(filename):

    with open(filename,"r") as fl:
        lines = fl.readlines()
    lines = map(lambda x: x.strip('\n'), lines)
     
    run = False
    annotation = False
    alignments = False
    
    match = None
    
    i = 1
    while i < len(lines):
        
        line = lines[i]
        
        if "Query" in line:
            match_name = line.split()[1]
            try:
                description = ' '.join(line.split()[2:])
            except (Exception, e):
                pass
            i+=1
            continue

        if 'Domain annotation for each' in line:
            
            i += 1
            run = True
            
            match = Matched_Sequence()
            match.name = match_name
            match.description = description

        found_domains = []
        
        while run:
            
            line = lines[i]
            
            if not line:
                i += 1
                continue
            
            if "No targets detected" in line:
                return []
                
             # Now we stop reading the file
            if line.startswith("//") or ("Internal pipeline statistics summary" in line):
                return match
            
            if line.startswith(">>"):
                
                if "No individual domain" in lines[i+1]:
                    i += 3
                    continue
                
                try:
                    query_hmm = line.split()[1]
                    
                except (Exception, e):
                    print(e)
                    
                i +=3
                annotation = True
                
                while annotation:
                    
                    line = lines[i]

                    if not line:
                        annotation = False
                        i +=1
                        break
                    
                    domain_annotation = line.strip('\n').split()
                    inclusion = domain_annotation[1]
                    
                    domain = Domain()
                    
                    domain.query_hmm = query_hmm
                    
                    domain.id = domain_annotation[0]
                    domain.inclusion = domain_annotation[1]
                    domain.bit_score = float(domain_annotation[2])
                    domain.bias = float(domain_annotation[3])
                    
                    ## conditional e-value
                    #domain.e_value = float(domain_annotation[4])
                    
                    # indipendent e-value
                    domain.e_value = float(domain_annotation[5])
                    
                    domain.alifrom = float(domain_annotation[9])
                    domain.alito = float(domain_annotation[10])
                    
                    found_domains.append(domain)
                    
                    i +=1
                
                line = lines[i]
                if 'Alignments for each domain' in line:
                    i +=1
                    alignments = True
                
                while alignments:

                    for j in xrange(len(found_domains)):
                        
                        line = lines[i]
                        domain_id = 'domain %i' % (j+1) 
                        
                        if domain_id in line:
                            
                            domain = next((d for d in found_domains if d.id==str(j+1)), None)
                            
                            while len(lines[i+1].split()) < 3:
                                i += 1
                            query_line = lines[i+1]
 
                            domain.matchfrom = int(query_line.split()[1])
                            domain.model_match = query_line.split()[2]
                            domain.matchto = int(query_line.split()[3])
                            
                            target_line = lines[i+3]
                            domain.alifrom = int(target_line.split()[1])
                            domain.sequence_match = target_line.split()[2]
                            domain.alito = int(target_line.split()[3])
                            
                            i += 6
                            
                    alignments = False
                
                match.found_domains.extend(found_domains)
                found_domains = []
        i += 1
        
    return match

def parse_hmmsearch_output(filename):

    matches = []
   
    with open(filename,"r") as fl:
        lines = fl.readlines()
    lines = map(lambda x: x.strip('\n'), lines)

    query_hmm = filename.split('/')[-1].split('.')[0]
     
    run = False
    annotation = False
    alignments = False
    
    i = 1
    while i < len(lines):
        
        line = lines[i]
        
        if 'Domain annotation for each' in line:
            i += 1
            run = True
        
        while run:
            
            line = lines[i]
            
            if not line:
                i += 1
                continue
            
            if "No targets detected" in line:
                return []
                
             # Now we stop reading the file
            if line.startswith("//") or ("Internal pipeline statistics summary" in line):
                return matches
            
            if line.startswith(">>"):
                
                if "No individual domain" in lines[i+1]:
                    i += 3
                    continue
                
                match = Matched_Sequence()
                
                try:
                    match.name = line.split()[1]
                except(Exception, e):
                    print(e)
                try:
                    match.description = ' '.join(line.split()[2:])
                except(Exception, e):
                    print(e)
                    
                i +=3
                annotation = True
                while annotation:
                    
                    line = lines[i]

                    if not line:
                        annotation = False
                        i +=1
                        break
                    
                    domain_annotation = line.strip('\n').split()
                    inclusion = domain_annotation[1]
                    
                    domain = Domain()
                    
                    domain.query_hmm = query_hmm
                    
                    domain.id = domain_annotation[0]
                    domain.inclusion = domain_annotation[1]
                    domain.bit_score = domain_annotation[2]
                    domain.e_value = float(domain_annotation[5])
                    domain.alifrom = float(domain_annotation[9])
                    domain.alito = float(domain_annotation[10])
                    
                    match.found_domains.append(domain)
                    
                    i +=1
                
                line = lines[i]
                if 'Alignments for each domain' in line:
                    i +=1
                    alignments = True
                
                while alignments:

                    for j in xrange(len(match.found_domains)):
                        
                        line = lines[i]
                        domain_id = 'domain %i' % (j+1) 
                        
                        if domain_id in line:
                            
                            domain = next((d for d in match.found_domains if d.id==str(j+1)), None)
                            
                            while len(lines[i+1].split()) < 3:
                                i += 1
                            query_line = lines[i+1]
 
                            domain.matchfrom = int(query_line.split()[1])
                            domain.model_match = query_line.split()[2]
                            domain.matchto = int(query_line.split()[3])
                            
                            target_line = lines[i+3]
                            domain.alifrom = int(target_line.split()[1])
                            domain.sequence_match = target_line.split()[2]
                            domain.alito = int(target_line.split()[3])
                            
                            i += 6
                            
                    alignments = False
            
                matches.append(match) 
        i += 1
        
    return matches

def convert(convert_funcs, seq):
    return [
        item if func is None else func(item)
        for func, item in zip(convert_funcs, seq)
        ]

def blocks(file_object, chunk_size=1024):
    """Lazy function (generator) to read a file piece by piece.
    Default chunk size: 1k."""
    while True:
        data = file_object.read(chunk_size)
        if not data:
            break
        yield data
        
