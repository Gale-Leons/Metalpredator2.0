#!/usr/bin/env python
# -*- coding: utf-8 -*-

class Matched_Sequence(object):
    
    def __init__(self):
        
        self.name = ''
        self.description = ''
        self.sequence = ''
        
        self.found_domains = []
        self.predicted_ligands = []
        self.single_pfamligands = []

    def get_predicted_ligands(self, ligands_positions):
        
        ligands = []
        for pos in ligands_positions:
            ligands.append(self.sequence[pos])
       
        return ligands

    def annotate_domain(self, domain_start, domain_end, predicted_ligands):
        
        annotation = []
        
        start_start = 0
        end_end = 0
        
        start = domain_start
        end = domain_end
        
        annotation.append(self.sequence[start_start:start])
        #annotation.append('<font class="match">')

        lstart = start
        for pos in predicted_ligands:
            annotation.append(self.sequence[lstart:pos])
            #annotation.append('<font class="ligand">')
            annotation.append(self.sequence[pos:pos+1])
            #annotation.append('</font>')
            annotation.append('*')   
            lstart=pos+1
                
        annotation.append(self.sequence[lstart:end+1])
        #annotation.append('</font>')

        end_end = end+1
        annotation.append(self.sequence[end_end:])
        annotated_sequence = ''.join(annotation)
                
        return annotated_sequence

    def annotate_domains_coverage(self, intervals, ligands):
        
        annotation = []

        start_start = 0
        end_end = 0
        for i in xrange(len(intervals)):

            start, end = intervals[i]
            ligs = ligands[i]
            
            annotation.append(self.sequence[start_start:start])
            #annotation.append('<font class="match">')

            lstart = start
            for pos in ligs:
                annotation.append(self.sequence[lstart:pos])
                #annotation.append('<font class="ligand">')
                annotation.append(self.sequence[pos:pos+1])
                #annotation.append('</font>')
                annotation.append('*')   
                lstart=pos+1

            annotation.append(self.sequence[lstart:end+1])
            #annotation.append('</font>')

            start_start = end+1
            end_end = end+1
            
        annotation.append(self.sequence[end_end:])
        annotated_sequence = ''.join(annotation)

        return annotated_sequence
    
class Domain(object):
    
    def __init__(self):
        
        self.id = None
        self.inclusion = None
        
        self.bit_score = None
        self.e_value = None

        self.matchfrom = None
        self.matchto = None
        self.model_match = ''
        
        self.alifrom = None
        self.alito = None
        self.sequence_match = ''
        
        self.query_hmm = ''

    def get_positions_of_ligands(self, ligands, ligands_positions):
        
        indiced_of_ligands_in_sequence = []
        
        for i in xrange(len(ligands_positions)):

            lig = ligands[i]
            pos = ligands_positions[i]
            
            if pos >= self.matchfrom and pos <= self.matchto:
                
                j = 0
                model_gaps = 0
                curr_pos = self.matchfrom
                
                while curr_pos <= pos:  
                    if self.model_match[j] != '.':
                        j += 1
                        curr_pos += 1
                    else:
                        j += 1
                        model_gaps += 1
                        
                match_ind = pos - self.matchfrom + model_gaps
                model_matched_res = self.sequence_match[match_ind]

                if model_matched_res == lig:
                    seq_gaps = self.sequence_match[:match_ind].count('-')
                    lig_ind_seq = self.alifrom+match_ind-seq_gaps-1
                    indiced_of_ligands_in_sequence.append(lig_ind_seq)
                else:
                    return []
                
        if len(ligands_positions) == len(indiced_of_ligands_in_sequence):  
            return indiced_of_ligands_in_sequence
        else:
            return []
