#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from operator import attrgetter
from itertools import chain
#import commands

from lib.data_structures import Fragment_profile, MFS_profile, Pfam_profile
from lib.utilities import parse_hmmscan_output, get_number_of_gaps, merge_ranges
from lib.utilities import read_fasta_as_dict, write_fasta_dict_to_file

import subprocess

class Search_Manager(object):

    def __init__(self):
        
        this_file_path = os.path.realpath(__file__)
        self.root = "/".join(this_file_path.split("/")[0:-2])
        self.tool = self.get_path_hmm()
        
        
    def set_metal(self, metal):
        self.metal = metal
        self.db_path = os.path.join(self.root, 'data/hmm_profiles_db/%s'%metal)
        
    def get_path_hmm(self):
       
        this_file_path = os.path.realpath(__file__)
        local_path = "/".join(this_file_path.split("/")[0:-2])
 
        file_path = open(os.path.join(local_path, "data/local_path_hmm.txt"), "r")
        line_path = file_path.readlines()
        file_path.close()
        
        return line_path[0].strip()

    def set_db_path(self, db_path):
        self.db_path = db_path

    def set_blacklist(self):
        
        self.blacklist_db = os.path.join(self.db_path, 'pfam_blacklist_db.hmm')

        blacklist_profiles_path = os.path.join(self.db_path, 'pfam_blacklist.txt')
        self.blacklist_profiles = []
        with open(blacklist_profiles_path, 'r') as fl:
            self.blacklist_profiles = map(lambda x: x.strip('\n').split('|'), fl.readlines())
            
        if self.blacklist_profiles != []:
            return True
        else:
            return False

    def load_input_sequences(self, file_path):        
        self.input_sequences = read_fasta_as_dict(file_path)

    def get_detected_sequences(self, userdir):
        
        # process hmmscan tool output
        
        #temp_sequence_file = os.path.join(userdir, '%s_seqtemp' % self.method)
        #temp_output_file = os.path.join(userdir, '%s_outtemp' % self.method)
        
        matches = []

        i = 0
        for sequence_name, sequence in self.input_sequences.items():

            temp_sequence_file = os.path.join(userdir, '%s_%s_seqtemp' % (self.method, str(i)))
            temp_output_file = os.path.join(userdir, '%s_%s_outtemp' % (self.method, str(i)))

            seqdict = {sequence_name: sequence}
            write_fasta_dict_to_file(temp_sequence_file, seqdict)

            run_opt = [self.tool, '--notextw', '-o', temp_output_file, self.db, temp_sequence_file]
            subprocess.Popen(run_opt, stdout=subprocess.PIPE).communicate()[0]

            match = parse_hmmscan_output(temp_output_file)
            if match:
                match.sequence = sequence
                matches.append(match)
            
            os.remove(temp_sequence_file)
            os.remove(temp_output_file)
            i+=1
            
        return matches

    def get_profile_by_name(self, hmm_name):
        
        for profile in self.profiles:
            if profile.name == hmm_name:
                return profile

    def get_blacklist_profile(self, name):

        profile = next((p for p in self.blacklist_profiles if p[0].strip()==name.strip()), None)
        if profile:
            return profile[0], profile[1]
        else:
            return None, None

    def filter_blacklisted_domains(self, userdir, detected_sequences):
        
        p = 0
        s = 0
        slen = len(detected_sequences)
        
        while p < slen:
            
            sequence = detected_sequences[s]
            
            temp_sequence_file = os.path.join(userdir, '%s_seqtemp' % self.method)
            temp_output_file = os.path.join(userdir, '%s_outtemp' % self.method)
            
            seqdict = {sequence.name: sequence.sequence}
            write_fasta_dict_to_file(temp_sequence_file, seqdict)

            run_opt = [self.tool, '--notextw', '-o', temp_output_file, self.blacklist_db, temp_sequence_file]
            subprocess.Popen(run_opt, stdout=subprocess.PIPE).communicate()[0]

            match = parse_hmmscan_output(temp_output_file)
            
            if match: 
                
                i = 0
                j = 0
                seq_domn = len(sequence.found_domains)
                
                while j < seq_domn:
                    
                    found_domain = sequence.found_domains[i]
                    positions_of_predicted_ligands = sequence.predicted_ligands[i]
                    
                    k = 0
                    l = 0
                    fdomn = len(match.found_domains)
                    
                    deleted = False
                    
                    while k < fdomn:

                        blacklist_domain = match.found_domains[l]
                        if blacklist_domain.e_value >= self.max_evalue:
                            del match.found_domains[l]
                            k += 1
                            continue
                        
                        # not an empty list for predictions using structural information
                        if positions_of_predicted_ligands != []:
                            
                            profile_name, profile_pattern = self.get_blacklist_profile(blacklist_domain.query_hmm)
                            pligands_aas = map(lambda x: x[0], profile_pattern.split('_'))
                            pligands_positions = map(lambda x: int(x[1:]), profile_pattern.split('_'))
                            
                            # get positions of ligands in domain alignment of sequence with a blacklisted profile
                            positions_of_ligands = blacklist_domain.get_positions_of_ligands(pligands_aas, pligands_positions)
                            
                            for position in positions_of_predicted_ligands:   
                                if set(positions_of_ligands).issubset(set(position)):
                                    del sequence.found_domains[i]
                                    j += 1
                                    deleted = True
                                    break
                            
                            if deleted:
                                break
                            
                            k += 1
                            l += 1
                                
                        # no structural information
                        else:
                            if max(blacklist_domain.alifrom,found_domain.alifrom) <= min(blacklist_domain.alito,found_domain.alito):
                                del sequence.found_domains[i]
                                j += 1
                                deleted = True
                                break
                            
                            k += 1
                            l += 1
                            
                    if not deleted:
                        i += 1
                        j += 1
                        
            if sequence.found_domains == []:
                del detected_sequences[s]
                if not detected_sequences:
                    return
            else:
                 s += 1
            p += 1
            
            os.remove(temp_sequence_file)
            os.remove(temp_output_file)

    def get_summary_annotation(self, sequence):
        
        best_domain_index = sequence.found_domains.index(min(sequence.found_domains,key=attrgetter('e_value')))
        bit_score = sequence.found_domains[best_domain_index].bit_score
        bias = sequence.found_domains[best_domain_index].bias
        evalue = sequence.found_domains[best_domain_index].e_value
        
        predicted_ligands = sorted(list(set(list(chain.from_iterable(list(chain.from_iterable(sequence.predicted_ligands)))))))
        intervals, ligands = self.get_maximum_coverage_for_domains(sequence.found_domains, predicted_ligands)
        annotated_sequence = sequence.annotate_domains_coverage(intervals, ligands)
       
        #annotation = "{hdr}@{seq}@{bit}@{bs}@{evl:.2E}\n".format(hdr = sequence.header,\
        #                                                        seq = annotated_sequence,\
        #                                                        bit = bit_score,\
        #                                                        bs = bias,\
        #                                                        evl = evalue)
        
        annotation = "%s,%s,%s,%s,%s\n"%(sequence.header, annotated_sequence, bit_score, bias, evalue)
        
        return annotation

    def get_maximum_coverage_for_domains(self, found_domains, predicted_ligands):

        intervals = []
        ligands = []

        start_ends = []            
        for domain in found_domains:
            st_en = [domain.alifrom-1,domain.alito-1]
            start_ends.append(filter(lambda a: a != None, st_en))
        
        for group in merge_ranges(start_ends):
            intervals.append(group)
            ligs = []
            for l in predicted_ligands:
                if l >= group[0] and l <= group[1]:
                    ligs.append(l)
            ligands.append(sorted(ligs))
        
        il = sorted(zip(intervals, ligands), key=lambda x: x[0][0])
        intervals, ligands = zip(*il)
        
        return intervals, ligands

    def get_sequences_annotations(self, matched_sequences, outfile_summary, outfile_detailed):
        
        summary_annotation = []
        detailed_annotation = []
        report_data = []
        
        header_summ_csv = "Seq header,Sequence,Bit score,Bias,Evalue\n"
        summary_annotation.append(header_summ_csv)
        
        header_dett_csv = self.get_header_csv()
        detailed_annotation.append(header_dett_csv)
        
        for msequence in matched_sequences:

            header, sequence = next(((h,s) for h,s in self.input_sequences.items() if msequence.name in h), None)

            header_clean = header.replace(",", ";")
            if "@" in header_clean:
               header_list = header_clean.split("@")[1:]
               header_clean = "".join(header_list)
            msequence.header = header_clean
            msequence.sequence = sequence

            annotated_sequence_summary = self.get_summary_annotation(msequence)
            summary_annotation.append(annotated_sequence_summary)
            
            annotated_sequence_detailes = self.get_detailed_annotation(msequence)
            detailed_annotation.extend(annotated_sequence_detailes)
        
        with open(outfile_summary, 'w') as fl:
            fl.writelines(summary_annotation)
            
        with open(outfile_detailed, 'w') as fl:
            fl.writelines(detailed_annotation)
            
class Fragments_Manager(Search_Manager):

    def __init__(self, max_evalue = float(1e-2)):
        
        Search_Manager.__init__(self)
        self.max_evalue = max_evalue
        
    def set_database(self, method):
        
        self.method = method
        self.db = os.path.join(self.db_path, 'fragments_db.hmm')
        self.profiles_path = os.path.join(self.db_path, 'fragments_patterns.txt')
        
    def create_library_of_profiles(self):
        
        # returns a list of MFS_profile objects (see data_structures.py library to refresh memory on MFS_profile class)
        
        self.profiles = []
        with open(self.profiles_path, 'r') as fl:
            profiles_info_notf = map(lambda x: x.strip('\n').split('|'), fl.readlines())
           

        if profiles_info_notf==[]:
           return False 
        
        ##VALE FOR FUNCTION TMP: Poi risostituisci 'profiles_info_notf' con 'profiles_info' e cancella il ciclo qui sotto!
        profiles_info = []
        for pinf in profiles_info_notf:
            profiles_info.append(pinf[:3])
        ##END VALE TMP

        
        for profile_name, pattern, ligands_positions in profiles_info:
            
            fields = profile_name.split('_')
            
            pdb_code = fields[0]
            chain_id = fields[1]
            site_id = '_'.join(fields[2:4])
            
            # Fragment profile is an actual searched HMM profile
            profile = Fragment_profile()
            profile.name = profile_name
            profile.pdb_chain = '%s_%s' % (pdb_code, chain_id)
            profile.site_id = site_id
            profile.id = int(fields[4].strip('F'))
            
            if pattern != '':
                
                profile.has_ligands = True
                ligands = pattern.split('_')
                profile.ligands_aas = map(lambda ligand: ligand[0], ligands)
                profile.ligands_resids = map(lambda ligand: int(ligand[1:]), ligands)
            
                profile.ligands_positions_in_model.extend(map(lambda x: int(x), ligands_positions.split('_')))
            
            self.profiles.append(profile)

        return True




    def filter_matched_domains(self, detected_sequences):

        k = 0
        s = 0
        seqn = len(detected_sequences)
        while s < seqn:
            
            sequence = detected_sequences[k]
        
            fragment_match = False

            i = 0
            j = 0
            domn = len(sequence.found_domains)
            
            while j < domn:
                
                domain = sequence.found_domains[i]
                if domain.e_value > self.max_evalue:
                    del sequence.found_domains[i]
                    j += 1
                    continue
                
                fragment_profile = self.get_profile_by_name(domain.query_hmm)

                if fragment_profile.has_ligands:
                    
                    ligands_aas = fragment_profile.get_ligands_aas()
                    ligands_positions = fragment_profile.get_ligands_positions_in_model()

                    positions_of_ligands = domain.get_positions_of_ligands(ligands_aas, ligands_positions)
                    
                    if positions_of_ligands:
                        sequence.predicted_ligands.append([positions_of_ligands])
                        fragment_match = True
                        i += 1
                        
                    else:
                        del sequence.found_domains[i]
                        if not sequence.found_domains:
                            break
                else:
                    sequence.predicted_ligands.append([])
                    i += 1
                    
                j += 1


            if self.metal=='molybdenum':
               fragment_match=True
            
            if not sequence.found_domains or not fragment_match:
                del detected_sequences[k]
                if not detected_sequences:
                    return []
            else:
                k += 1
            s += 1
            
    def get_header_csv(self):
        
        head = "Seq header,Sequence,Pdb code,Site id,Profile id,Fragment ligands,Bit score,Bias,Evalue,Predicted ligands,Fragment range\n"
        return head

    def get_detailed_annotation(self, sequence):
        
        annotations = []
   
        for i in xrange(len(sequence.found_domains)):
            
            domain = sequence.found_domains[i]
            domain_start = domain.alifrom - 1
            domain_end = domain.alito - 1
            fragment_range = '%s-%s' % (domain.alifrom, domain.alito)
            
            profile = self.get_profile_by_name(domain.query_hmm)
            
            ligands_positions = []
            predicted_ligands = ''
            fragment_ligands = ''
            
            if profile.has_ligands:
                
                fragment_ligands_list = map(lambda x: '%s%i'%(x[0], x[1]), zip(profile.ligands_aas, profile.ligands_resids))
                fragment_ligands = '-'.join(fragment_ligands_list)

                ligands_positions = sequence.predicted_ligands[i][0]
                predicted_ligands_aas = sequence.get_predicted_ligands(ligands_positions)
                predicted_ligands_list = map(lambda x: '%s%i'%(x[0], x[1]), zip(predicted_ligands_aas, map(lambda y: y+1,ligands_positions)))
                predicted_ligands = '-'.join(predicted_ligands_list)
            
            annotated_sequence = sequence.annotate_domain(domain_start, domain_end, ligands_positions)

            #annotation = "{hdr}@{seq}@{pdb}@{st}@{pid}@{lig}@{bit}@{bs}@{evl:.2E}@{plig}@{rng}\n".format(hdr = sequence.header,\
            #                                                                                        seq = annotated_sequence,\
            #                                                                                        pdb = profile.pdb_chain,\
            #                                                                                        st = profile.site_id,\
            #                                                                                        pid = profile.id,\
            #                                                                                        lig = fragment_ligands,\
            #                                                                                        bit = domain.bit_score,\
            #                                                                                        bs = domain.bias,\
            #                                                                                        evl = domain.e_value,\
            #                                                                                        plig = predicted_ligands,\
            #                                                                                        rng = fragment_range)
            
            annotation = "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%(sequence.header,annotated_sequence, profile.pdb_chain, profile.site_id, profile.id, fragment_ligands, domain.bit_score, domain.bias, domain.e_value, predicted_ligands, fragment_range)
            
            
            annotations.append(annotation)

        return annotations
    
class MFS_Manager(Fragments_Manager):

    def __init__(self, max_evalue=float(1e-2)):
        
        Fragments_Manager.__init__(self)
        self.max_evalue = max_evalue
        
    def create_library_of_profiles(self):
        
        # returns a list of MFS_profile objects (see data_structures.py library to refresh memory on MFS_profile class)
        
        self.profiles = []
        with open(self.profiles_path, 'r') as fl:
            profiles_info_notf = map(lambda x: x.strip('\n').split('|'), fl.readlines())
        
        if profiles_info_notf==[]:
           return False
    
        ##VALE FOR FUNCTION TMP: Poi risostituisci 'profiles_info_notf' con 'profiles_info' e cancella il ciclo qui sotto!
        profiles_info = []
        for pinf in profiles_info_notf:
            profiles_info.append(pinf[:3])
        ##END VALE TMP
        
        for profile_name, pattern, ligands_positions in profiles_info:
            
            fields = profile_name.split('_')
            
            pdb_code = fields[0]
            chain_id = fields[1]
            site_id = '_'.join(fields[2:4])
            
            name = '%s_%s_%s' % (pdb_code, chain_id, site_id)
            
            mfs_profile = next((p for p in self.profiles if p.name==name), None)
            
            if mfs_profile == None:
                
                # MFS profile is composed of profiles of fragments and has some additional metadata
                mfs_profile = MFS_profile()
                
                mfs_profile.name = '%s_%s_%s' % (pdb_code, chain_id, site_id)
                mfs_profile.pdb_chain = '%s_%s' % (pdb_code, chain_id)
                mfs_profile.site_id = site_id
                
                self.profiles.append(mfs_profile)
            
            # Fragment profile is an actual searched HMM profile
            fragment_profile = Fragment_profile()
            fragment_profile.name = profile_name
            fragment_profile.pdb_chain = '%s_%s' % (pdb_code, chain_id)
            fragment_profile.site_id = site_id
            fragment_profile.id = int(fields[4].strip('F'))

            if pattern != '':
                
                fragment_profile.has_ligands = True
                ligands = pattern.split('_')
                fragment_profile.ligands_aas = map(lambda ligand: ligand[0], ligands)
                fragment_profile.ligands_resids = map(lambda ligand: int(ligand[1:]), ligands)
            
                fragment_profile.ligands_positions_in_model.extend(map(lambda x: int(x), ligands_positions.split('_')))

            mfs_profile.fragments.append(fragment_profile)
            
            mfs_profile.ligands_aas.append(fragment_profile.ligands_aas)
            mfs_profile.ligands_resids.append(fragment_profile.ligands_resids)
            mfs_profile.ligands_positions_in_model.append(fragment_profile.ligands_aas)
    
        return True



    def group_domains_matched_by_profile(self, domains):
        
        domains_matched_by_mfs_profile = {}
        
        for i in xrange(len(domains)):
            
            domain = domains[i]

            query_hmm = domain.query_hmm
            profile_name = '_'.join(query_hmm.split('_')[:4])
            
            pname = next((pn for pn in domains_matched_by_mfs_profile.keys() if pn == profile_name), None)
            if pname:
                domains_matched_by_mfs_profile[pname].append(domain)
            else:
                domains_matched_by_mfs_profile[profile_name]=[domain]
        
        return domains_matched_by_mfs_profile

    def filter_matched_domains(self, detected_sequences):
        
        k = 0
        s = 0
        seqn = len(detected_sequences)
        while s < seqn:
            
            sequence = detected_sequences[k]

            matches = self.group_domains_matched_by_profile(sequence.found_domains)
            
            found_domains = []
            predicted_ligands = []
            for mfs_profile_name, matched_domains in matches.items():

                mfs_profile = self.get_profile_by_name(mfs_profile_name)
                
                site_match = True
                predicted_ligands_fragments = []
                matched_domains_fragments = []
                
                for i in xrange(len(mfs_profile.fragments)):
                    
                    fid = i+1
                    fragment_profile = mfs_profile.get_fragment_by_id(fid)
                    
                    if not fragment_profile:
                        continue
                    
                    domains = [d for d in matched_domains if d.query_hmm==fragment_profile.name]
                    if fragment_profile.has_ligands and not domains:
                        site_match = False
                        matched_domains_fragments = []
                        predicted_ligands_fragments = []
                        break
                    
                    elif fragment_profile.has_ligands:
                        
                        ligands_aas = fragment_profile.get_ligands_aas()
                        ligands_positions = fragment_profile.get_ligands_positions_in_model()
                        
                        ligands_match = False
                        for domain in domains:
                            if domain.e_value > self.max_evalue:
                                continue
                            positions_of_ligands = domain.get_positions_of_ligands(ligands_aas, ligands_positions)
                            if positions_of_ligands:
                                ligands_match = True
                                matched_domains_fragments.append(domain)
                                predicted_ligands_fragments.append([positions_of_ligands])
                        if not ligands_match:
                            site_match = False
                            matched_domains_fragments = []
                            predicted_ligands_fragments = []
                            
                    elif ((fragment_profile.has_ligands is False) and (domains != [])):
                        for domain in domains:
                            matched_domains_fragments.append(domain)
                            predicted_ligands_fragments.append([])
                
                if site_match:
                    found_domains.extend(matched_domains_fragments)
                    predicted_ligands.extend(predicted_ligands_fragments)
            
            if found_domains:
                sequence.found_domains = found_domains
                sequence.predicted_ligands = predicted_ligands
                k += 1
            else:
                del detected_sequences[k]
                if not detected_sequences:
                    return []
            s += 1
            
            
    def get_header_csv(self):
        
        head = "Seq header,Sequence,Pdb code,Site id,Ligands,Bit score,Bias,Evalue,Predicted ligands,Ranges\n"
        return head
    
    def get_detailed_annotation(self, sequence):
        
        annotations = []
        
        matches = self.group_domains_matched_by_profile(sequence.found_domains)
        
        for profile_name, domains in matches.items():
            
            profile = self.get_profile_by_name(profile_name)

            ligands_aas, ligands_resids = profile.get_ligands()
            if not ligands_aas:
               if self.metal!='molybdenum':
                  continue
               else:
                  mo_nolig = True
            else:
               mo_nolig = False

            if mo_nolig:
               ligands_positions = []
            else:
               ligands_list = map(lambda x: '%s%i'%(x[0], x[1]), zip(ligands_aas, ligands_resids))
               ligands = '-'.join(ligands_list)
            
               ligands_positions = []
               for domain in domains:
                   dom_ind = sequence.found_domains.index(domain)
                   predicted_ligands_domain = sequence.predicted_ligands[dom_ind]
                   ligands_positions.extend(predicted_ligands_domain)
               ligands_positions = [x for l in ligands_positions for x in l]
               ligands_positions = sorted(list(set(ligands_positions)))
            
            domains_intervals, ligands_intervals = self.get_maximum_coverage_for_domains(domains, ligands_positions)
            annotated_sequence = sequence.annotate_domains_coverage(domains_intervals, ligands_intervals)
            
            # best e-value
            best_domain_index = domains.index(min(domains,key=attrgetter('e_value')))
            bit_score = domains[best_domain_index].bit_score
            bias = domains[best_domain_index].bias
            evalue = domains[best_domain_index].e_value
           
            if mo_nolig:
               ligands = ''
               predicted_ligands = ''
            else:
               # predicted ligands
               predicted_ligands_aas = sequence.get_predicted_ligands(ligands_positions)
               predicted_ligands_list = map(lambda x: '%s%i'%(x[0], x[1]), zip(predicted_ligands_aas, map(lambda y: y+1,ligands_positions)))
               predicted_ligands = '-'.join(predicted_ligands_list)
            
            # domains ranges
            domains_ranges = map(lambda x: '-'.join(map(lambda y: str(y+1), x)), domains_intervals)
            ranges = ','.join(domains_ranges)

            #annotation = "{hdr}@{seq}@{pdb}@{st}@{lig}@{bit}@{bs}@{evl:.2E}@{plig}@{rng}\n".format(hdr = sequence.header,\
            #                                                                                        seq = annotated_sequence,\
            #                                                                                        pdb = profile.pdb_chain,\
            #                                                                                        st = profile.site_id,\
            #                                                                                        lig = ligands,\
            #                                                                                        bit = bit_score,\
            #                                                                                        bs = bias,\
            #                                                                                        evl = evalue,\
            #                                                                                        plig = predicted_ligands,\
            #                                                                                        rng = ranges)
            
            annotation = "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%(sequence.header, annotated_sequence, profile.pdb_chain, profile.site_id, ligands, bit_score, bias, evalue, predicted_ligands, ranges)
            
            annotations.append(annotation)

        return annotations
    
class Pfam_Manager(Search_Manager):

    def __init__(self, max_evalue=float(1e-5)):
        
        Search_Manager.__init__(self)
        self.max_evalue = max_evalue
        
    def set_database(self, method='pfam'):
        
        self.method = method

        self.db = os.path.join(self.db_path, '%s_db.hmm' % method)
        self.profiles_path = os.path.join(self.db_path, '%s_patterns.txt' % method)
        self.blacklist_profiles_path = os.path.join(self.db_path, 'pfam_blacklist.txt')
        
    def create_library_of_profiles(self):
        
        self.profiles = []
        with open(self.profiles_path, 'r') as fl:
            profiles_info_notf = map(lambda x: x.strip('\n'), fl.readlines())
           
        if profiles_info_notf==[]:
           return False
 
        ##VALE FOR FUNCTION TMP: Poi risostituisci 'profiles_info_notf' con 'profiles_info' e cancella il ciclo qui sotto!
        profiles_info = []
        for pinf in profiles_info_notf:
            pinf_list = pinf.split("|")
            profiles_info.append("|".join(pinf_list[:2]))
        ##END VALE TMP
        
        profiles_names = list(set(map(lambda x: x.split('|')[0], profiles_info)))

        for profile_name in profiles_names:
            
            patterns = []
            for line in profiles_info:
                if profile_name in line:
                    name, pattern = line.strip('\n').split('|')
                    if pattern != '':
                        patterns.append(pattern)
            
            profile = Pfam_profile()
            profile.name = profile_name
           
            if patterns != []:
                
                profile.has_structure = True
                
                lig_aas = []
                lig_positions_in_model = []
                for pattern in patterns:
        

                    plig_aas = map(lambda x: x[1], pattern.split('_'))
                    plig_positions = map(lambda x: int(x[2:]), pattern.split('_'))
                    
                    lig_aas.append(plig_aas)
                    lig_positions_in_model.append(plig_positions)
                
                il = sorted(zip(lig_aas, lig_positions_in_model), key=lambda x: len(x[0]), reverse=True)
                profile.ligands_aas, profile.ligands_positions_in_model = zip(*il)
        
            self.profiles.append(profile)

        return True


    def get_profile_by_name(self, hmm_name):
        
        for profile in self.profiles:
            if profile.name == hmm_name:
                return profile

    def filter_matched_domains(self, detected_sequences):
        
        k = 0
        s = 0
        seqn = len(detected_sequences)
        while s < seqn:
            
            sequence = detected_sequences[k]
            
            i = 0
            j = 0
            domn = len(sequence.found_domains)
            while j < domn:
                
                domain = sequence.found_domains[i]

                if domain.e_value > self.max_evalue:
                    del sequence.found_domains[i]
                    j += 1
                    continue

                profile = self.get_profile_by_name(domain.query_hmm)
                if not profile:
                   del sequence.found_domains[i]
                   j += 1
                   continue               
 
                if profile.has_structure:

                    patterns_ligands_aas = profile.get_ligands_aas()
                    patterns_ligands_positions = profile.get_ligands_positions_in_model()
                    ligands_are_conserved = False
                    
                    posi_save = []

                    for p in xrange(len(patterns_ligands_aas)):
                        
                        pligands_aas = patterns_ligands_aas[p]
                        pligands_positions = patterns_ligands_positions[p]
                        
                        positions_of_ligands = domain.get_positions_of_ligands(pligands_aas, pligands_positions)
                        
                        if positions_of_ligands: 
                            if len(sequence.predicted_ligands) == i+1:
                                pattern_flag = True
                                for pattern in sequence.predicted_ligands[i]:
                                    if set(positions_of_ligands).issubset(set(pattern)):
                                        pattern_flag = False
                                        break
                                if pattern_flag:
                                    sequence.predicted_ligands[i].append(positions_of_ligands)
                                    posi_save.append(map(str, positions_of_ligands))
                            else:
                                sequence.predicted_ligands.append([positions_of_ligands])
                                posi_save.append(map(str, positions_of_ligands))
                            
                            ligands_are_conserved = True
 
                    if not ligands_are_conserved:
                        del sequence.found_domains[i]
                        if not sequence.found_domains:
                            break
                    else:
                        for psv in posi_save:
                            sequence.single_pfamligands.append("_".join(psv))
                        i += 1
                else:
                    sequence.predicted_ligands.append([])
                    i += 1
                j += 1
            
            if not sequence.found_domains:
                del detected_sequences[k]
                if not detected_sequences:
                    return []
            else:
                k += 1
            s += 1
            
            
    def get_header_csv(self):
        
        head = "Seq header,Sequence,Domain query,Bit score,Bias,Evalue,Predicted ligands,Domain range\n"
        return head


    def get_detailed_annotation(self, sequence):
        
        annotations = []
        
        for i in xrange(len(sequence.found_domains)):
            
            domain = sequence.found_domains[i]
            domain_start = domain.alifrom - 1
            domain_end = domain.alito - 1
            domain_range = '%s-%s' % (domain.alifrom, domain.alito)            

            if sequence.predicted_ligands[i]:
                
                predicted_ligands_positions = sequence.predicted_ligands[i]
                for j in xrange(len(predicted_ligands_positions)):
                    
                    ligands_positions = predicted_ligands_positions[j]
                    predicted_ligands_aas = sequence.get_predicted_ligands(ligands_positions)
                    predicted_ligands_list = map(lambda x: '%s%i'%(x[0], x[1]), zip(predicted_ligands_aas, map(lambda y: y+1,ligands_positions)))
                    predicted_ligands = '-'.join(predicted_ligands_list)
                    
                    annotated_sequence = sequence.annotate_domain(domain_start, domain_end, ligands_positions)

                    posi_setted = list(set(sequence.single_pfamligands))
                    single_ligand_pattern = ";".join(posi_setted)

                    #annotation = "{hdr}@{seq}@{hmm}@{bit}@{bs}@{evl:.2E}@{plig}@{rng}\n".format(hdr = sequence.header,\
                    #                                                                            seq = annotated_sequence,\
                    #                                                                            hmm = domain.query_hmm,\
                    #                                                                            bit = domain.bit_score,\
                    #                                                                            bs = domain.bias,\
                    #                                                                            evl = domain.e_value,\
                    #                                                                            plig = predicted_ligands,\
                    #                                                                            rng = domain_range)
                    
                    annotation = "%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%(sequence.header, annotated_sequence, domain.query_hmm, domain.bit_score, domain.bias, domain.e_value, predicted_ligands, domain_range, single_ligand_pattern)
                    
                    annotations.append(annotation)
            else:
                annotated_sequence = sequence.annotate_domain(domain_start, domain_end, [])
                
                #annotation = "{hdr}@{seq}@{hmm}@{bit}@{bs}@{evl:.2E}@{rng}\n".format(hdr = sequence.header,\
                #                                                                    seq = annotated_sequence,\
                #                                                                    hmm = domain.query_hmm,\
                #                                                                    bit = domain.bit_score,\
                #                                                                    bs = domain.bias,\
                #                                                                    evl = domain.e_value,\
                #                                                                    rng = domain_range)
                
                annotation = "%s,%s,%s,%s,%s,%s,-,%s\n"%(sequence.header, annotated_sequence, domain.query_hmm, domain.bit_score, domain.bias, domain.e_value, domain_range)
                
                annotations.append(annotation)

        return annotations
