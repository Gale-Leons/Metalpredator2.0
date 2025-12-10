#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
from lib.utilities import get_number_of_gaps

class Profile_HMM(object):
    
    def __init__(self):
        
        self.name = ''
        self.has_structure = False
        
        self.ligands_aas = []
        self.ligands_resids = []
        self.ligands_positions_in_model = []

    def get_ligands_aas(self):
        return self.ligands_aas

    def get_ligands_positions_in_model(self):
        return self.ligands_positions_in_model
    
class Fragment_profile(Profile_HMM):

    def __init__(self):

        Profile_HMM.__init__(self)
        
        self.id = None
        self.pdb_chain = ''
        self.site_id = ''
        
        self.has_structure = True
        self.has_ligands = False

class MFS_profile(Profile_HMM):

    def __init__(self):

        Profile_HMM.__init__(self)

        self.pdb_chain = ''
        self.site_id = ''
        
        self.has_structure = True
        self.fragments = []
    
    def get_fragment_by_id(self, fid):
        return next((f for f in self.fragments if f.id==fid), None)
    
    def get_ligands(self):
        
        ligands_aas = []
        ligands_resids = []
        for i in range(len(self.ligands_aas)):
            fragment_ligands_aas = self.ligands_aas[i]
            ligands_aas.extend(fragment_ligands_aas)
            fragment_ligands_resids = self.ligands_resids[i]
            ligands_resids.extend(fragment_ligands_resids)

        if ligands_aas==[] and ligands_resids==[]:
           return False, False

        il = sorted(zip(ligands_aas, ligands_resids), key=lambda x: x[1])
        ligands_aas, ligands_resids = zip(*il)
        
        return ligands_aas, ligands_resids

class Pfam_profile(Profile_HMM):
    
    def __init__(self):
        
        Profile_HMM.__init__(self)
    
    def get_match_by_name(self, mid):
        return next((m for m in self.matched_sequences if m.name==mid), None)
    
    def get_patterns(self):
        
        #path = os.path.join(os.path.expanduser('~'), 'mfs_profiling/iron_sulfur/HMM_profiles_of_domains_patterns.txt')
        with open(path, 'r') as fl:
            lines = fl.readlines()
        patterns = []
        for line in lines:
            if self.name in line:
                name, pattern = line.strip('\n').split('|')
                patterns.append(pattern)

        return patterns
    
