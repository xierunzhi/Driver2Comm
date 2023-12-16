import random
import time

import pandas as pd

from graph import *
from external import External
import os
from association_test import AssociationTest
class shuffle_evaluation(object):
    def __init__(self,c2c_networks,patient_info,driver):
        self.c2c_networks = c2c_networks
        self.patient_info = patient_info
        self.driver = driver
    def shuffle_c2c_networks(self):
        """
        shuffle the cell-cell communication network in order to compute the empirical p value
        :return:
        """
        c2c_networks = self.c2c_networks
        shuffled_c2c_networks = dict()
        for patient_cnt, network in c2c_networks.items():
            shuffled_network = Graph(network.gid, network.is_undirected, network.eid_auto_increment)
            # vertices_num = len(network.vertices)
            set_of_vlb = set()
            for vid, vertex in network.vertices.items():
                set_of_vlb.add(vertex.vlb)
            list_of_vlb = list(set_of_vlb)
            random.shuffle(list_of_vlb)  # key step of this function
            for vid, vertex in network.vertices.items():
                # randvid = random.randint(0,len(list_of_vlb)-1)     #will choose the same vertex label for a vertex
                randvid = int(vid)
                shuffled_network.add_vertex(vid, list_of_vlb[randvid])
            # there are two types of edges in network elb
            for elb, set_of_edge in network.set_of_elb.items():
                set_of_edge = list(set_of_edge)
                for frm, to in set_of_edge:
                    shuffled_network.add_edge(AUTO_EDGE_ID, frm, to, elb)
            shuffled_c2c_networks[patient_cnt] = shuffled_network
        return shuffled_c2c_networks

    def filter_celltype(self,gene:str):
        if 'Macrophages' in gene:
            return gene.replace('__Macrophages','')
        elif 'Cd8+Tcells' in gene:
            return gene.replace('__Cd8+Tcells','')
        elif 'Tumor' in gene:
            return gene.replace('__Tumor','')
        else:
            raise Exception('Something wrong in cell type')
    def generate_candidate_targets(self,shuffled_fp):
        ret = list()
        targets_set = set()
        for i in range(len(shuffled_fp)):
            candidate_graph = shuffled_fp[i].to_graph()
            for vid,vertex in candidate_graph.vertices.items():
                external_target = self.filter_celltype(vertex.vlb)
                if external_target.upper() not in targets_set:
                    targets_set.add(external_target.upper())
                    ret.append((self.driver,external_target.upper()))
        return ret
    def run(self):
        """
        generate shuffled frequent pattern matrix to post analysis using R
        :return:
        """
        shuffled_candidate_targets = list()
        candidate_targets_set = set()
        start = time.time()
        for i in range(1000):
            shuffled_c2c_networks = self.shuffle_c2c_networks()
            external = External(shuffled_c2c_networks, min_support=3, patient_info=self.patient_info, output2screen=False)
            external.run()
            shuffled_fp = external.get_frequent_pattern()
            if shuffled_fp is None:
                continue
            candidate_targets = self.generate_candidate_targets(shuffled_fp)
            for target in candidate_targets:
                if target not in candidate_targets_set:
                    candidate_targets_set.add(target)
                    shuffled_candidate_targets.append(target)
            if(len(shuffled_candidate_targets)>5000):
                break
                #association_test_ret = self.association_test(shuffled_external_matrix)
            if (i + 1) % 5 == 0:
                print('It has been shuffle {} times'.format(i + 1))
                print('length of shuffle candidate_targets : {}'.format(len(shuffled_candidate_targets)))
                end = time.time()
                print('shuffle for {} times takes {}'.format(i,end - start))
        return shuffled_candidate_targets

