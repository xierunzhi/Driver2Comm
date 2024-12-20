"""
 Driver2Comm is designed for identifying cell-type communication patterns associated with genetic variation in cancer.

"""

import pandas as pd
import numpy as np
import os
from .external import External
from .association_test import AssociationTest
import collections

class Driver2Comm(object):

    def __init__(self,c2c_network:dict,patient_driver_info:pd.DataFrame,minsup,association_test_threshold = 0.05,outputPATH = './result',
                 cancer_type = 'Brca'):
        """

        :param c2c_network: the formulated cell-cell communication network data
        :param patient_driver_info: the Driver information of patients, columns of this dataframe
               looks like [Patient_id,Driver1,Driver2], the Driver information should be encode as a one-hot vector
        :param minsup: Hyperparameter, the minimal suppport of gSpan algorithm
        :param outputPATH: the path where output files place
        :param cancer_type: type of cancer
        """

        self.minsup = minsup
        self.outputPATH = outputPATH
        self.c2c_network = c2c_network
        self.patient_info = dict()      # patient_info[Graph id] = patient label
        self.external_matrix = None
        self.internal_matrix = None
        self.patient_driver = dict()        # patient_driver[patient name] = patient driver
        self.frequent_subgraphs = None
        self.candidate_targets = None
        #self.drivers = dict()       #example: drivers[BRAF] = 0
        self.association_test_ret = None
        self.cancer_type = cancer_type
        self.association_test_threshold = association_test_threshold
        self.patient_driver_info = patient_driver_info
        if not os.path.exists(self.outputPATH):
            os.makedirs(self.outputPATH)
        self.get_patient_info(patient_driver_info)


    def get_patient_info(self,patient_metadata):
        """
        extract patient label and driver information from patient_metadata
        :param patient_metadata: pd.Dataframe
        :return:
        """
        assert 'Patient_id' in patient_metadata.columns, 'Please check again if the column name of patient id is correct!'
        # assert 'Driver' in patient_metadata.columns, 'Please check again if the column name of Driver is correct!'
        for i in range(patient_metadata.shape[0]):
            self.patient_info[i] = patient_metadata.loc[i,'Patient_id']
            # self.patient_driver[patient_metadata.loc[i,'Patient_id']] = patient_metadata.loc[i,'Driver']
        return self

    def run(self,visualize = False):
        """
        main process of Driver2Comm
        :return:
        """
        print('modeling internal and external factor!')
        self.model_internal_and_external_factor(visualize)
        print('Start Associating testing!')
        self.association_test(threshold = self.association_test_threshold)
        print('Finishing Associating testing!')
        return self

    def model_external_factors(self,visualize = False):
        external = External(self.c2c_network,min_support=self.minsup,patient_info = self.patient_info,visualize = visualize)
        print('Start frequent subnetwork mining!')
        external.run()
        print('frequent subnetwork mining finishing!')
        self.frequent_subgraphs = external.get_frequent_pattern()
        self.external_matrix = external.output()
        return self

    def model_internal_factor(self):
        internal_matrix = self.patient_driver_info.iloc[:,1:].T
        internal_matrix.columns = self.patient_driver_info['Patient_id']
        self.internal_matrix = internal_matrix
        return self
    def model_internal_and_external_factor(self,visualize= False):

        self.model_external_factors(visualize)
        self.model_internal_factor()
        #ret = pd.concat([self.external_matrix,self.internal_matrix],axis=0)
        self.output_internal_and_external_matrix()
        return self

    def association_test(self,threshold):
        a = AssociationTest(self.internal_matrix,self.external_matrix,threshold)
        self.association_test_ret = a.association_test()
        return self
    def display_associated_FP(self,save = False):
        for i in range(len(self.association_test_ret)):
            association_test_ret = self.association_test_ret[i]
            internal_factor = association_test_ret['internal']
            tested_FP_list = association_test_ret['idx of passed FP']
            p_value_list = association_test_ret['pvalue of passed FP']
            sorted_idx = np.argsort(p_value_list)
            for i in range(len(sorted_idx)):
                print(internal_factor)
                print('Rank {}'.format(i+1))
                g = self.frequent_subgraphs[tested_FP_list[sorted_idx[i]]].to_graph(gid = tested_FP_list[sorted_idx[i]])
                g.display(output2screen = True)
                print("p values of FP {} : {}".format(tested_FP_list[sorted_idx[i]], p_value_list[sorted_idx[i]]))
                graph_annotation = {'title':(f"{internal_factor}-associated CCC signature {i+1}" ),'P':p_value_list[sorted_idx[i]]}
                g.plot(annotation = graph_annotation,save=save,path=os.path.join(self.outputPATH,internal_factor+'_associated_CCC_signature_'+str(i+1)+'.pdf'))

        return self
    def output_internal_and_external_matrix(self):
        external_name = self.cancer_type+'ExternalSup'+str(self.minsup)+'.csv'
        self.external_matrix.to_csv(os.path.join(self.outputPATH,external_name))
        internal_name = self.cancer_type+'Internal.csv'
        self.internal_matrix.to_csv(os.path.join(self.outputPATH,internal_name))
        return None
    def output_Frequent_subgraph(self,outputPATH = None):
        if outputPATH is None:
            outputPATH = self.outputPATH
        try:
            import codecs
        except Exception as e:
            print('Can not output result: {}'.format(e))
            return
        if not os.path.isdir(outputPATH):
            os.makedirs(outputPATH)
        output_file = codecs.open(os.path.join(outputPATH,'Frequent_Pattern_'+str(self.minsup)+'.txt'),'w','utf-8')
        for gid in range(len(self.frequent_subgraphs)):
            g = self.frequent_subgraphs[gid].to_graph(gid)
            output_file.write(g.display(output2screen=False))
        output_file.close()
        return self

    def output_association_FP(self,outputPATH = None):
        if outputPATH is None:
            outputPATH = self.outputPATH
        try:
            import codecs
        except Exception as e:
            print('Can not output result: {}'.format(e))
            return
        if not os.path.isdir(outputPATH):
            os.makedirs(outputPATH)
        output_file = codecs.open(os.path.join(outputPATH,'AssociationTestResult.txt'),'w','utf-8')
        for i in range(len(self.association_test_ret)):
            association_test_ret = self.association_test_ret[i]
            internal_factor = association_test_ret['internal']
            tested_FP_list = association_test_ret['idx of passed FP']
            p_value_list = association_test_ret['pvalue of passed FP']
            sorted_idx = np.argsort(p_value_list)
            for i in range(len(sorted_idx)):
                output_file.write('internal : {0}\n'.format(internal_factor))
                output_file.write('Rank {}\n'.format(i+1))
                g = self.frequent_subgraphs[tested_FP_list[sorted_idx[i]]].to_graph(gid = tested_FP_list[sorted_idx[i]])
                output_file.write(g.display(output2screen = False))
                output_file.write("adjust p values of FP {} : {}\n".format(tested_FP_list[sorted_idx[i]], p_value_list[sorted_idx[i]]))
        output_file.close()
        return self

    def identify_co_occurrence_Fp(self,association_test_ret:dict):
        internal_vec = association_test_ret['internal vec']
        external_matrix = association_test_ret['external matrix']
        tested_FP_list = association_test_ret['idx of passed FP']
        co_occur_list = []
        for i in range(external_matrix.shape[0]):
            #external_vec = external_matrix.loc[:,tested_FP_list[i]]
            co_occur_cnt = sum((external_matrix.loc[:,tested_FP_list[i]] + internal_vec) == 2)
            sup = sum(external_matrix.loc[:,tested_FP_list[i]])
            if(co_occur_cnt*1.0/sup)>=0.5:
                co_occur_list.append(i)
        return co_occur_list

    def get_candidate_target(self,g,internal):
        return self.get_candidate_target_driver(g,internal)

    def filter_celltype(self,node:str):
        gene,celltype = node.split('__')
        return gene

    def get_candidate_target_driver(self,g,internal_gene):
        candidate_targets = list()
        for vid in g.vertices:
            candidate_targets.append((internal_gene, self.filter_celltype(g.vertices[vid].vlb)))
        return candidate_targets


    def generate_combination_targets_of_an_internal_factor(self,frequent_pattern,gid_list,internal):
        """
        generate all the combination targets of a subtype or an internal gene
        :param frequent_pattern: the DFScode of all the frequent subgraphs
        :param gid_list: list of tested Frequetn subgraph id
        :param internal: the internal factor of combination therapy : can be a subtype or a driver
        :return:
        """
        candidate_targets = list()
        candidate_targets_dict = dict()
        for i in gid_list:
            g = frequent_pattern[i].to_graph(gid=i)
            candidate_targets_g = self.get_candidate_target(g,internal)
            for candidate_target in candidate_targets_g:
                if not candidate_targets_dict.__contains__(candidate_target):
                    candidate_targets_dict[candidate_target] = 1
                    candidate_targets.append(candidate_target)
        return candidate_targets

    def generate_combination_targets(self):
        candidate_targets = collections.defaultdict(list)
        for i in range(len(self.association_test_ret)):
            association_test_ret = self.association_test_ret[i]
            internal_factor = association_test_ret['internal']
            tested_FP_list = association_test_ret['idx of passed FP']
            candidate_targets[internal_factor] = self.generate_combination_targets_of_an_internal_factor(
                self.frequent_subgraphs,
                tested_FP_list,
                internal_factor,
               )
        self.candidate_targets = candidate_targets

    def output_candidated_targets(self):
        try:
            import os
            import codecs
        except Exception as e:
            print("importError with os and codecs")
            return
        outputPATH = self.outputPATH
        candidate_targets = self.candidate_targets
        #check if directory exists
        if not os.path.isdir(outputPATH):
            os.makedirs(outputPATH)
        output_file = codecs.open(os.path.join(outputPATH,'candidate_targets.txt'),'w','utf-8')
        output_file.write('{0} {1}\n'.format('internal', 'exteranl'))
        for _,candidate_targets_internal in candidate_targets.items():
            for target_pair in candidate_targets_internal:
                output_file.write('{0} {1}\n'.format(target_pair[0].upper(),target_pair[1].upper()))
        output_file.close()
        return self

    def __load_lrp_human(self,):
        return pd.read_csv(os.path.join( '../data', 'lrp_human.csv'))

    def formulate_associated_ret(self,associated_file, lrp_human = None):
        if lrp_human is None:
            lrp_human = self.__load_lrp_human()
        node1 = []
        node2 = []
        node1_type = []
        node2_type = []
        is_ct_edge = []
        internal = []  # internal factor
        rank = []
        adjust_p = []
        i = 0
        while (i < len(associated_file)):
            if associated_file[i].startswith('internal'):
                internal_t = associated_file[i].split(' ')[2].rstrip()
                rank_t = int(associated_file[i + 1].split(' ')[1].rstrip())
                nodes_t = {}
                i += 3
                edge_num = 0
                while not associated_file[i].startswith('adjust'):
                    if associated_file[i].startswith('v'):
                        line_split = associated_file[i].split(' ')
                        nodes_t[line_split[1]] = line_split[2].rstrip()
                    elif associated_file[i].startswith('e'):
                        line_split = associated_file[i].split(' ')
                        node1_i, node1_type_i = nodes_t[line_split[1]].split('__')
                        node2_i, node2_type_i = nodes_t[line_split[2]].split('__')
                        if node1_type_i != node2_type_i:
                            if node1_i in lrp_human['ligand'].values:
                                node1.append(node1_i)
                                node1_type.append(node1_type_i)
                                node2.append(node2_i)
                                node2_type.append(node2_type_i)
                            else:
                                node1.append(node2_i)
                                node1_type.append(node2_type_i)
                                node2.append(node1_i)
                                node2_type.append(node1_type_i)
                        else:
                            if node1_i < node2_i:
                                node1.append(node1_i)
                                node1_type.append(node1_type_i)
                                node2.append(node2_i)
                                node2_type.append(node2_type_i)
                            else:
                                node1.append(node2_i)
                                node1_type.append(node2_type_i)
                                node2.append(node1_i)
                                node2_type.append(node1_type_i)
                        if line_split[3].rstrip() == '1':
                            is_ct_edge.append(True)
                        else:
                            is_ct_edge.append(False)
                        rank.append(rank_t)
                        internal.append(internal_t)
                        edge_num += 1
                    i += 1
                adjust_p_t = float(associated_file[i].split(':')[1])
                adjust_p += [adjust_p_t] * edge_num
            i += 1
        associated_mat = pd.DataFrame({'node1': node1, 'node2': node2, 'node1_type': node1_type,
                                       'node2_type': node2_type, 'is_ct_edge': is_ct_edge,
                                       'internal': internal, 'rank': rank, 'adjust_p': adjust_p}, )
        return associated_mat
