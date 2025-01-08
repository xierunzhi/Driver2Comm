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
        :param association_test_threshold : Hyperparameter, the p value threshold of association testing
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
        self.association_test_ret_mat = None
        self.ct_color_mat = {}


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
        self.association_test_ret_mat = self._formulated_associated_FP_mat()
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
        self.output_internal_and_external_matrix()
        return self

    def association_test(self,threshold):
        a = AssociationTest(self.internal_matrix,self.external_matrix,threshold)
        self.association_test_ret = a.association_test()
        return self
    def display_associated_FP(self,save = False):
        """

        :param save:
        :return:
        """
        self._get_ct_color_map()
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
                g.plot(ct_color_mat=self.ct_color_mat,annotation = graph_annotation,save=save,path=os.path.join(self.outputPATH,internal_factor+'_associated_CCC_signature_'+str(i+1)+'.pdf'))

        return self
    def output_internal_and_external_matrix(self):
        external_name = self.cancer_type+'ExternalSup'+str(self.minsup)+'.csv'
        self.external_matrix.to_csv(os.path.join(self.outputPATH,external_name))
        internal_name = self.cancer_type+'Internal.csv'
        self.internal_matrix.to_csv(os.path.join(self.outputPATH,internal_name))
        return None

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
    def _formulated_associated_FP_mat(self,outputPATH = None, lrp_human = None):
        """

        :param outputPATH:
        :param lrp_human: ligand receptor annotation used to sorted edge
        :return:
        """
        if outputPATH is None:
            outputPATH = self.outputPATH
        if lrp_human is None:
            lrp_human = self.__load_lrp_human()
        if not os.path.exists(outputPATH):
            os.makedirs(outputPATH)
        node1 = []
        node2 = []
        node1_type = []
        node2_type = []
        is_ct_edge = []
        internal = []  # internal factor
        rank = []
        adjust_p = []
        for i in range(len(self.association_test_ret)):
            association_test_ret = self.association_test_ret[i]
            internal_factor = association_test_ret['internal']
            tested_FP_list = association_test_ret['idx of passed FP']
            p_value_list = association_test_ret['pvalue of passed FP']
            sorted_idx = np.argsort(p_value_list)
            for i in range(len(sorted_idx)):
                g = self.frequent_subgraphs[tested_FP_list[sorted_idx[i]]].to_graph(gid=tested_FP_list[sorted_idx[i]])
                edge_info = g.display(output2screen=False).rstrip().split('\n')
                nodes_t = {}
                rank_t = i+1
                internal_t = internal_factor
                edge_num = 0

                for j in range(1,len(edge_info)):
                    if edge_info[j].startswith('v'):
                        line_split = edge_info[j].split(' ')
                        nodes_t[line_split[1]] = line_split[2].rstrip()
                    elif edge_info[j].startswith('e'):
                        line_split = edge_info[j].split(' ')
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
                adjust_p_t = p_value_list[sorted_idx[i]]
                adjust_p += [adjust_p_t] * edge_num
        associated_mat = pd.DataFrame({'node1': node1, 'node2': node2, 'node1_type': node1_type,
                                       'node2_type': node2_type, 'is_ct_edge': is_ct_edge,
                                       'internal': internal, 'rank': rank, 'adjust_p': adjust_p}, )
        return associated_mat
    def _get_ct_color_map(self):
        associate_ret_mat = self.association_test_ret_mat
        ct_list = sorted(list(set(associate_ret_mat['node1_type'])|set(associate_ret_mat['node2_type'])))
        import matplotlib as mpl
        colors = mpl.colors.TABLEAU_COLORS
        pal = [color for name, color in colors.items()]
        assert len(ct_list)<=10, "the number of celltype is larger than the maximum setting of Driver2Comm"
        for i,ct in enumerate(ct_list):
            self.ct_color_mat[ct] = pal[i]
        return self
    def output_assocaited_mat(self):
        return self.association_test_ret
    def __load_lrp_human(self):
        current_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(current_dir, '../data', 'lrp_human.csv')
        return pd.read_csv(file_path)


