"""
 Driver2Comm is designed for identifying cell-type communication patterns associated with genetic variation in cancer.

"""

import pandas as pd
import numpy as np
import os
from graph import *
from external import External
from association_test import AssociationTest
import codecs
import collections

class Driver2Comm(object):

    def __init__(self,minsup,inputPATH,patient_metadata:pd.DataFrame,num_celltype=3,mode = 'driver',outputPATH = './result',
                 subtype_annotation = '/input/brca_annotation.csv',visualize = False,cancer_type = 'Lung'):

        self.minsup = minsup
        self.num_celltype = num_celltype
        self.outputPATH = outputPATH
        self.inputPATH = inputPATH
        self.c2c_network = dict()
        self.patient_info = dict()
        self.external_matrix = None
        self.internal_matrix = None
        self.patient_driver = dict()        #patient_driver[patient name] = patient driver
        self.frequent_subgraphs = None
        self.candidate_targets = None
        self.drivers = dict()       #example: drivers[BRAF] = 0
        self.mode = mode
        self.association_test_ret = None
        self._visualize = visualize
        self.cancer_type = cancer_type
        if mode == 'subtype':
            self.subtype2gene = self.get_subtype_gene(subtype_annotation)
        if not os.path.exists(self.outputPATH):
            os.makedirs(self.outputPATH)

        self.get_patient_info(patient_metadata)
        self.formulate_c2c_network()
        self.read_c2c_network()

    def read_c2c_network(self):
        """
        read the formulated c2c_network data
        :return:
        """
        inputPATH = os.path.join(self.outputPATH,'c2c_network.data')
        with codecs.open(inputPATH, 'r', 'utf-8') as f:
            lines = [line.strip() for line in f.readlines()]
            tgraph, graph_cnt = None, 0
            for i, line in enumerate(lines):
                cols = line.split(' ')
                if cols[0] == 't':
                    if tgraph is not None:
                        self.c2c_network[graph_cnt] = tgraph
                        graph_cnt += 1
                        tgraph = None
                    if cols[-1] == '-1':
                        break
                    tgraph = Graph(graph_cnt)
                elif cols[0] == 'v':
                    tgraph.add_vertex(cols[1], cols[2])
                elif cols[0] == 'e':
                    tgraph.add_edge(AUTO_EDGE_ID, cols[1], cols[2], cols[3])
            # adapt to input files that do not end with 't # -1'
            if tgraph is not None:
                self.c2c_network[graph_cnt] = tgraph
        return self
    def formulate_c2c_network(self):
        input_file_path = self.inputPATH
        dir_list = os.listdir(input_file_path)
        dir_list.sort()
        output_file = os.path.join(self.outputPATH,'c2c_network.data')
        output = codecs.open(output_file, 'w', 'utf-8')
        count = 0
        for sample_name in dir_list:
            celltype_list = os.listdir(os.path.join(input_file_path, sample_name))
            output.write("t # " + str(count))
            output.write("\n")
            vertex_dict = {}
            vertex_count = 0
            for celltype_name in celltype_list:
                network_path = os.path.join(input_file_path, sample_name, celltype_name, "FinalNetwork.txt")
                network = pd.read_table(network_path, sep='\t')
                for i in range(network.shape[0]):
                    edge = network.iloc[i,]
                    node1 = edge['node1'] + '__' + edge['node1_type']
                    node2 = edge['node2'] + '__' + edge['node2_type']
                    if (not vertex_dict.__contains__(node1)):
                        vertex_dict[node1] = vertex_count
                        output.write('v {0} {1}'.format(str(vertex_count), node1))
                        output.write('\n')
                        vertex_count += 1
                    if (not vertex_dict.__contains__(node2)):
                        vertex_dict[node2] = vertex_count
                        output.write('v {0} {1}'.format(str(vertex_count), node2))
                        output.write('\n')
                        vertex_count += 1
                    output.write(
                        'e {0} {1} {2}'.format(vertex_dict[node1], vertex_dict[node2], int(edge['is_ct_edge'])))
                    output.write('\n')
                    # finish graph edge writing
            count = count + 1  # turn to the next graph
        # after all graphs done,write the end marker
        output.write('t # -1')
        output.close()
        return self


    def get_patient_info(self,patient_metadata):
        for i in range(patient_metadata.shape[0]):
            self.patient_info[i] = patient_metadata.iloc[i,0]
            self.patient_driver[patient_metadata.iloc[i,0]] = patient_metadata.iloc[i,1]
        return self

    def model_external_factors(self):
        external = External(self.c2c_network,min_support=self.minsup,patient_info = self.patient_info,visualize = self._visualize)
        external.run()
        self.frequent_subgraphs = external.get_frequent_pattern()
        self.external_matrix = external.output()
        return self

    def model_internal_factor(self):
        patients = [patient for patient in self.patient_driver.keys()]
        patients.sort()
        n_patient = len(patients)
        driver_set = list(set([driver for driver in self.patient_driver.values()]))
        driver_set.sort()
        for i in range(len(driver_set)):
            self.drivers[driver_set[i]] = i
        n_driver = len(driver_set)
        internal_dataframe = np.zeros((n_driver, n_patient))
        for patient_name, patient_driver in self.patient_driver.items():
            internal_dataframe[driver_set.index(patient_driver), patients.index(patient_name)] = 1
        # print(internal_dataframe)
        self.internal_matrix = pd.DataFrame(internal_dataframe, index=driver_set, columns=patients)
        return self
    def model_internal_and_external_factor(self):
        self.model_external_factors()
        self.model_internal_factor()
        #ret = pd.concat([self.external_matrix,self.internal_matrix],axis=0)
        self.output_internal_and_external_matrix()
        return self

    def run(self):
        """
        main process of Driver2Comm
        :return:
        """
        self.model_internal_and_external_factor()
        self.association_test(threshold = 0.05)
        self.generate_combination_targets()

    def association_test(self,threshold):
        a = AssociationTest(self.internal_matrix,self.external_matrix,threshold)
        self.association_test_ret = a.association_test()
        #self.display_associated_FP()
        return self
    def display_associated_FP(self):
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
                g.display()
                print("p values of FP {} : {}".format(tested_FP_list[sorted_idx[i]], p_value_list[sorted_idx[i]]))
                graph_annotation = {'title':('Communication signature '+ str(i+1)),'P':p_value_list[sorted_idx[i]]}
                g.plot(annotation = graph_annotation,save=True,path=os.path.join('./output/result/',internal_factor+str(i+1)+'.pdf'))
                #candidate_targets_g = self.get_candidate_target(g,self.external_matrix,self.drivers,self.internal_matrix)
                #for candidate_target in candidate_targets_g:
                #    self.candidate_targets.append(candidate_target)
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
        output_file.write('Mode : {0} \n'.format(self.mode.capitalize()))
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
                output_file.write("p values of FP {} : {}\n".format(tested_FP_list[sorted_idx[i]], p_value_list[sorted_idx[i]]))
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


    def get_subtype_gene(self,annotatedPATH):
        subtype_gene_matrix = pd.read_csv(annotatedPATH)
        subtype2gene = collections.defaultdict(list)
        for i in range(subtype_gene_matrix.shape[0]):
            subtype2gene[subtype_gene_matrix.iloc[i,0]].append(subtype_gene_matrix.iloc[i,1])
        return subtype2gene
    def get_candidate_target(self,g,internal,mode = 'driver'):
        if mode == 'driver':
            return self.get_candidate_target_driver(g,internal)
        elif mode == 'subtype':
            return self.get_candidate_target_subtype(g,internal)
        else:
            raise Exception('Please correct the mode to driver or subtype')

    def filter_celltype(self,gene:str):
        if 'Macrophages' in gene:
            return gene.replace('__Macrophages','')
        elif 'Cd8+Tcells' in gene:
            return gene.replace('__Cd8+Tcells','')
        elif 'Tumor' in gene:
            return gene.replace('__Tumor','')
        else:
            raise Exception('Something wrong in cell type')
    def get_candidate_target_driver(self,g,internal_gene):
        candidate_targets = list()
        for vid in g.vertices:
            candidate_targets.append((internal_gene, self.filter_celltype(g.vertices[vid].vlb)))
        return candidate_targets

    def get_candidate_target_subtype(self,g,internal):
        internal_gene_list = self.subtype2gene[internal]
        candidate_targets = list()
        for vid in g.vertices:
            for internal_gene in internal_gene_list:
                if internal_gene != 'NOGENE':
                    candidate_targets.append((internal_gene,self.filter_celltype(g.vertices[vid].vlb)))
        return candidate_targets

    def generate_combination_targets_of_an_internal_factor(self,frequent_pattern,gid_list,internal,mode = 'driver'):
        """
        generate all the combination targets of a subtype or an internal gene
        :param frequent_pattern: the DFScode of all the frequent subgraphs
        :param gid_list: list of tested Frequetn subgraph id
        :param internal: the internal factor of combination therapy : can be a subtype or a driver
        :param mode: driver or subtype
        :return:
        """
        candidate_targets = list()
        candidate_targets_dict = dict()
        for i in gid_list:
            g = frequent_pattern[i].to_graph(gid=i)
            candidate_targets_g = self.get_candidate_target(g,internal,mode)
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
                self.mode)
        self.candidate_targets = candidate_targets

    def output_candidated_targets(self,outputPATH,candidate_targets:dict):
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



