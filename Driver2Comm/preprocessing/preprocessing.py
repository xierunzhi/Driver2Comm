import os
import pandas as pd
import codecs
from ..model.graph import Graph,AUTO_EDGE_ID
def formulate_c2c_network(input_dir,output_dir,patient_list):
    input_file_path = input_dir
    output_file = os.path.join(output_dir, 'c2c_network.data')
    output = codecs.open(output_file, 'w', 'utf-8')
    count = 0
    for sample_name in patient_list:
        print(sample_name)
        celltype_list = os.listdir(os.path.join(input_file_path, sample_name))
        output.write("t # " + str(count))
        output.write("\n")
        vertex_dict = {}
        vertex_count = 0
        for celltype_name in celltype_list:
            network_path = os.path.join(input_file_path, sample_name, celltype_name, "FinalNetwork.txt")
            if not os.path.exists(network_path):
                continue
            network = pd.read_table(network_path, sep='\t')
            for i in range(network.shape[0]):
                edge = network.iloc[i,]
                node1 = edge['node1'].upper() + '__' + edge['node1_type']
                node2 = edge['node2'].upper() + '__' + edge['node2_type']
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
    return
def read_c2c_network(inputPATH):
        """
        read the formulated c2c_network data
        :return:
        """
        c2c_network = dict()
        #inputPATH = os.path.join(inputPATH,'c2c_network.data')
        with codecs.open(inputPATH, 'r', 'utf-8') as f:
            lines = [line.strip() for line in f.readlines()]
            tgraph, graph_cnt = None, 0
            for i, line in enumerate(lines):
                cols = line.split(' ')
                if cols[0] == 't':
                    if tgraph is not None:
                        c2c_network[graph_cnt] = tgraph
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
                c2c_network[graph_cnt] = tgraph
        return c2c_network
