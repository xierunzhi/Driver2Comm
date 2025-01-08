import pandas as pd
import numpy as np
import codecs
import os
from .CytoTalk_network  import c2c_network
from .CytoTalk_network import VACANT_EDGE_ID
from .CytoTalk_network import VACANT_VERTEX_LABEL
from .CytoTalk_network import VACANT_EDGE_COST
from queue import PriorityQueue
import collections
def read_graph_from_cytotalk_output(patientpath):
    """
    read communication networks constructed by cytotalk outputs with edge cost
    :param patientpath: PATH that store cytotalk 's output .
    This directory should contain subdirectory named celltypeA-cellTypeB
    :return:
    """
    original_network = c2c_network()
    pcsf_network = c2c_network()
    # read the final network
    celltype_list = os.listdir(patientpath)
    for celltype in celltype_list:
        pcsf_matrix = pd.read_table(os.path.join(patientpath, celltype, 'FinalNetwork.txt'))
        for i in range(pcsf_matrix.shape[0]):
            node1 = pcsf_matrix.loc[i, 'node1'].upper() + '__' + pcsf_matrix.loc[i, 'node1_type']
            node2 = pcsf_matrix.loc[i, 'node2'].upper() + '__' + pcsf_matrix.loc[i, 'node2_type']
            if not pcsf_network.vertices.__contains__(node1):
                pcsf_network.add_vertex(node1)
            if not pcsf_network.vertices.__contains__(node2):
                pcsf_network.add_vertex(node2)
            pcsf_network.add_edge(VACANT_EDGE_ID, node1, node2, int(pcsf_matrix.loc[i, 'is_ct_edge']),
                                  pcsf_matrix.loc[i, 'cost'])
    for celltype in celltype_list:
        original_matrix = pd.read_table(os.path.join(patientpath, celltype, 'IntegratedEdges.txt'))
        for i in range(original_matrix.shape[0]):
            node1 = original_matrix.loc[i, 'node1']
            node2 = original_matrix.loc[i, 'node2']
            gene1, celltype1 = node1.split('__')
            gene2, celltype2 = node2.split('__')
            # update annotation
            node1 = gene1.upper() + '__' + celltype1
            node2 = gene2.upper() + '__' + celltype2
            if not original_network.vertices.__contains__(node1):
                original_network.add_vertex(node1)
            if not original_network.vertices.__contains__(node2):
                original_network.add_vertex(node2)
            original_network.add_edge(VACANT_EDGE_ID, node1, node2, int(celltype1 == celltype2),
                                      original_matrix.loc[i, 'cost'])
    # construct an edge priority
    pcsf_network.sort_edges()
    original_network.sort_edges()
    return pcsf_network, original_network
def _inequal_celltype(node1, node2):
    gene1, ct1 = node1.split('__')
    gene2, ct2 = node2.split('__')
    if ct1 == ct2:
        ret = True
    else:
        ret = False
    return ret
def get_ie_pathway(internal_gene, external_gene, pcsf_network, original_network):
    """
    identify shortest pathway between internal_gene and external_gene in extend_network
    :param internal_gene:
    :param external_gene:
    :param pcsf_network:
    :param original_network:
    :return:
    """
    if internal_gene not in original_network.set_of_vlb or external_gene not in pcsf_network.set_of_vlb:
        return pd.DataFrame()
    internal_network = c2c_network()
    internal_network.add_vertex(internal_gene)
    if (internal_gene in pcsf_network.set_of_vlb):
        internal_network = _combine_network(internal_network, pcsf_network)
    else:
        original_network.sort_edges()
        internal_network = _extend_network(pcsf_network.set_of_vlb, internal_network, original_network)
        internal_network = _combine_network(internal_network, pcsf_network)
    external_network = c2c_network()
    if (external_gene in pcsf_network.set_of_vlb):
        external_network = _combine_network(internal_network, external_network)
    else:
        external_network.add_vertex(external_gene)
        original_network.sort_edges()
        external_network = _extend_network(pcsf_network.set_of_vlb, external_network, original_network)
    combined_network = _combine_network(internal_network, external_network)
    ret = _identify_ie_pathway_util(internal_gene, external_gene, combined_network)
    return ret
def _extend_network(geneset, internal_network, original_network):
    while len(geneset & internal_network.set_of_vlb) < 1:
        internal_geneset = internal_network.set_of_vlb.copy()
        for node1 in internal_geneset:
            while not original_network.vertices[node1].edges_queue.empty():
                cost, (_, node2) = original_network.vertices[
                    node1].edges_queue.get()  # get a shortest edge connected to node1
                if node2 not in internal_geneset:
                    internal_network.add_vertex(node2)
                    internal_network.add_edge(VACANT_EDGE_ID, node1, node2,
                                              original_network.vertices[node1].edges[node2].edge_type,
                                              original_network.vertices[node1].edges[node2].edge_cost)
    # once jump out of circle, we find the intersect part between internal and external
    return internal_network

def _combine_network(internal_network: c2c_network, external_network):
    # combined the internal network and external network to internal network
    combined_network = internal_network
    for label, vertex in external_network.vertices.items():
        if not combined_network.vertices.__contains__(label):
            combined_network.add_vertex(label)
        for label2, edge in vertex.edges.items():
            if not combined_network.vertices.__contains__(label2):
                combined_network.add_vertex(label2)
            combined_network.add_edge(VACANT_EDGE_ID, label, label2, edge.edge_type, edge.edge_cost)
    return combined_network
def _get_shortest_path(internal_gene, exteranl_node, prior_vertex: dict, combination_network):
    shortest_path = collections.defaultdict(list)
    node1 = prior_vertex[exteranl_node]
    node2 = exteranl_node
    while node2 != internal_gene:
        gene1,ct1 = node1.split('__')
        gene2, ct2 = node2.split('__')
        shortest_path['node1'].append(gene1)
        shortest_path['node2'].append(gene2)
        shortest_path['node1_type'].append(ct1)
        shortest_path['node2_type'].append(ct2)
        shortest_path['is_ct_edge'].append(bool(combination_network.vertices[node1].edges[node2].edge_type))
        shortest_path['cost'].append(combination_network.vertices[node1].edges[node2].edge_cost)
        if (node2 == exteranl_node):
            shortest_path['edge_role'].append('end')
        elif (node1 == internal_gene):
            shortest_path['edge_role'].append('start')
        else:
            shortest_path['edge_role'].append('intermediate')
        node2 = node1
        node1 = prior_vertex[node2]
    shortest_path_matrix = pd.DataFrame(shortest_path)
    shortest_path_matrix = shortest_path_matrix.reindex(index=shortest_path_matrix.index[::-1])
    return shortest_path_matrix
def _identify_ie_pathway_util(internal_gene, exteranl_gene, combination_network):
    '''
    utils:identify the shorest path from internal gene to external gene
    :param internal_gene:
    :param exteranl_gene:
    :param combination_network:
    :return:
    '''
    set_of_searched_vertices = set()
    set_of_unsearched_vertices = set()
    combination_network.sort_edges()
    set_of_searched_vertices.add(internal_gene)
    distance = {internal_gene: 0}
    prior_vertex = {internal_gene: VACANT_VERTEX_LABEL}
    closest_vertex_queue = PriorityQueue()
    for vlb, vertex in combination_network.vertices.items():
        if vlb != internal_gene:
            set_of_unsearched_vertices.add(vlb)
            distance[vlb] = float('inf')
            prior_vertex[vlb] = None
    while not combination_network.vertices[internal_gene].edges_queue.empty():
        cost, (node1, node2) = combination_network.vertices[internal_gene].edges_queue.get()
        distance[node2] = cost
        prior_vertex[node2] = node1
        closest_vertex_queue.put((cost, (node1, node2)))
    # start searching shortest path from internal to external
    while (not closest_vertex_queue.empty()) and (exteranl_gene in set_of_unsearched_vertices):
        _, (_, node2) = closest_vertex_queue.get()
        if node2 not in set_of_searched_vertices:
            set_of_searched_vertices.add(node2)
            set_of_unsearched_vertices.remove(node2)
            # update distance dict and prior dict
            while not combination_network.vertices[node2].edges_queue.empty():
                cost, (_, node3) = combination_network.vertices[node2].edges_queue.get()
                if node3 not in set_of_searched_vertices:
                    if cost + distance[node2] < distance[node3]:
                        distance[node3] = cost + distance[node2]
                        prior_vertex[node3] = node2
                    closest_vertex_queue.put((distance[node3], (node2, node3)))
    shortest_path_matrix =_get_shortest_path(internal_gene, exteranl_gene, prior_vertex,
                                                  combination_network)
    return shortest_path_matrix
def _get_frequent_subgraph_matrix(k,association_test_ret,frequent_subgraphs):
    """
    generate edge for frequent subgraph as a dataframe format and then concat to the shortest path
    :param k:
    :param association_test_ret:  assoicated result of interested Driver
    :return:
    """
    set_of_edge = set()
    fp_info = collections.defaultdict(list)
    for i in range(k):
        fpid = association_test_ret['idx of passed FP'][i]
        g = frequent_subgraphs[fpid].to_graph(gid = fpid)
        for vid,vertex in g.vertices.items():
            for toid in vertex.edges.keys():
                node1 = vertex.vlb
                node2 = g.vertices[toid].vlb
                if (node1,node2) not in set_of_edge:
                    set_of_edge.add((node1,node2))
                    set_of_edge.add((node2,node1))
                    gene1, ct1 = node1.split('__')
                    gene2, ct2 = node2.split('__')
                    fp_info['node1'].append(gene1)
                    fp_info['node2'].append(gene2)
                    fp_info['node1_type'].append(ct1)
                    fp_info['node2_type'].append(ct2)
                    fp_info['is_ct_edge'].append(ct1!=ct2)
                    fp_info['cost'].append(VACANT_EDGE_COST)
                    fp_info['edge_role'].append('frequent pathway')
    return  pd.DataFrame(fp_info)
def generate_top_k_ie_pathway(association_test_ret,frequent_subgraphs,internal_node, k, patient_path):
    """

    :param internal_gene:
    :param k:
    :param patient_path:
    :param association_test_ret:
    :return:
    """
    print(patient_path)
    driver,ct = internal_node.split('__')
    for i in range(len(association_test_ret)):
        if association_test_ret[i]['internal'] == driver:
            association_test_ret = association_test_ret[i]
            break

    tested_FP_list = association_test_ret['idx of passed FP']
    shortest_path_network = pd.DataFrame()
    pcsf_network,original_network = read_graph_from_cytotalk_output(patient_path)
    for i in range(k):
        # if i == 9:
        #     print(1)
        g = frequent_subgraphs[tested_FP_list[i]].to_graph(
            gid=tested_FP_list[i])
        for vid, vertex in g.vertices.items():
            shortest_path_matrix = get_ie_pathway(internal_node, vertex.vlb, pcsf_network,
                                                          original_network)
            if shortest_path_matrix.shape[0] > 0:
                print(f'identifying pathway to {vertex.vlb} in {driver}-associated CCC signature {i + 1} ')
            shortest_path_network = pd.concat((shortest_path_network, shortest_path_matrix), axis=0)
    frequent_subgraph_matrix = _get_frequent_subgraph_matrix(k,association_test_ret,frequent_subgraphs)
    shortest_path_network = pd.concat((frequent_subgraph_matrix,shortest_path_network),axis=0)
    shortest_path_network = shortest_path_network.reset_index(drop=True)
    return shortest_path_network
