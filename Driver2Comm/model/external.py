from .gSpan import *
from .gSpan import record_timestamp
import collections
import pandas as pd
import numpy as np

class External(gSpan):
    """

    """
    def __init__(self,graphs,min_support,patient_info,visualize = False,output2screen = False):
        super(External, self).__init__(graphs=graphs,min_support=min_support,visualize = visualize,output2screen=output2screen)
        self.subgraph_to_sample = dict()
        self.sample_info = patient_info


    @record_timestamp
    def run(self):
        """Run the gSpan algorithm."""
        self._generate_1edge_frequent_subgraphs()
        if self._max_num_vertices < 2:
            return
        root = collections.defaultdict(Projected)
        for gid, g in self.graphs.items():
            for vid, v in g.vertices.items():
                edges = self._get_forward_root_edges(g, vid)
                for e in edges:
                    root[(v.vlb, e.elb, g.vertices[e.to].vlb)].append(
                        PDFS(gid, e, None)
                    )

        for vevlb, projected in root.items():
            self._DFScode.append(DFSedge(0, 1, vevlb))
            self._subgraph_mining(projected)
            self._DFScode.pop()

    def _report(self, projected):
        sample_info = self.sample_info
        self._frequent_subgraphs.append(copy.copy(self._DFScode))
        if self._DFScode.get_num_vertices() < self._min_num_vertices:
            return
        g = self._DFScode.to_graph(gid=next(self._counter),
                                   is_undirected=self._is_undirected)
        display_str = g.display(self.output2screen)
        display_str += '\nSupport: '+str(self._support)
        if self.output2screen:
            print('\nSupport: {}'.format(self._support))
        if self._visualize:
            g.plot()
        if self._where:
            output = list(set([sample_info[p.gid] for p in projected])) # p.gid means frequent_subgraph belongs to which original graph (patient)
            output_idx = list(set([p.gid for p in projected]))
            if(not self.subgraph_to_sample.__contains__(g.gid)):
                self.subgraph_to_sample[g.gid] = output_idx
            display_str += '\nwhere: {}'.format(output)
            if self.output2screen:
                print('where: {}'.format(output))
        if self.output2screen:
            print('\n-----------------\n')
        # Add some report info to pandas dataframe "self._report_df".
        self._report_df = pd.concat(
            [self._report_df,
             pd.DataFrame(
                 {
                     'support': [self._support],
                     'description': [display_str],
                     'num_vert': self._DFScode.get_num_vertices()
                 },
                 index=[int(repr(self._counter)[6:-1])]
             )
             ],
            axis=0
        )
    def output(self):
        """

        :param patient_path:
        :return: ret is a Dataframe whose shape is (n_fp,n_patient)
        """
        #self.get_sample_name(patient_path)
        nrow = next(self._counter)        #get the num of frequent subgraph
        ncol = len(self.sample_info)
        heatmap = np.zeros((nrow,ncol))
        sample_names = [self.sample_info[key] for key in self.sample_info]
        if len(self.subgraph_to_sample)==0:
            return None
        for key in self.subgraph_to_sample:
            for idx in self.subgraph_to_sample[key]:
                heatmap[key,idx] = 1
        ret = pd.DataFrame(heatmap,columns=sample_names)
        return ret

    def get_frequent_pattern(self):
        return self._frequent_subgraphs