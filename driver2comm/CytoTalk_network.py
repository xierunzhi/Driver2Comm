VACANT_EDGE_ID = -1
VACANT_VERTEX_ID = -1
VACANT_EDGE_TYPE = -1
VACANT_VERTEX_LABEL = 'NOGENE'
VACANT_GRAPH_ID = -1
AUTO_EDGE_ID = -1
VACANT_EDGE_COST = -1
VACANT_VERTEX_TYPE = -1
from queue import PriorityQueue
#from multiprocessing import queues
import itertools
class Edge(object):
    """Edge class."""

    def __init__(self,
                 eid=VACANT_EDGE_ID,
                 node1=VACANT_VERTEX_ID,
                 node2=VACANT_VERTEX_ID,
                 edge_type=VACANT_EDGE_TYPE,
                 edge_cost= VACANT_EDGE_COST):
        """
        Initialize Edge instance.
        :param eid: edge id
        :param node1: gene 1 of edge
        :param node2: gene 2 of edge
        :param edge_type: edge type
        :param edge_cost: edge cost
        """
        self.eid = eid
        self.node1 = node1
        self.node2 = node2
        self.edge_type = edge_type
        self.edge_cost = edge_cost

class Vertex(object):
    """Vertex class."""

    def __init__(self,
                 vlb=VACANT_VERTEX_LABEL):
        """Initialize Vertex instance.

        Args:
            vid: id of this vertex.
            vlb: label of this veterx :gene name in this case eg. EGFR_Tumor,APOE_Macrophage
        """
        self.vlb = vlb
        self.edges = dict()
        self.edges_queue = PriorityQueue()

    def add_edge(self, eid, node1, node2, edge_type,edge_cost):
        """Add an outgoing edge."""
        self.edges[node2] = Edge(eid, node1, node2, edge_type,edge_cost)
    def fix_edge_cost(self,node2,edge_cost):
            self.edges[node2].edge_cost = (self.edges[node2].edge_cost + edge_cost) / 2
    def set_priority_Queue(self):
        self.edges_queue = PriorityQueue()
        for node2,edge in self.edges.items():
            self.edges_queue.put((self.edges[node2].edge_cost,(self.vlb,node2)))

class c2c_network(object):
    """Graph class."""

    def __init__(self,
                 gid=VACANT_GRAPH_ID,
                 eid_auto_increment=True):
        """Initialize Graph instance.

        Args:
            gid: id of this graph.
            eid_auto_increment: whether to increment edge ids automatically.
        """
        self.gid = gid
        self.vertices = dict()
        self.set_of_elb = set()  # the default elements of the dictionary is a set
        self.set_of_vlb = set()
        self.eid_auto_increment = eid_auto_increment
        self.counter = itertools.count()

    def get_num_vertices(self):
        """Return number of vertices in the graph."""
        return len(self.vertices)

    def add_vertex(self,vlb):
        """Add a vertex to the graph."""
        if vlb in self.vertices:
            return self
        self.vertices[vlb] = Vertex(vlb)
        self.set_of_vlb.add(vlb)
        return self

    def add_edge(self, eid, node1, node2, edge_type,edge_cost):
        """Add an edge to the graph."""
        if (node1 in self.vertices and
                node2 in self.vertices and
                node2 in self.vertices[node1].edges):
            if(edge_cost != self.vertices[node1].edges[node2].edge_cost):
                self.vertices[node1].fix_edge_cost(node2,edge_cost)
                self.vertices[node2].fix_edge_cost(node1,edge_cost)
            return self
        if self.eid_auto_increment:
            eid = next(self.counter)
        self.vertices[node1].add_edge(eid, node1, node2, edge_type,edge_cost)
        self.set_of_elb.add((node1, node2))
        self.vertices[node2].add_edge(eid, node2, node1, edge_type,edge_cost)
        self.set_of_elb.add((node2, node1))
        return self
    def sort_edges(self):
        for _,vertex in self.vertices.items():
            vertex.set_priority_Queue()