import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

class MERA:

    def _system_len(self, fname):
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i

    def _get_cut_list(self):
        cut_list = [0]*(self.L-1)
        for U in self.U_list:
            for i in range(U[0], U[2]):
                cut_list[i] += 1
        return cut_list

    def _get_levels(self):
        level_list = {i+1:[] for i in range(self.depth)}
        for idx in range(self.L):
            U = "U{}".format(idx)
            U_conns = [self.G.nodes["U{}".format(j)]['level'] for j in range(idx) \
                    if "U{}".format(j) in self.shortest_paths[U].keys() ]
            level = 1 if len(U_conns) == 0 else max(U_conns) + 1
            self.G.nodes[U]['level'] = level
            self.G.nodes[idx]['level'] = 0
            level_list[level].append(U)
        return level_list

    def __init__(self, fname):

        # create database
        L = self._system_len(fname)
        U_list = [None]*L
        p_list = [None]*L

        # read in unitaries from file
        fp = open(fname, 'r')
        fp.readline()
        i = 0
        for line in fp:
            data = list(map(int, line.split()[:-1]))
            U_list[i] = data
            i += 1

        # construct graph
        G = nx.DiGraph()
        for i in range(L):
            G.add_node(i)
            p_list[i] = i

        for k, U in enumerate(U_list):
            lbl = "U" + str(k)
            G.add_node(lbl)
            for i in range(U[0], U[2] + 1):
                if not p_list[i] == None:
                    G.add_edge(lbl, p_list[i])
                p_list[i] = lbl
            p_list[U[1]] = None

        self.G = G
        self.U_list = U_list
        self.p_list = p_list
        self.L = L
        self.depth = nx.dag_longest_path_length(self.G)
        self.cut_list = self._get_cut_list()
        self.shortest_paths =  dict(nx.all_pairs_shortest_path(G))
        self.level_list = self._get_levels()

    def draw_graph(self, func, **kwargs):
        plt.figure()
        func(self.G, **kwargs)

    def unitaries_cut(self, i):
        '''
        computes the number of unitaries cut
        if the the system is cut between site
        i and i + 1
        '''
        return self.cut_list[i]

    def get_path(self, i, j):
        G = nx.Graph(self.G)
        return self.shortest_paths[i][j]

    def get_largest_level_from_path(self, i, j):
        i_set = nx.dag.ancestors(self.G, i)
        j_set = nx.dag.ancestors(self.G, j)
        intersect = i_set.intersection(j_set)
        return min([self.G.nodes[U]['level'] for U in intersect])

    def get_depth(self):
        return self.depth

    def get_level_map(self):
        level_map = np.zeros((self.L, self.L))
        for i in range(self.L):
            for j in range(self.L):
                level_map[i,j] = self.get_largest_level_from_path(i,j)
        return level_map

    def get_mean_cut_number(self):
        return np.mean(self.cut_list)

    def get_avg_dist_depth(self, d):
        depths = [self.get_largest_level_from_path(i, i + d) for i in range(self.L - d)]
        return np.mean(depths)

    def get_avg_dist_depth_list(self):
        return [self.get_avg_dist_depth(i) for i in range(self.L)]

    def get_balance_metric(self):
        shortest_path = min([len(self.shortest_paths["U{}".format(self.L - 1)][i]) for i in range(self.L)])
        return self.depth - shortest_path



