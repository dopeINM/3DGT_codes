# -*- coding: utf-8 -*-
"""
Created on Wed Nov  5 17:46:43 2025

@author: domin
"""

# PYTHON Part 4: Calculate assortativity and current-flow betweenness centrality

# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 09:42:25 2023

@author: domin
"""

import os
import networkx as nx
import scipy.io

os.chdir(r'URL') # put here your file location from "TestData_GTA_CentralityAnalysis_Subgraph.mat"

mat = scipy.io.loadmat("TestData_GTA_CentralityAnalysis_Subgraph.mat")
AdjacMat = mat.get("adjac")
G = nx.Graph(AdjacMat)
r = nx.degree_assortativity_coefficient(G)
C = nx.current_flow_betweenness_centrality(G)