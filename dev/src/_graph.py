#!/usr/bin/env python3

import pandas as pd
import numpy as np
import networkx as nx


def get_node_attributes_as_df(nxgraph, filters=dict()):
	"""
	Arguments:
		nxgraph (networkx.Graph): any instance of networkx.Graph
	Returns:
		pd.DataFrame: nodes from the input graph and their attributes as a dataframe
	"""
	node_attributes = pd.DataFrame.from_dict(dict(nxgraph.nodes(data=True)), orient='index')
	for key,val in filters.items(): 
		node_attributes = node_attributes[node_attributes[key] == val]

	return node_attributes


def get_edge_attributes_as_df(nxgraph):
	"""
	Arguments:
		nxgraph (networkx.Graph): any instance of networkx.Graph
	Returns:
		pd.DataFrame: edges from the input graph and their attributes as a dataframe
	"""
	return nx.to_pandas_edgelist(nxgraph, 'protein1', 'protein2')



class Forest: 

	def __init__(self, forest): 

		self.forest = forest

	@property
	def snps(self): 
		return self.node_attributes.index[self.node_attributes["type"] == "SNP"]

	@property
	def proteins(self): 
		return self.node_attributes.index[self.node_attributes["type"] == "protein"]

	@property
	def node_attributes(self): 
		return pd.DataFrame.from_dict(dict(self.forest.nodes(data=True)), orient='index')

	@property
	def edge_attributes(self): 
		return nx.to_pandas_edgelist(self.forest, "protein1", "protein2")

	def get_rich_node_attributes(self): 
		df = pd.concat([
			self.node_attributes, 
			self.degree_by_edge_type()
		], axis=1)
		return df.reindex(self.snps), df.reindex(self.proteins).drop(columns=["chr", "bp", "b", "p"])

	def degree_by_edge_type(self): 
		col1_values = self.edge_attributes.rename(columns={"protein1": "node"}).groupby(["node", "edge"]).size().unstack().fillna(0)
		col2_values = self.edge_attributes.rename(columns={"protein2": "node"}).groupby(["node", "edge"]).size().unstack().fillna(0)
		degrees = col1_values.add(col2_values, fill_value=0).astype(int)
		degrees.columns = [ "degree_" + col for col in degrees.columns ]
		return degrees














