import networkx as nx

from solpolpy.polarizers import bpb_to_mzp, mzp_to_bpb


graph = nx.DiGraph()
graph.add_edge("MZP", "BpB", func=mzp_to_bpb)
graph.add_edge("BpB", "MZP", func=bpb_to_mzp)


if __name__ == "__main__":
    print(graph)
