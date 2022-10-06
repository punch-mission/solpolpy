import networkx as nx

from solpolpy.polarizers import bpb_to_mzp, mzp_to_bpb, bpb_to_btbr


transform_graph = nx.DiGraph()
transform_graph.add_edge("MZP", "BpB", func=mzp_to_bpb)
transform_graph.add_edge("BpB", "MZP", func=bpb_to_mzp)
transform_graph.add_edge("BpB", "BtBr", func=bpb_to_btbr)


if __name__ == "__main__":
    print(transform_graph)
    print(nx.shortest_path(transform_graph, "MZP", "BtBr"))
    # print(nx.shortest_path(transform_graph, "BtBr", "MZP"))
    print(transform_graph.get_edge_data("MZP", "BpB"))

