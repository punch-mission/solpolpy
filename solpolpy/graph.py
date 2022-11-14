import networkx as nx

from solpolpy.polarizers import any_to_mzp, bpb_to_mzp, mzp_to_bpb, bpb_to_btbr


transform_graph = nx.DiGraph()
transform_graph.add_edge("any", "MZP",
                         func=any_to_mzp,
                         requires_alpha=False)
transform_graph.add_edge("MZP", "BpB",
                         func=mzp_to_bpb,
                         requires_alpha=True)
transform_graph.add_edge("BpB", "MZP",
                         func=bpb_to_mzp,
                         requires_alpha=True)
transform_graph.add_edge("BpB", "BtBr",
                         func=bpb_to_btbr,
                         requires_alpha=False)
