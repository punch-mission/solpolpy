import networkx as nx

from solpolpy.polarizers import npol_to_mzp, bpb_to_mzp, mzp_to_bpb, \
    bpb_to_btbr, btbr_to_bpb, mzp_to_stokes, stokes_to_mzp, \
    mzp_to_bp3, bp3_to_mzp, btbr_to_mzp, bp3_to_bthp, btbr_to_npol, \
    fourpol_to_stokes


transform_graph = nx.DiGraph()
transform_graph.add_edge("npol", "MZP",
                         func=npol_to_mzp,
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
transform_graph.add_edge("BtBr", "BpB",
                         func=btbr_to_bpb,
                         requires_alpha=False)
transform_graph.add_edge("MZP", "Stokes",
                         func=mzp_to_stokes,
                         requires_alpha=False)
transform_graph.add_edge("Stokes", "MZP",
                         func=stokes_to_mzp,
                         requires_alpha=False)
transform_graph.add_edge("MZP", "Bp3",
                         func=mzp_to_bp3,
                         requires_alpha=True)
transform_graph.add_edge("Bp3", "MZP",
                         func=bp3_to_mzp,
                         requires_alpha=True)
transform_graph.add_edge("BtBr", "MZP",
                         func=btbr_to_mzp,
                         requires_alpha=True)
transform_graph.add_edge("Bp3", "Bthp",
                         func=bp3_to_bthp,
                         requires_alpha=True)
transform_graph.add_edge("BtBr", "npol",
                         func=btbr_to_npol,
                         requires_alpha=True)
transform_graph.add_edge("fourpol", "Stokes",
                         func=fourpol_to_stokes,
                         requires_alpha=False)
