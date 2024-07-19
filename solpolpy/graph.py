"""Constructs transformation graph"""
import networkx as nx

from solpolpy.polarizers import (
    bp3_to_bthp,
    bp3_to_mzp,
    bpb_to_btbr,
    bpb_to_mzp,
    btbr_to_bpb,
    btbr_to_mzp,
    btbr_to_npol,
    fourpol_to_stokes,
    mzp_to_bp3,
    mzp_to_bpb,
    mzp_to_npol,
    mzp_to_stokes,
    npol_to_mzp,
    stokes_to_mzp,
)

transform_graph = nx.DiGraph()
transform_graph.add_edge("npol", "mzp",
                         func=npol_to_mzp,
                         requires_alpha=False)
transform_graph.add_edge("mzp", "bpb",
                         func=mzp_to_bpb,
                         requires_alpha=True)
transform_graph.add_edge("bpb", "mzp",
                         func=bpb_to_mzp,
                         requires_alpha=True)
transform_graph.add_edge("bpb", "btbr",
                         func=bpb_to_btbr,
                         requires_alpha=False)
transform_graph.add_edge("btbr", "bpb",
                         func=btbr_to_bpb,
                         requires_alpha=False)
transform_graph.add_edge("mzp", "stokes",
                         func=mzp_to_stokes,
                         requires_alpha=False)
transform_graph.add_edge("stokes", "mzp",
                         func=stokes_to_mzp,
                         requires_alpha=False)
transform_graph.add_edge("mzp", "bp3",
                         func=mzp_to_bp3,
                         requires_alpha=True)
transform_graph.add_edge("bp3", "mzp",
                         func=bp3_to_mzp,
                         requires_alpha=True)
transform_graph.add_edge("btbr", "mzp",
                         func=btbr_to_mzp,
                         requires_alpha=True)
transform_graph.add_edge("bp3", "bthp",
                         func=bp3_to_bthp,
                         requires_alpha=True)
transform_graph.add_edge("btbr", "npol",
                         func=btbr_to_npol,
                         requires_alpha=True)
transform_graph.add_edge("fourpol", "stokes",
                         func=fourpol_to_stokes,
                         requires_alpha=False)
transform_graph.add_edge("mzp", "npol",
                         func=mzp_to_npol,
                         requires_alpha=False)
