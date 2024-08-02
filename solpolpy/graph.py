"""Constructs transformation graph."""
from enum import StrEnum
from inspect import getmembers, isfunction

import astropy.units as u
import networkx as nx

from solpolpy import transforms

System = StrEnum("System", ["bpb", "npol", "stokes", "mzp", "btbr", "bthp", "fourpol", "bp3"])
SYSTEM_REQUIRED_KEYS = {System.bpb: {"B", "pB"},
                        System.npol: set(),
                        System.stokes: {"I", "Q", "U"},
                        System.mzp: {"M", "Z", "P"},
                        System.btbr: {"Bt", "Br"},
                        System.bp3: {"B", "pB", "pBp"},
                        System.bthp: {"B", "theta", "p"},
                        System.fourpol: {str(q) for q in [0.0, 45.0, 90.0, 135.0] * u.degree},
                        }

transform_graph = nx.DiGraph()
polarizer_functions = getmembers(transforms, isfunction)
for function_name, function in polarizer_functions:
    if "_to_" in function_name:
        source, destination = function_name.split("_to_")
        transform_graph.add_edge(System[source.lower()],
                                 System[destination.lower()],
                                 func=function)
