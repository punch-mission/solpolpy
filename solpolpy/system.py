
from enum import StrEnum

from ndcube import NDCollection

System = StrEnum("System", ["BpB", "angles", "Stokes"])
SYSTEM_REQUIRED_KEYS = {System.BpB: ["B", "pB"],
                        System.Stokes: ["I", "Q", "U"]}


def transform_a(input: NDCollection) -> NDCollection:
    return input


def use_alpha(transform):
    def wrapped(*args, **kwargs) -> None:
        pass
    wrapped.uses_alpha = True
    return wrapped


@use_alpha
def transform_b(input: NDCollection) -> NDCollection:
    assert "alpha" in input
    return input


def check_alpha(func):
    return getattr(func, "uses_alpha", False)


if __name__ == "__main__":
    from inspect import getmembers, isfunction

    import networkx as nx

    from solpolpy import polarizers

    transform_graph = nx.DiGraph()
    polarizer_functions = getmembers(polarizers, isfunction)
    for function_name, function in polarizer_functions:
        if "_to_" in function_name:
            source, destination = function_name.split("_to_")
            transform_graph.add_edge(source, destination, func=function)
