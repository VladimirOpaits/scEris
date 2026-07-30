import networkx as nx
import obonet

CL_URL = "http://purl.obolibrary.org/obo/cl/cl-basic.obo"

TARGETS = [
    ("T_cell", "CL:0000084"),
    ("NK_cell", "CL:0000623"),
    ("ILC", "CL:0001065"),
    ("B_plasma", "CL:0000945"),
    ("dendritic", "CL:0000451"),
    ("mast", "CL:0000097"),
    ("microglia", "CL:0000129"),
    ("myeloid", "CL:0000766"),
    ("erythroid", "CL:0000764"),
    ("megakaryocyte", "CL:0000556"),
    ("epithelial", "CL:0000066"),
    ("endothelial", "CL:0000115"),
    ("pericyte", "CL:0000669"),
    ("fibroblast", "CL:0000057"),
    ("muscle", "CL:0000187"),
    ("neuron", "CL:0000540"),
    ("glia", "CL:0000125"),
    ("neoplastic", "CL:0001063"),
]


def load_graph(path=CL_URL):
    return obonet.read_obo(path)


def build_mapper(graph, targets=TARGETS):
    name2id = {d["name"]: i for i, d in graph.nodes(data=True) if "name" in d}
    order = [t[1] for t in targets]
    label = {t[1]: t[0] for t in targets}

    def compartment(cell_type=None, cl_id=None):
        cid = cl_id if (cl_id and cl_id in graph) else name2id.get(cell_type)
        if cid is None or cid not in graph:
            return "other"
        anc = {cid} | nx.descendants(graph, cid)
        for tid in order:
            if tid in anc:
                return label[tid]
        return "other"

    return compartment
