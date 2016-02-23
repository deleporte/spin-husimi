from apgl.graph import *

def buildHusimiTree(skeleton)
    """Builds a Husimi tree with given 3-regular skeleton.

    Takes an apgl DictGraph tree as an argument, returns another DictGraph tree"""

    husimi = DictGraph()
    for vertex in skeleton.getAllVertexIds():
        for nb in skeleton.neighbours(vertex):
            husimi.addEdge(*deepcopy(skeleton.neighbours(vertex)).remove(nb))
    # for each neighbour of that vector, we create an edge corresponding to the two others neighbours
    return husimi
