from apgl.graph import *

def buildHusimiTree(skeleton)
    """Builds a Husimi tree with given 3-regular skeleton.

    Takes an apgl DictGraph tree as an argument, returns another DictGraph tree"""

    husimi = apgl.graph.DictGraph()
    for vertex in skeleton.getAllVertexIds():
        for nb in skeleton.neighbours(vertex):
            others=deepcopy(skeleton.neighbours(vertex)).remove(nb)
            edges=[]
            for item in others:
                if str(item) < str(vertex):
                    edges.append('('+str(item)+','+str(vertex)+')')
                else
                    edges.append('('+str(vertex)+','+str(item)+')')
            husimi.addEdge(*edges)
    # for each neighbour of that vector, we create an edge corresponding to the two others neighbours. The vertex of the husimi is on the edge of the skeleton, we label it (a,b)
    return husimi
