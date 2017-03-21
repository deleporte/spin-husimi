def buildTwoHex():
    graph= buildPolygon(6)
    graph.add_edge('12','6')
    graph.add_edge('12','26')
    graph.add_edge('26','6')
    graph.add_edge('6','67')
    graph.add_edge('6','7')
    graph.add_edge('67','7')
    graph.add_edge('7','78')
    graph.add_edge('7','8')
    graph.add_edge('78','8')
    graph.add_edge('8','23')
    graph.add_edge('8','28')
    graph.add_edge('28','23')
    return graph

def buildHex():
    return buildPolygon(6)

def buildPolygon(N=4):
    graph = networkx.Graph()
    for k in range(N):
        after=(k+1)%N
        graph.add_edge(str(k),str(after))
        graph.add_edge(str(k),str(k)+str(after))
        graph.add_edge(str(after),str(k)+str(after))
    return graph

def buildHusimiTree(skeleton):
    """Builds a Husimi tree with given 3-regular skeleton.

    Takes an apgl DictGraph tree as an argument, returns another DictGraph tree"""

    husimi = networkx.Graph()
    for vertex in skeleton.nodes():
        if networkx.degree(skeleton,vertex)==3:
            nblist= list(skeleton.neighbors(vertex))
            for nb in list(skeleton.neighbors(vertex)):
                #print nb
                others = [x for x in nblist if x != nb]
                #print others
                edges=[]
                for item in others:
                    if str(item) < str(vertex):
                        edges.append('('+str(item)+','+str(vertex)+')')
                    else:
                        edges.append('('+str(vertex)+','+str(item)+')')
                print edges
                husimi.add_edge(*edges)
        # for each neighbour of that vector, we create an edge corresponding to the two others neighbours. The vertex of the husimi is on the edge of the skeleton, we label it (a,b)
    return husimi

def checkHusimi(graph): #deprecated
    """Checks if given Dictgraph tree is a Husimi subgraph.

    This function seems to take a lot of time..."""
    isHusimi=True
    for edge in graph.getAllEdges():
        n0=graph.neighbours(edge[0])
        n1=graph.neighbours(edge[1])
        if len(n0) != 4 or len(n1) != 4:
            if len(n0)!=2 or len(n1) != 2:
                isHusimi=False
                print "Not enough neighbours !"
        else:
            if len(set(n0).intersection(set(n1))) == 0:
                isHusimi=False
                print "No triangle !"
            if len(set(n0).intersection(set(n1))) > 1:
                print "Many triangles."
    return isHusimi
