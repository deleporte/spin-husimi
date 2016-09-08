# We use the DictGraph library and an algorithm by Bollobas to create a 3-regular random graph with n vertices
import numpy as np
import networkx

def randomreggraph(N=100):
    tries = 0
    success = False

    while (success == False):
        graph= networkx.Graph()
        graph.add_nodes_from(range(0,N))
        print graph.number_of_edges()
        tries += 1
        success = True
        # First select 3N points and make a random permutation
        perm = np.random.permutation(3*N)
        # Then glue them together by groups of 3
        
        for i in range(0, 3*N/2):
            vertex1 = perm[2*i]/3
            vertex2 = perm[2*i+1]/3
            if vertex2 in graph.neighbors(vertex1) or (vertex1 == vertex2): #We want a true graph
                success = False
                print "Graph failed because of double edges."
            graph.add_edge(vertex1, vertex2)

        # we then check whether the graph is connected (unlikely at large N)
        if not networkx.is_connected(graph):
            success= False
            print "Graph failed because not connected."

    print "Created graph after", tries, "tries."
    return graph
    
def halftree(N=5):
    if N<1:
        graph=networkx.Graph()
        graph.add_edge(0,1)
        return graph
    else:
        graph=halftree(N-1)
        for i in range(2**(N-1),2**N):
            graph.add_edge(i,2*i)
            graph.add_edge(i,2*i+1)
        return graph

def hexacomb(N=5):
    graph=networkx.Graph()
    for i in range(-N,N):
        for j in range(-N,N):
            graph.add_edge((i,j),(i,j+1))
            if (i+j)%2:
                graph.add_edge((i,j),(i+1,j))
    return graph
