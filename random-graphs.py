# We use the DictGraph library and an algorithm by Bollobas to create a 3-regular random graph with n vertices
from apgl.graph import *
from numpy import *

N=10000 # must be even

tries = 0
success = False

while (success == False):
    graph= DictGraph()
    tries += 1
    success = True
# First select 3N points and make a random permutation
    perm = random.permutation(3*N)
# Then glue them together by groups of 3
    for i in range(0, 3*N/2):
        vertex1 = str(perm[2*i]/3)
        vertex2 = str(perm[2*i+1]/3)
        if graph.edgeExists(vertex1, vertex2) or (vertex1 == vertex2): #We want a true graph
            success = False
            print "Graph failed because of double edges."
        graph.addEdge(vertex1, vertex2)

    # we then check whether the graph is connected (unlikely at large N)
    if (len(graph.depthFirstSearch(0)) != N):
        success = False
        print "Graph failed because not connected"
    print tries

print "Created graph after", tries, "tries."
#print graph.getWeightMatrix()
#print perm    
