import numpy as np
import apgl.graph
import scipy.linalg
import scipy.sparse
import random

# to each vertex of the graph are attached five values
# the first value is a boolean to see in which hemisphere our spin vector is
# the second value is a complex double denoting the stereographical projection
# the third is a complex double denoting the holomorphic derivative of the hamiltonian
# the last two are complex doubles for the hessian wrt variables on the same sphere (dzdz and dzdbarz)

# besides, to each edge of the graph is attached two complex doubles, denoting the various second derivatives involving two variables linked by an edge

def getpositions(graph)
    """Returns the list of the vertex positions in the graph"""
    positions=[]
    for vertex in graph.getAllvertexIds():
        upz,z,__=graph.getvertex(vertex)
        if upz==0:
            positions.append(z)
        else:
            positions.append(4./z)
    return positions

def under(z,w,powz,poww)
    return (1.+abs(z)**2./4.)**powz*(1.+abs(w)**2./4.)**poww

def jet1(graph):
    """computes the hamiltonian and its gradient, for a given position (aka first three values attached do the edge)

    Its input is a graph with two labels on each vertex. The function transforms the graph to add a third label on each vertex, then returns the value of the hamiltonian.

    If more than two labels are assigned to a vertex beforehand, this function will erase them."""
    ham = 0
    for edge in graph.getAllEdge():
        upz,z,*__ = graph.getvertex(edge[0])
        upw,w,*__ = graph.getvertex(edge[1])
        zbar=np.conjugate(z)
        wbar=np.conjugate(w)
        if upz == upw:
            ham += 1.5-abs(z-w)**2./(2.*under(z,w,1,1))
            dz = (wbar-zbar)*(1.+zbar*w/4.)/(2.*under(z,w,2,1))
            dw = (zbar-wbar)*(1.+wbar*z/4.)/(2.*under(z,w,1,2))
        else:
            ham += 1.5 - abs(z*w/4.-1.)**2./(2.*under(z,w,1,1))
            dz = -(zbar*wbar/4.-1.)*(zbar+w)/(2.*under(z,w,2,1))
            dw = -(wbar*zbar/4.-1.)*(wbar+z)/(2.*under(z,w,1,2))
        graph.setvertex(edge[0])[upz,z,dz]
        graph.setvertex(edge[1])[upw,w,dw]

    return ham

def jet2(graph):
    """computes the second jet, for a given position (aka first three values attached do the edge)

    Its input is a graph with three labels on each vertex. The function transforms the graph to add two labels on each vertex and two on each edge.

    If more than three labels are assigned to a vertex beforehand, this function will erase them."""
    
    # The computations should be checked
    for edge in graph.getAllEdges():
        upz,z,dz,*__ = graph.getvertex(edge[0])
        upw,w,dw,*__ = graph.getvertex(edge[1])
        zbar=np.conjugate(z)
        wbar=np.conjugate(w)
        if upz == upw:
            dzdz = zbar/4.*(zbar-wbar))*(1.+w*zbar)/4.)/under(z,w,3,1)
            dzdbarz = -((1.-abs(w)**2./4.)*(1.-abs(z)**2./4.)+(z*wbar).real)/(2.*under(z,w,3,1))
            dwdw = wbar/4.*(wbar-zbar)*(1.+z*wbar/4.)/under(z,w,1,3)
            dwdbarw = -((1.-abs(z)**2./4.)*(1.-abs(w)**2./4.)+(w*zbar).real)/under(z,w,1,3)
            dzdw = -(zbar-wbar)**2./(8.*under(z,w,2,2))
            dzdbarw = (1+zbar*w/4.)**2/(2.*under(z,w,2,2))
        else:
            dzdz = zbar/4.*(zbar*wbar/4.-1.)*(w+zbar)/under(z,w,3,1)
            dzdbarz = ((1.-abs(w)**2./4.)*(1.-abs(z)**2./4.)-(z*w).real)/(2.*under(z,w,3,1))
            dwdw = wbar/4.*(wbar*zbar/4.-1.)*(z+wbar)/under(z,w,1,3)
            dwdbarw = ((1.-abs(z)**2./4.)*(1.-abs(w)**2./4.)-(w*z).real)/(2.*under(z,w,1,3))
            dzdw = -(1.+zbar*wbar/4.)**2./(2.*under(z,w,2,2))
            dzdbarw = (w+zbar)**2./(8.*under(z,w,2,2))
        graph.setvertex(edge[0],[upz,z,dz,dzdz,dzdbarz])
        graph.setvertex(edge[1],[upw,w,dw, dwdw, dwdbarw])
        graph.addEdge(edge[0],edge[1],[dzdw,dzdbarw]) #the addEdge method erases the previous label (undocumented)

def hess(graph)
    """computes the Hessian with given 2-jet"""

    dbgraph=apgl.graph.DictGraph(False)
    for vertex in graph.getAllVertexIds():
        __,__,__,dzdz,dzdbarz=graph.getVertex(vertex)
        dbgraph.addEdge(str(vertex)+'x',str(vertex)+'x',dzdbarz+dzdz.real)
        dbgraph.addEdge(str(vertex)+'y',str(vertex)+'y',dzdbarz-dzdz.real)
        dbgraph.addEdge(str(vertex)+'x',str(vertex)+'y',-dzdz.imag)
        dbgraph.addEdge(str(vertex)+'y',str(vertex)+'x',-dzdz.imag)
    for e in graph.getAllEdges():
        dzdw,dzdbarw=graph.getEdge(*edge)
        dbgraph.addEdge(str(e[0])+'x',str(e[1])+'x',(dzdbarw+dzdw).real)
        dbgraph.addEdge(str(e[0])+'x',str(e[1])+'y',(dzdbarw+dzdw).imag)
        dbgraph.addEdge(str(e[0])+'y',str(e[1])+'x',(dzdw-dzdbarw).imag)
        dbgraph.addEdge(str(e[0])+'y',str(e[1])+'y',(dzdbarw-dzdw).real)
        dbgraph.addEdge(str(e[1])+'x',str(e[0])+'x',(dzdbarz+dzdw).real)
        dbgraph.addEdge(str(e[1])+'x',str(e[0])+'y',(dzdw-dzdbarw).imag)
        dbgraph.addEdge(str(e[1])+'y',str(e[0])+'x',(dzdbarz+dzdw).imag)
        dbgraph.addEdge(str(e[1])+'y',str(e[0])+'y',(dzdbarw-dzdw).real)
 
    #implementing the (absent) weighted graph 
    Jhessian=np.matrix(np.zeros((dbgraph.getNumVertices(),dbgraph.getNumVertices())))
    i=-1
    for vertex1 in dbgraph.getAllVertexIds():
        i++
        j=-1
        for vertex2 in dbgraph.getAllVertexIds():
            j++
            if dbgraph.edgeExists(vertex1,vertex2):
                Jhessian[i,j]=dbgraph.getEdge(vertex1,vertex2)
    return J*Jhessian

def charac(hessian)
    """Computes the characteristic, given the hessian"""
    J=np.matrix(np.zeros(hessian.shape))
    for i in range(0.5*sqrt(hessian.shape[1])):
        J[2*i,2*i+1]=-1
        J[2*i+1,2*i]=1
    
    eigs=scipy.linalg.eigvals(J*hessian).imag #We want ALL the eigenvalues

    #There must be a better python way to do the following :
    mu=0
    for value in eigs:
        if value >0:
            mu += value

    return mu

def initialize(graph)
    """Initializes the spins at random positions. Flushes the labels."""
    for vertex in graph.getAllVertexIds():
        graph.setVertex(vertex, [random.randint(0,1),random.uniform(-2,2)+random.uniform(-2,2)j])

def find_well(graph)
    """Uses a gradient descent to find a local minimum of the hamiltonian.

    Returns a boolean who tells whether the algorithm converged or not."""
    THRESHOLD = 0.01
    DESCENT = 0.1
    STOP = 0.000001
    ham = jet1(graph)
    N=graph.getNumVertices()
    oldham = ham + graph.getNumVertices()
    while (ham > THRESHOLD*N) and (abs(ham-oldham) > STOP*N):
        for vertex in graph.getAllVertexIds():
            upz,z,dz = graph.getVertex(vertex)
            z -= DESCENT*np.conjugate(dz)
            if abs(z) > 4:
                upz=1.-upz
                z=4./z
            graph.setVertex(vertex,[upz,z])
        oldham,ham = ham,jet1(graph)
    if ham > THRESHOLD*N:
        print "Algorithm stopped at local minimum ", ham
        return False
    else:
        return True

def find_miniwell(graph)
    """uses two gradient descents to find the minimum of the characteristics among points where the hamiltonian is minimal """
    random.seed()
    at_minimum = False
    THRESHOLD=0.01
    DIFF_STEP=0.0001
    DESCENT=0.1
    N=graph.getNumVertices()
    while not at_minimum:
        print "Searching for minimum..."
        initialize(graph)
        at_minimum=find_well(graph)
    print "Found minimum."
    jet2(graph)
    at_miniwell = False
    while not at_miniwell:
        hessian = hess(graph)
        values,vectors = scipy.linalg.eigvals(hessian)
        mu = charac(hessian)
        smallvectors=[]
        diff=[]
        for index , value in enumerate(values):
            if value < sqrt(THRESHOLD):
                smallvectors.append(vectors[index])
        dim=smallvectors.size
        if  dim==0:
            print "Found well."
        else:
            print "Found submanifold of dimension ", dim,"."
            print "Searching for miniwell..."
        at_miniwell = True
        #first we compute the derivative of mu wrt each "zero" direction
        for tangent in smallvectors:
            for i,vertex in enumerate(graph.getAllVertexIds()):
                upz,z,__=graph.getVertex(vertex)
                z += DIFF_STEP*(vectors[2*i]+vectors[2*i+1]j)
                graph.setVertex(vertex,[upz,z])
            jet1(graph)
            jet2(graph)
            muvar = charac(hess(graph))
            if abs(muvar-mu)/mu > sqrt(THRESHOLD)*DIFF_STEP*(N-dim):
                at_miniwell = False
            diff.append((muvar-mu)/(DIFF_STEP*mu)))
            for i,vertex in enumerate(graph.getAllVertexIds()):
                upz,z,__=graph.getVertex(vertex)
                z -= DIFF_STEP*(vectors[2*i]+vectors[2*i+1]j)
                graph.setVertex(vertex,[upz,z])
        #then we move in the good direction
        for tangent in smallvectors:
            for i,vertex in enumerate(graph.getAllVertexIds()):
                upz,z,__=graph.getVertex(vertex)
                z -= diff[index]*(vectors[2*i]+vectors[2*i+1]j)*DESCENT
                graph.setVertex(vertex,[upz,z])
        ham = jet1(graph)
        if ham > THRESHOLD*N:
            at_minimum = False
            print "Got out. Searching again for minimum..."
            while not at_minimum:
                at_minimum=find_well(graph)
            print "Found minimum again."
        jet2(graph)
    if dim>0:
        Print "Found miniwell."
    return graph
