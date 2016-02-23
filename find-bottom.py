from numpy import *
from apgl.graph import *

# to each vertex of the graph are attached  values
# the first value is a boolean to see in which hemisphere our spin vector is
# the second value is a complex double denoting the stereographical projection
# the third is a complex double denoting the holomorphic derivative of the hamiltonian
# the last two are complex doubles for the hessian wrt variables on the same sphere (dzdz and dzdbarz)

# besides, to each edge of the graph is attached two complex doubles, denoting the various second derivatives involving two variables linked by an edge

def jet1(graph):
    """computes the hamiltonian and its gradient, for a given position (aka first three values attached do the edge)

    Its input is a graph with two labels on each vertex. The function transforms the graph to add a third label on each vertex, then returns the value of the hamiltonian.

    If more than two labels are assigned to a vertex beforehand, this function will erase them."""
    ham = 0
    for edge in graph.getAllEdge():
        upz=graph.getvertex(edge[0])[0]
        upw=graph.getvertex(edge[1])[0]
        z=graph.getvertex(edge[0])[1]
        w=graph.getvertex(edge[1])[1]
        if upz == upw:
            ham -= abs(z-w)**2./(2.*(abs(z)**2./4.+1.)*(abs(w)**2./4.+1.))
            dz = -(z-w).conjugate()*(1.+z.conjugate()*w/4.)/(2.*(abs(z)**2./4.+1.)**2.*(abs(w)**2./4.+1.))
            dw = -(w-z).conjugate()*(1.+w.conjugate()*z/4.)/(2.*(abs(w)**2./4.+1.)**2.*(abs(z)**2./4.+1.))
        else:
            ham -= abs(z*w/4.-1.)**2./(2.*(abs(z)**2./4.+1.)*(abs(w)**2./4.+1.))
            dz = -(z*w/4.-1.).conjugate()*(z.conjugate()+w)/(2.*(abs(z)**2./4.+1.)**2.*(abs(w)**2./4.+1.))
            dw = -(w*z/4.-1.).conjugate()*(w.conjugate()+z)/(2.*(abs(w)**2./4.+1.)**2.*(abs(z)**2./4.+1.))
        graph.setvertex(edge[0])[upz,z,dz]
        graph.setvertex(edge[1])[upw,w,dw]

    return ham
    
def hess(graph):
    """computes the hessian, for a given position (aka first three values attached do the edge)

    Its input is a graph with three labels on each vertex. The function transforms the graph to add two labels on each vertex and two on each edge.

    If more than three labels are assigned to a vertex beforehand, this function will erase them."""
    
    # To be done
    for edge in graph.getAllEdge():
        upz=graph.getvertex(edge[0])[0]
        upw=graph.getvertex(edge[1])[0]
        z=graph.getvertex(edge[0])[1]
        w=graph.getvertex(edge[1])[1]
        dz=graph.getvertex(edge[0])[2]
        dw=graph.getvertex(edge[1])[2] #we don't need them, except for stability at rewriting
        if upz=upw:
            dzdz=z.conjugate()/4.*(z-w).conjugate()*(1.+w*z.conjugate()/4.)/((abs(z)**2./4.+1.)**3.*(abs(w)**2./4.+1.))
            dzdbarz= #todo
            dzdw= #todo
            dzdbarw=#todo
        else:
            dzdz=#todo
            dzdbarz=#todo
            dzdw=#todo
            dzdbarw=#todo
        graph.setvertex(edge[0],[upz,z,dz,dzdz,dzdbarz])
        graph.setvertex(edge[1],[upw,w,dw, dwdw, dwdbarw])
        graph.addEdge(edge[0],edge[1],[dzdw,dzdbarw]) #the addEdge method erases the previous label (undocumented)
        
def charac(graph)
    """computes the Toeplitz characteristic value of the Hessian. Its entry is a graph with fully computed 2-jet (with jet1 and hess functions), it outputs a real double"""
    mu=0
    #todo
    return mu

def initialize(graph)
    """Initializes the spins at random positions. Flushes the labels."""


def findBottom(graph)
    """The findBottom function finds an
