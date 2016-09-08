import numpy as np
import scipy.linalg
import scipy.sparse
import random

def miniwell(graph,magn=0):
    """
    Input: an initialized graph
    Output: The graph at a miniwell

    Using a simple gradient descent, this function slowly finds a miniwell for a given graph. Now it works on two triangles.
    """
    count=0
    THRESHOLD=0.000002
    DIFF_STEP=0.000001
    EPSILON=0.001
    DESCENT=0.3
    N=networkx.number_of_nodes(graph)
    #initialize(graph)
    at_miniwell= False
    ham = jet1(graph,magn)
    jet2(graph,magn)
    M=hess(graph)
    mu = charac(M)
    while not at_miniwell:
        if not count%10:
            z2u(graph)
            s = distrib(graph)
            print "ham=",ham,"and mu=",mu
            #print "s=",s
            count =0
        count += 1
        #computing the gradient
        grad=[]
        gradnorm=0
        ham = jet1(graph,magn)
        jet2(graph,magn)
        M=realhess(graph)
        stab = max(0,-1.1*min(scipy.linalg.eigvalsh(M)))
        mu = charac(hess(graph,stab))
        for vertex in graph.nodes():
            z=graph.node[vertex]['z']
            graph.add_node(vertex,z=z+ DIFF_STEP)
            varham= jet1(graph,magn)
            jet2(graph,magn)
            varmu = charac(hess(graph,stab))
            diff=(varham-ham+EPSILON*(varmu-mu))/DIFF_STEP
            grad.append(diff)
            gradnorm+=diff**2
            graph.add_node(vertex,z=z+DIFF_STEP*1j)
            varham= jet1(graph,magn)
            jet2(graph,magn)
            varmu = charac(hess(graph,stab))
            diff=(varham-ham+EPSILON*(varmu-mu))/DIFF_STEP
            grad.append(diff)
            graph.add_node(vertex,z=z)
            gradnorm+=diff**2
        if np.sqrt(gradnorm) < THRESHOLD*np.sqrt(N):
            at_miniwell=True
        #moving
        for (i,vertex) in enumerate(graph.nodes()):
            upz=graph.node[vertex]['upz']
            z=graph.node[vertex]['z']
            step=grad[2*i]+grad[2*i+1]*1j
            step*=DESCENT
            z-=step
            if abs(z) > 4:
                upz=1.-upz
                z=4./z
                graph.add_node(vertex,upz=upz)
            graph.add_node(vertex,z=z)

def well(graph):
    """
    Input: an initialized graph
    Output: The graph at a well

    Using a simple gradient descent, this function slowly finds a well for a given graph.
    """
    count=0
    THRESHOLD=0.000002
    DIFF_STEP=0.000001
    EPSILON=0#.001
    DESCENT=1.2
    N=networkx.number_of_nodes(graph)
    #initialize(graph)
    at_miniwell= False
    ham = jet1(graph)
    #jet2(graph)
    #mu = charac(hess(graph))
    mu=0
    while not at_miniwell:
        if not count%10:
            z2u(graph)
            s = distrib(graph)
            #print "ham=",ham#,"and mu=",mu
            #print "s=",s
            count =0
        count += 1
        #computing the gradient
        grad=[]
        gradnorm=0
        ham = jet1(graph)
        mu=0
        varmu=0
        #jet2(graph)
        #mu = charac(hess(graph))
        for vertex in graph.nodes():
            z=graph.node[vertex]['z']
            #[upz,z] = graph.getVertex(vertex)[:2]
            graph.add_node(vertex,z=z+ DIFF_STEP)
            #graph.setVertex(vertex,[upz,z+ DIFF_STEP])
            varham= jet1(graph)
            #jet2(graph)
            #varmu = charac(hess(graph))
            diff=(varham-ham+EPSILON*(varmu-mu))/DIFF_STEP
            grad.append(diff)
            gradnorm+=diff**2
            graph.add_node(vertex,z=z+DIFF_STEP*1j)
            #graph.setVertex(vertex,[upz,z+ DIFF_STEP*1j])
            varham= jet1(graph)
            #jet2(graph)
            #varmu = charac(hess(graph))
            diff=(varham-ham+EPSILON*(varmu-mu))/DIFF_STEP
            grad.append(diff)
            graph.add_node(vertex,z=z)
            #graph.setVertex(vertex,[upz,z])
            gradnorm+=diff**2
        if np.sqrt(gradnorm) < THRESHOLD*np.sqrt(N):
            at_miniwell=True
        #moving
        for (i,vertex) in enumerate(graph.nodes()):
            upz=graph.node[vertex]['upz']
            z=graph.node[vertex]['z']
            #[upz,z] = graph.getVertex(vertex)[:2]
            step=grad[2*i]+grad[2*i+1]*1j
            step*=DESCENT
            z-=step
            if abs(z) > 2.1:
                upz=1.-upz
                z=4./z
                graph.add_node(vertex,upz=upz)
            graph.add_node(vertex,z=z)
            #graph.setVertex(vertex,[upz,z])    

def show(hus):
    z2u(hus)
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111, projection='3d')
    X=[]
    Y=[]
    Z=[]
    for vertex in hus.nodes():
        X.append(hus.node[vertex]['ux'])
        Y.append(hus.node[vertex]['uy'])
        Z.append(hus.node[vertex]['uz'])
    ax.scatter(X,Y,Z,marker='o',color='b')
    X=[0]
    Y=[0]
    Z=[0]
    ax.scatter(X,Y,Z,marker='o',color='r')
    for edge in hus.edges():
        X=[]
        Y=[]
        Z=[]
        X.append(hus.node[edge[0]]['ux'])
        X.append(hus.node[edge[1]]['ux'])
        Y.append(hus.node[edge[0]]['uy'])
        Y.append(hus.node[edge[1]]['uy'])
        Z.append(hus.node[edge[0]]['uz'])
        Z.append(hus.node[edge[1]]['uz'])
        ax.plot(X,Y,Z, color='g')
    matplotlib.pyplot.show()
