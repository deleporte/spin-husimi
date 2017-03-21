import numpy as np
import scipy.linalg
import scipy.sparse
import random
import networkx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def gooddirections(graph):
    """
    From a graph configuration (at a minimum), gives a basis for the tangent space to the minimal set.
    """
    S=realhess(graph)
    values,vectors=scipy.linalg.eigh(S)
    threshold=4.*numpy.sqrt(hamil(graph))
    directions=[]
    for (i,val) in enumerate(values):
        if val < threshold:
            directions.append(vectors[i])
    return directions

def miniwell_sub(graph,magn=0):
    """
    Using a gradient descent, this function finds the minimal set for a given graph. Then, it uses a gradient descent on the tangent space of the minimal set.
    """
    THRESHOLD=0.000002
    DIFF_STEP=0.000001
    DESCENT=0.0005
    at_well=False
    at_miniwell=False
    nodes=sorted(graph.nodes(),key=str)
    N=len(nodes)
    while not at_miniwell:
        if not at_well:
            well(graph)
            at_well=True
            print "Minimal set found."
            ham=hamil(graph)
        directions = gooddirections(graph)
        print len(directions), "good directions."
        mu = charac(hess(graph))
        print mu
        #computing the gradient
        grad=[]
        gradnorm=0.
        for vector in directions:
            zmem=[]
            for (i,vertex) in enumerate(nodes):
                z=graph.node[vertex]['z']
                zmem.append(z)
                z+=DIFF_STEP*(1.+abs(z)**2./4.)*(vector[2*i]+vector[2*i+1]*1j)
                graph.add_node(vertex,z=z)
            jet2(graph)
            varmu=charac(hess(graph))
            grad.append((varmu-mu)/DIFF_STEP)
            gradnorm += ((varmu-mu)/DIFF_STEP)**2.
            for vertex in nodes[::-1]:
                z=zmem.pop()
                graph.add_node(vertex,z=z)
        if np.sqrt(gradnorm) < THRESHOLD*np.sqrt(N):
            at_miniwell=True
        #descent
        move = numpy.empty_like(directions[0])
        for (i,vector) in enumerate(directions):
            move += DESCENT*grad[i]*vector
        for (i,vertex) in enumerate(nodes):
            z=graph.node[vertex]['z']
            z -= (1.+abs(z)**2./4.)*(move[2*i] + move[2*i+1]*1j)
            graph.add_node(vertex,z=z)
            if abs(z)>2.1:
                graph.add_node(vertex,upz=1-graph.node[vertex]['upz'])
                graph.add_node(vertex,z=4./z)
        if hamil(graph) > 0.01:
            print "Got out !"
            at_well=False
        
        

def miniwell_gradient(graph,magn=0):
    """
    Input: an initialized graph
    Output: The graph at a miniwell

    Using a simple gradient descent, this function slowly finds a miniwell for a given graph. Now it works on two triangles.
    """
    if not graph.graph['init']:
        print "Graph not initialized !"
        return -1
    count=0
    THRESHOLD=0.000002
    DIFF_STEP=0.00001
    EPSILON=0.01
    DESCENT=0.04
    N=networkx.number_of_nodes(graph)
    #initialize(graph)
    at_miniwell= False
    ham = hamil(graph,magn)
    jet2(graph,magn)
    M=hess(graph)
    mu = charac(M)
    plt.ion()
    fig=plt.figure()
    ax=fig.add_subplot(111, projection='3d')
    points=ax.plot([],[],[],marker="o",color="b")[0]
    points.set_linewidth(0)
    zero=ax.scatter([0.],[0.],[0.],marker="o",color="r")
    lines=[]
    for k in range(len(graph.edges())):
        lines.append(ax.plot([0,0],[0,0],[0,0],color="g")[0])
    fig.canvas.draw()
    plt.show(block=False)
    image=0
    #while not at_miniwell:
    while count<1000:
        if not count%10:
            z2u(graph)
            s = distrib(graph)
            print "ham=",ham,"and mu=",mu
            #print "s=",s
        count += 1
        #computing the gradient
        grad=[]
        gradnorm=0
        ham = hamil(graph,magn)
        jet2(graph,magn)
        M=hess(graph)
        stab = 0.#max(0,-1.1*min(scipy.linalg.eigvalsh(M)))
        mu = charac(M)
        for vertex in graph.nodes():
            z=graph.node[vertex]['z']
            graph.add_node(vertex,z=z+ DIFF_STEP*(1.+abs(z)**2./4.))
            graph.graph['jetcalc']=False
            varham= hamil(graph,magn)
            jetpart(graph,vertex,magn)
            varmu = charac(hesspart(graph,vertex,M,stab))
            #varmu = charac(hess(graph,stab))
            diff=(varham-ham+EPSILON*(varmu-mu))/DIFF_STEP
            grad.append(diff)
            gradnorm+=diff**2
            graph.add_node(vertex,z=z+DIFF_STEP*(1.+abs(z)**2./4.)*1j)
            graph.graph['jetcalc']=False
            varham= hamil(graph,magn)
            jetpart(graph,vertex,magn)
            varmu = charac(hesspart(graph,vertex,M,stab))
            #varmu = charac(hess(graph,stab))
            diff=(varham-ham+EPSILON*(varmu-mu))/DIFF_STEP
            #print varmu
            grad.append(diff)
            graph.add_node(vertex,z=z)
            graph.graph['jetcalc']=False
            gradnorm+=diff**2
        if np.sqrt(gradnorm) < THRESHOLD*np.sqrt(N):
            at_miniwell=True
        #moving
        for (i,vertex) in enumerate(graph.nodes()):
            upz=graph.node[vertex]['upz']
            z=graph.node[vertex]['z']
            step=grad[2*i]+grad[2*i+1]*1j
            step*=DESCENT
            step*=(1.+abs(z)**2./4.)
            z-=step
            if abs(z) > 2.1:
                upz=1.-upz
                z=4./z
                graph.add_node(vertex,upz=upz)
            graph.add_node(vertex,z=z)
        graph.graph['jetcalc']=False
        #every 100th step we plot
        if not count%100:
            z2u(hus)
            Xpoints=[]
            Ypoints=[]
            Zpoints=[]
            for vertex in graph.nodes():
                Xpoints.append(graph.node[vertex]['ux'])
                Ypoints.append(graph.node[vertex]['uy'])
                Zpoints.append(graph.node[vertex]['uz'])
            points.remove()
            points=ax.scatter(Xpoints,Ypoints,Zpoints,marker="o",color="b")
            for (k,edge) in enumerate(graph.edges()):
                Xlines=[]
                Ylines=[]
                Zlines=[]
                Xlines.append(graph.node[edge[0]]['ux'])
                Xlines.append(graph.node[edge[1]]['ux'])
                Ylines.append(graph.node[edge[0]]['uy'])
                Ylines.append(graph.node[edge[1]]['uy'])
                Zlines.append(graph.node[edge[0]]['uz'])
                Zlines.append(graph.node[edge[1]]['uz'])
                lines[k].remove()
                lines[k],=ax.plot(Xlines,Ylines,Zlines)
                lines[k].set_color('g')
            fig.canvas.draw()
            count=0
            fig.savefig('record/image'+str(image).zfill(3)+'.png', dpi=fig.dpi)
            image += 1

def well(graph,magn=0):
    """
    Input: an initialized graph
    Output: The graph at a well

    Using a simple gradient descent, this function slowly finds a miniwell for a given graph. Now it works on two triangles.
    """
    if not graph.graph['init']:
        print "Graph not initialized !"
        return -1
    count=0
    THRESHOLD=0.000002
    DIFF_STEP=0.000001
    DESCENT=0.2
    N=networkx.number_of_nodes(graph)
    #initialize(graph)
    at_miniwell= False
    ham = hamil(graph,magn)
    while not at_miniwell:
        #computing the gradient
        grad=[]
        gradnorm=0
        ham = hamil(graph,magn)
        #print ham
        for vertex in graph.nodes():
            z=graph.node[vertex]['z']
            graph.add_node(vertex,z=z+ DIFF_STEP*(1.+abs(z)**2./4.))
            graph.graph['jetcalc']=False
            varham= hamil(graph,magn)
            diff=(varham-ham)/DIFF_STEP
            grad.append(diff)
            gradnorm+=diff**2
            graph.add_node(vertex,z=z+DIFF_STEP*(1.+abs(z)**2./4.)*1j)
            graph.graph['jetcalc']=False
            varham= hamil(graph,magn)
            diff=(varham-ham)/DIFF_STEP
            grad.append(diff)
            graph.add_node(vertex,z=z)
            graph.graph['jetcalc']=False
            gradnorm+=diff**2
        if np.sqrt(gradnorm) < THRESHOLD*np.sqrt(N):
            at_miniwell=True
        #moving
        for (i,vertex) in enumerate(graph.nodes()):
            upz=graph.node[vertex]['upz']
            z=graph.node[vertex]['z']
            step=grad[2*i]+grad[2*i+1]*1j
            step*=DESCENT
            step*=(1.+abs(z)**2./4.)
            z-=step
            if abs(z) > 2.1:
                upz=1.-upz
                z=4./z
                graph.add_node(vertex,upz=upz)
            graph.add_node(vertex,z=z)
        graph.graph['jetcalc']=False



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
