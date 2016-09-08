def buildChain(n=2):
    if n<= 1:
        graph=networkx.Graph()
        graph.add_edge(0,1)
        graph.add_edge(0,2)
        graph.add_edge(1,2)
        return graph
    else:
        graph=buildChain(n-1)
        graph.add_edge(2*n-2,2*n-1)
        graph.add_edge(2*n-2,2*n)
        graph.add_edge(2*n-1,2*n)
        return graph

def rot(u,v,theta):
    au=np.array(u)
    av=np.array(v)
    n=np.sqrt(np.dot(np.cross(au,av),np.cross(au,av)))
    aw=1./n*np.cross(au,av)
    ax=-np.cross(au,aw)
    return list(-0.5*au+np.sqrt(3.)/2.*np.cos(theta)*ax+np.sqrt(3.)/2.*np.sin(theta)*aw)
    

def u2z(graph):
    for vertex in graph.nodes():
        ux=graph.node[vertex]['ux']
        uy=graph.node[vertex]['uy']
        uz=graph.node[vertex]['uz']
        re=ux*2./(1.-uz)
        im=uy*2./(1.-uz)
        graph.add_node(vertex,z=re+im*1j)
        graph.add_node(vertex,upz=0)
        

def initChain(T=[]):
    graph=buildChain(len(T)+1)
    graph.add_node(0,ux=1.)
    graph.add_node(0,uy=0.)
    graph.add_node(0,uz=0.)
    graph.add_node(1,ux=-0.5)
    graph.add_node(1,uy=-np.sqrt(3)/2.)
    graph.add_node(1,uz=0.)
    graph.add_node(2,ux=-0.5)
    graph.add_node(2,uy=np.sqrt(3)/2.)
    graph.add_node(2,uz=0.)
    for k,theta in enumerate(T):
        u=[]
        v=[]
        u.append(graph.node[2*k+2]['ux'])
        u.append(graph.node[2*k+2]['uy'])
        u.append(graph.node[2*k+2]['uz'])
        v.append(graph.node[2*k+1]['ux'])
        v.append(graph.node[2*k+1]['uy'])
        v.append(graph.node[2*k+1]['uz'])
        w=rot(u,v,theta)
        graph.add_node(2*k+3,ux=w[0])
        graph.add_node(2*k+3,uy=w[1])
        graph.add_node(2*k+3,uz=w[2])
        v=[]
        v.append(graph.node[2*k]['ux'])
        v.append(graph.node[2*k]['uy'])
        v.append(graph.node[2*k]['uz'])
        w=rot(u,v,theta)
        graph.add_node(2*k+4,ux=w[0])
        graph.add_node(2*k+4,uy=w[1])
        graph.add_node(2*k+4,uz=w[2])
    u2z(graph)
    return graph
