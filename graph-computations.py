import time

def vect(z):
    """Returns the sphere vector associated to a complex number with north stereographical projection."""
    #if type(z) is np.complex128:
    uz = (abs(z)**2./4.-1.)/(abs(z)**2./4.+1.)
    ux = z.real/(abs(z)**2./4.+1.)
    uy = z.imag/(abs(z)**2./4.+1.)
    #else:
        #uz = 1
        #ux = 0
        #uy = 0
    return np.array([ux,uy,uz])

def z2u(graph):
    for vertex in graph.nodes():
        upz=graph.node[vertex]['upz']
        z=graph.node[vertex]['z']
        if upz:
            z = 4./z
        u=vect(z)
        graph.add_node(vertex,ux=u[0])
        graph.add_node(vertex,uy=u[1])
        graph.add_node(vertex,uz=u[2])

def distrib(graph):
    M=np.zeros((3,networkx.number_of_nodes(graph)))
    for (i,vertex) in enumerate(graph.nodes()):
        ux=graph.node[vertex]['ux']
        uy=graph.node[vertex]['uy']
        uz=graph.node[vertex]['uz']
        M[0][i]=ux
        M[1][i]=uy
        M[2][i]=uz
    (U,s,Vh)=scipy.linalg.svd(M)
    #print U
    return s
    
def under(z,w,powz,poww):
    return (1.+abs(z)**2./4.)**powz*(1.+abs(w)**2./4.)**poww

def jet1(graph, magn=0):
    """computes the hamiltonian and its gradient, for a given position (aka first three values attached do the edge)

    Its input is a graph with two labels on each vertex. The function transforms the graph to add a third label on each vertex, then returns the value of the hamiltonian.

    If more than two labels are assigned to a vertex beforehand, this function will erase them."""
    t = time.time()
    ham = 0.
    #for vertex in graph.nodes():
        #(upz,z) = graph.getVertex(vertex)[:2]
        #graph.setVertex(vertex,[upz,z,0])
    for edge in graph.edges():
        upz=graph.node[edge[0]]['upz']
        z=graph.node[edge[0]]['z']
        upw=graph.node[edge[1]]['upz']
        w=graph.node[edge[1]]['z']
        #(upz,z,dz) = graph.getVertex(edge[0])[:3]
        #(upw,w,dw) = graph.getVertex(edge[1])[:3]
        zbar=np.conjugate(z)
        wbar=np.conjugate(w)
        if upz == upw:
            ham += 1.5 -abs(z-w)**2./(2.*under(z,w,1,1))
            #dz += (wbar-zbar)*(1.+zbar*w/4.)/(2.*under(z,w,2,1))
            #dw += (zbar-wbar)*(1.+wbar*z/4.)/(2.*under(z,w,1,2))
        else:
            ham += 1.5 - abs(z*w/2.-2.)**2./(2.*under(z,w,1,1))
            #dz += -(zbar*wbar/4.-1.)*(zbar+w)/(2.*under(z,w,2,1))
            #dw += -(wbar*zbar/4.-1.)*(wbar+z)/(2.*under(z,w,1,2))
        if upz:
            ham -= magn*(1.-abs(z)**2./4.)*(1.+abs(z)**2./4.)/4.
        else:
            ham += magn*(1.-abs(z)**2./4.)*(1.+abs(z)**2./4.)/4.
        if upw:
            ham -= magn*(1.-abs(w)**2./4.)*(1.+abs(w)**2./4.)/4.
        else:
            ham += magn*(1.-abs(w)**2./4.)*(1.+abs(w)**2./4.)/4.
        #graph.add_node(edge[0],'dz'=dz)
        #graph.add_node(edge[1],'dz'=dw)
    #print "jet1 took ", time.time() - t
    return ham

def jetpart(graph, v0):
    t = time.time()
    ham = 0.
    modified = graph.neighbors(v0)
    modified.append(v0)
    for vertex in modified:
        upz=graph.node[vertex]['upz']
        z=graph.node[vertex]['z']
        zbar=np.conjugate(z)
        for nb in graph.neighbors(vertex):
            upw=graph.node[nb]['upz']
            w=graph.node[nb]['z']
            #(upw,w)= graph.getVertex(nb)[:2]
            wbar=np.conjugate(w)
            if upz==upw:
                ham += 1.5-abs(z-w)**2./(2.*under(z,w,1,1))
                dzdw = -(zbar-wbar)**2./(8.*under(z,w,1,1))
                dzdbarw = (1+zbar*w/4.)**2./(2.*under(z,w,1,1))
            else:
                ham += 1.5 - abs(z*w/2.-2.)**2./(2.*under(z,w,1,1))
                dzdw= -w**2.*(1.-zbar*wbar/4.)**2./(8.*under(z,w,1,1))
                dzdbarw=(wbar/4.)**2.*(w+zbar)**2./(2.*under(z,w,1,1))
            graph.add_edge(vertex,nb,dzdw=dzdw)
            graph.add_edge(vertex,nb,dzdbarw=dzdbarw)
            graph.add_edge(vertex,nb,origin=vertex)
    return [ham,mu]

def jet2(graph,magn=0):
    """computes the second jet, for a given position (aka first three values attached do the edge)

    Its input is a graph with three labels on each vertex. The function transforms the graph to add two labels on each vertex and two on each edge.

    If more than three labels are assigned to a vertex beforehand, this function will erase them."""
    t=time.time()
    for vertex in graph.nodes():
        graph.add_node(vertex,dzdz=0)
        graph.add_node(vertex,dzdbarz=0)
    for edge in graph.edges():
        upz=graph.node[edge[0]]['upz']
        z=graph.node[edge[0]]['z']
        dzdz=graph.node[edge[0]]['dzdz']
        dzdbarz=graph.node[edge[0]]['dzdbarz']
        upw=graph.node[edge[1]]['upz']
        w=graph.node[edge[1]]['z']
        dwdw=graph.node[edge[1]]['dzdz']
        dwdbarw=graph.node[edge[1]]['dzdbarz']
        zbar=np.conjugate(z)
        wbar=np.conjugate(w)
        if upz == upw:
            dzdz -= (wbar-zbar)*(1.+zbar*w/4.)*zbar/(4.*under(z,w,1,1))
            dzdbarz -= ((1.-abs(w)**2./4.)*(1.-abs(z)**2./4.)+(z*wbar).real)/(2.*under(z,w,1,1))
            dwdw -= (zbar-wbar)*(1.+wbar*z/4.)*wbar/(4.*under(z,w,1,1))
            dwdbarw -= ((1.-abs(w)**2./4.)*(1.-abs(z)**2./4.)+(z*wbar).real)/(2.*under(z,w,1,1))
            dzdw = -(zbar-wbar)**2./(8.*under(z,w,1,1)) #OK
            dzdbarw = (1+zbar*w/4.)**2./(2.*under(z,w,1,1)) #OK
        else:
            dzdz -= zbar*(1-zbar*wbar/4)*(w+zbar)/(4.*under(z,w,1,1))
            dzdbarz -= (-(1.-abs(w)**2/4.)*(1.-abs(z)**2/4.)+(z*w).real)/(2.*under(z,w,1,1))
            dwdw -= wbar*(1-wbar*zbar/4)*(z+wbar)/(4.*under(z,w,1,1))
            dwdbarw -= (-(1.-abs(z)**2/4.)*(1.-abs(w)**2/4.)+(w*z).real)/(2.*under(z,w,1,1))
            dzdw = (1.-zbar*wbar/4.)**2./(2.*under(z,w,1,1)) #OK
            dzdbarw = -(w+zbar)**2./(8.*under(z,w,1,1))
        if upz:
            dzdz+=magn*zbar*zbar/4.*(1.+abs(z)**2./4.)/4.
            dzdbarz+=-magn*(1.-abs(z)**2./4.)*(1.+abs(z)**2./4.)/4.
        else:
            dzdz+=-magn*zbar*zbar/4.*(1.+abs(z)**2./4.)/4.
            dzdbarz+=magn*(1.-abs(z)**2./4.)*(1.+abs(z)**2./4.)/4.
        if upw:
            dwdw+=magn*wbar*wbar/4.*(1.+abs(w)**2./4.)/4.
            dwdbarw+=-magn*(1.-abs(w)**2./4.)*(1.+abs(w)**2./4.)/4.
        else:
            dwdw+=-magn*wbar*wbar/4.*(1.+abs(w)**2./4.)/4.
            dwdbarw+=magn*(1.-abs(w)**2./4.)*(1.+abs(w)**2./4.)/4.
        graph.add_edge(edge[0],edge[1],dzdw=dzdw)
        graph.add_edge(edge[0],edge[1],dzdbarw=dzdbarw)
        graph.add_edge(edge[0],edge[1],origin=edge[0])
        graph.add_node(edge[0],dzdz=dzdz)
        graph.add_node(edge[0],dzdbarz=dzdbarz)
        graph.add_node(edge[1],dzdz=dwdw)
        graph.add_node(edge[1],dzdbarz=dwdbarw)

def hess(graph,stab=0.):
    t=time.time()
    """computes the Hessian with given 2-jet"""

    dbgraph=networkx.DiGraph()
    ham=jet1(graph)
    for vertex in graph.nodes():
        dzdz=graph.node[vertex]['dzdz']
        dzdbarz=graph.node[vertex]['dzdbarz']+stab
        #print dzdz
        #print dzdbarz
        #(dzdz,dzdbarz)=graph.getVertex(vertex)[3:]
        dbgraph.add_edge(str(vertex)+'x',str(vertex)+'x',val=-dzdz.imag)
        dbgraph.add_edge(str(vertex)+'y',str(vertex)+'y',val=dzdz.imag)
        dbgraph.add_edge(str(vertex)+'x',str(vertex)+'y',val=-dzdbarz-dzdz.real)
        dbgraph.add_edge(str(vertex)+'y',str(vertex)+'x',val=dzdbarz-dzdz.real)
    for e in graph.edges():
        dzdw=graph.edge[e[0]][e[1]]['dzdw']
        dzdbarw=graph.edge[e[0]][e[1]]['dzdbarw']
        if graph.edge[e[0]][e[1]]['origin']==e[1]:
            dzdbarw = np.conjugate(dzdbarw)
        dbgraph.add_edge(str(e[0])+'x',str(e[1])+'y',val=-(dzdbarw+dzdw).real)
        dbgraph.add_edge(str(e[0])+'x',str(e[1])+'x',val=(-dzdbarw+dzdw).imag)
        dbgraph.add_edge(str(e[0])+'y',str(e[1])+'y',val=-(dzdw+dzdbarw).imag)
        dbgraph.add_edge(str(e[0])+'y',str(e[1])+'x',val=(dzdbarw-dzdw).real)
        dbgraph.add_edge(str(e[1])+'x',str(e[0])+'y',val=-(dzdbarw+dzdw).real)
        dbgraph.add_edge(str(e[1])+'x',str(e[0])+'x',val=(dzdw+dzdbarw).imag)
        dbgraph.add_edge(str(e[1])+'y',str(e[0])+'y',val=(dzdbarw-dzdw).imag)
        dbgraph.add_edge(str(e[1])+'y',str(e[0])+'x',val=(dzdbarw-dzdw).real)
    
    hessian = networkx.adjacency_matrix(dbgraph,weight='val').toarray()
    return hessian

def realhess(graph):
    t=time.time()
    """computes the Hessian with given 2-jet"""
    dbgraph=networkx.DiGraph()
    for vertex in graph.nodes():
        dzdz=graph.node[vertex]['dzdz']
        dzdbarz=graph.node[vertex]['dzdbarz']
        #(dzdz,dzdbarz)=graph.getVertex(vertex)[3:]
        dbgraph.add_edge(str(vertex)+'x',str(vertex)+'x',val=dzdbarz+dzdz.real)
        dbgraph.add_edge(str(vertex)+'y',str(vertex)+'y',val=dzdbarz-dzdz.real)
        dbgraph.add_edge(str(vertex)+'x',str(vertex)+'y',val=-dzdz.imag)
        dbgraph.add_edge(str(vertex)+'y',str(vertex)+'x',val=-dzdz.imag)
    for e in graph.edges():
        dzdw=graph.edge[e[0]][e[1]]['dzdw']
        dzdbarw=graph.edge[e[0]][e[1]]['dzdbarw']
        if graph.edge[e[0]][e[1]]['origin']==e[1]:
            dzdbarw = np.conjugate(dzdbarw)
        dbgraph.add_edge(str(e[0])+'x',str(e[1])+'x',val=(dzdbarw+dzdw).real)
        dbgraph.add_edge(str(e[0])+'x',str(e[1])+'y',val=(-dzdbarw+dzdw).imag)
        dbgraph.add_edge(str(e[0])+'y',str(e[1])+'x',val=(dzdw+dzdbarw).imag)
        dbgraph.add_edge(str(e[0])+'y',str(e[1])+'y',val=(dzdbarw-dzdw).real)
        dbgraph.add_edge(str(e[1])+'x',str(e[0])+'x',val=(dzdbarw+dzdw).real)
        dbgraph.add_edge(str(e[1])+'x',str(e[0])+'y',val=(dzdw+dzdbarw).imag)
        dbgraph.add_edge(str(e[1])+'y',str(e[0])+'x',val=-(dzdbarw-dzdw).imag)
        dbgraph.add_edge(str(e[1])+'y',str(e[0])+'y',val=(dzdbarw-dzdw).real)
    
    hessian = networkx.adjacency_matrix(dbgraph,weight='val').toarray()
    return hessian


def charac(hessian):
    """Computes the characteristic, given the hessian"""
    #J=np.zeros(hessian.shape)
    #for i in range(int(0.5*hessian.shape[1])):
        #J[2*i,2*i+1]=-1
        #J[2*i+1,2*i]=1
    eigs=scipy.linalg.eigvals(hessian).imag #We want ALL the eigenvalues

    #There must be a better python way to do the following :
    mu=0
    for value in eigs:
        if value >0:
            mu += value
    return mu

def penalty(hessian):
    eigs=scipy.linalg.eigvals(hessian).real
    p=0
    for value in eigs:
        if value >0:
            p += value*value
    return p

#def charac(graph):
    #mu=0
    #for vertex in graph.nodes():
        #mu += graph.node[vertex]['dzdbarz']
    #print "mu=", mu
    #return mu

def initialize(graph):
    """Initializes the spins at random positions. Flushes the labels."""
    for vertex in graph.nodes():
        graph.add_node(vertex, upz=random.randint(0,1))
        graph.add_node(vertex, z=random.uniform(-2,2)+random.uniform(-2,2)*1j)
        #graph.setVertex(vertex, [random.randint(0,1),random.uniform(-2,2)+random.uniform(-2,2)*1j])
