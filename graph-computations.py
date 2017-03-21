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
    """From stereographic representation to vectors on a sphere"""
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
    """Studies the space distribution of the spins.

    Returns the length of the three principal axes.
    If the last one is small, then the spins are coplanar.
    """
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

def hamil(graph, magn=0):
    """computes the hamiltonian, for a given position
    """
    if not graph.graph['init']:
        print "Graph not initialized !"
        return Nan
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

def jet2(graph,magn=0):
    """computes the second jet, for a given position

    """
    if not graph.graph['init']:
        print "Graph not initialized !"
        return Nan
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
    graph.graph['jetcalc']=True
    
def jetpart(graph,vertex,magn=0):
    """Updates the 2-jet when one spin is changed
    """
    effects=[vertex]
    for nb in graph.neighbors(vertex):
        effects.append(nb)
    for v in effects:
        graph.add_node(v,dzdz=0)
        graph.add_node(v,dzdbarz=0)
    fulledges = graph.subgraph(effects).edges()
    halfedges = []
    for e in graph.edges():
        if (e[0] in effects and e[1] not in effects):
            halfedges.append(e)
        if (e[1] in effects and e[0] not in effects):
            halfedges.append(e[1::-1])
    alledges = [e for e in fulledges]
    alledges.extend(halfedges)
    for edge in alledges:
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
        if edge in fulledges:
            graph.add_node(edge[1],dzdz=dwdw)
            graph.add_node(edge[1],dzdbarz=dwdbarw)
    graph.graph['jetcalc']=True
    
def hess(graph,stab=0.):
    """Computes J*the hessian"""
    if not graph.graph['init']:
        print "Graph not initialized !"
        return Nan
    if not graph.graph['jetcalc']:
        jet2(graph)
    t=time.time()
    """computes the Hessian with given 2-jet"""
    dbgraph=networkx.DiGraph()
    for vertex in graph.nodes():
        dzdz=graph.node[vertex]['dzdz']
        dzdbarz=graph.node[vertex]['dzdbarz']+stab
        dbgraph.add_edge(str(vertex)+'!x',str(vertex)+'!x',val=-dzdz.imag)
        dbgraph.add_edge(str(vertex)+'!y',str(vertex)+'!y',val=dzdz.imag)
        dbgraph.add_edge(str(vertex)+'!x',str(vertex)+'!y',val=-dzdbarz-dzdz.real)
        dbgraph.add_edge(str(vertex)+'!y',str(vertex)+'!x',val=dzdbarz-dzdz.real)
    for e in graph.edges():
        dzdw=graph.edge[e[0]][e[1]]['dzdw']
        dzdbarw=graph.edge[e[0]][e[1]]['dzdbarw']
        if graph.edge[e[0]][e[1]]['origin']==e[1]:
            dzdbarw = np.conjugate(dzdbarw)
        dbgraph.add_edge(str(e[0])+'!x',str(e[1])+'!y',val=-(dzdbarw+dzdw).real)
        dbgraph.add_edge(str(e[0])+'!x',str(e[1])+'!x',val=(-dzdbarw+dzdw).imag)
        dbgraph.add_edge(str(e[0])+'!y',str(e[1])+'!y',val=-(dzdw+dzdbarw).imag)
        dbgraph.add_edge(str(e[0])+'!y',str(e[1])+'!x',val=(dzdbarw-dzdw).real)
        dbgraph.add_edge(str(e[1])+'!x',str(e[0])+'!y',val=-(dzdbarw+dzdw).real)
        dbgraph.add_edge(str(e[1])+'!x',str(e[0])+'!x',val=(dzdw+dzdbarw).imag)
        dbgraph.add_edge(str(e[1])+'!y',str(e[0])+'!y',val=(dzdbarw-dzdw).imag)
        dbgraph.add_edge(str(e[1])+'!y',str(e[0])+'!x',val=(dzdbarw-dzdw).real)
        
    hessian = networkx.adjacency_matrix(dbgraph,nodelist=sorted(dbgraph.nodes()),weight='val').toarray()
    return hessian

def hesspart(graph,vertex,hessprev,stab=0.):
    """Updates the symplectic hessian when one spin is changed"""
    hess=np.empty_like(hessprev)
    hess[:] = hessprev
    effects=[vertex]
    for nb in graph.neighbors(vertex):
        effects.append(nb)
    vertices = sorted(graph.nodes(),key=str)
    for v in effects:
        i = vertices.index(v)
        dzdz=graph.node[v]['dzdz']
        dzdbarz=graph.node[v]['dzdbarz']+stab
        hess[2*i,2*i]=-dzdz.imag
        hess[2*i+1,2*i+1]=dzdz.imag
        hess[2*i,2*i+1]=-dzdbarz-dzdz.real
        hess[2*i+1,2*i]=dzdbarz-dzdz.real
    for e in graph.subgraph(effects).edges():
        i = vertices.index(e[0])
        j = vertices.index(e[1])
        dzdw=graph.edge[e[0]][e[1]]['dzdw']
        dzdbarw=graph.edge[e[0]][e[1]]['dzdbarw']
        if graph.edge[e[0]][e[1]]['origin']==e[1]:
            dzdbarw = np.conjugate(dzdbarw)
        hess[2*i,2*j+1]=-(dzdbarw+dzdw).real
        hess[2*i,2*j]=(-dzdbarw+dzdw).imag    
        hess[2*i+1,2*j+1]=-(dzdw+dzdbarw).imag
        hess[2*i+1,2*j]=(dzdbarw-dzdw).real
        hess[2*j,2*i+1]=-(dzdbarw+dzdw).real
        hess[2*j,2*i]=(dzdw+dzdbarw).imag
        hess[2*j+1,2*i+1]=(dzdbarw-dzdw).imag
        hess[2*j+1,2*i]=(dzdbarw-dzdw).real
    return hess

def realhess(graph):
    """Computes the standard hessian (definite symmetric)"""
    if not graph.graph['jetcalc']:
        jet2(graph)
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
    
    hessian = networkx.adjacency_matrix(dbgraph,nodelist=sorted(dbgraph.nodes()),weight='val').toarray()
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
    """Computes how much the hessian is not positive"""
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
    """Initializes the spins at random positions."""
    for vertex in graph.nodes():
        graph.add_node(vertex, upz=random.randint(0,1))
        graph.add_node(vertex, z=random.uniform(-2,2)+random.uniform(-2,2)*1j)
        #graph.setVertex(vertex, [random.randint(0,1),random.uniform(-2,2)+random.uniform(-2,2)*1j])
    graph.graph['init']=True
