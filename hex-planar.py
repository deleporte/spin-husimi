def inout():
    graph=buildHex()
    for vertex in graph.nodes():
        graph.add_node(vertex,upz=0)
    j=-1.+np.sqrt(3)*1j
    jj=-1.-np.sqrt(3)*1j
    graph.add_node(1,z=2)
    graph.add_node(23,z=2)
    graph.add_node(4,z=2)
    graph.add_node(56,z=2)
    graph.add_node(12,z=j)
    graph.add_node(3,z=j)
    graph.add_node(45,z=j)
    graph.add_node(6,z=j)
    graph.add_node(2,z=jj)
    graph.add_node(34,z=jj)
    graph.add_node(5,z=jj)
    graph.add_node(61,z=jj)
    return graph

def windmill():
    graph=buildHex()
    for vertex in graph.nodes():
        graph.add_node(vertex,upz=0)
    j=-1.+np.sqrt(3)*1j
    jj=-1.-np.sqrt(3)*1j
    graph.add_node(1,z=2)
    graph.add_node(2,z=j)
    graph.add_node(3,z=2)
    graph.add_node(4,z=j)
    graph.add_node(5,z=2)
    graph.add_node(6,z=j)
    graph.add_node(12,z=jj)
    graph.add_node(23,z=jj)
    graph.add_node(34,z=jj)
    graph.add_node(45,z=jj)
    graph.add_node(56,z=jj)
    graph.add_node(61,z=jj)
    return graph

def para():
    graph=buildHex()
    for vertex in graph.nodes():
        graph.add_node(vertex,upz=0)
    j=-1.+np.sqrt(3)*1j
    jj=-1.-np.sqrt(3)*1j
    graph.add_node(1,z=2)
    graph.add_node(23,z=2)
    graph.add_node(34,z=2)
    graph.add_node(5,z=2)
    graph.add_node(2,z=j)
    graph.add_node(4,z=j)
    graph.add_node(56,z=j)
    graph.add_node(61,z=j)
    graph.add_node(12,z=jj)
    graph.add_node(3,z=jj)
    graph.add_node(45,z=jj)
    graph.add_node(6,z=jj)
    return graph

def hex345():
    graph=buildHex()
    for vertex in graph.nodes():
        graph.add_node(vertex,upz=0)
    j=-1.+np.sqrt(3)*1j
    jj=-1.-np.sqrt(3)*1j
    graph.add_node(1,z=2)
    graph.add_node(2,z=j)
    graph.add_node(12,z=jj)
    graph.add_node(23,z=jj)
    graph.add_node(3,z=2)
    graph.add_node(34,z=j)
    graph.add_node(4,z=jj)
    graph.add_node(45,z=j)
    graph.add_node(5,z=2)
    graph.add_node(56,z=jj)
    graph.add_node(6,z=j)
    graph.add_node(61,z=jj)
    return graph

