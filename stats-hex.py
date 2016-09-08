exec open("build-husimi.py")
exec open("hex-planar.py")
exec open("graph-computations.py")
exec open("find-bottom.py")
import networkx
import numpy
import pylab as plt
import matplotlib.pyplot
from mpl_toolkits.mplot3d import Axes3D

def statshex(samples):
    hus=buildHex()
    with open('stats.csv','a') as f:
        for j in range(1,samples):
            initialize(hus)
            well(hus)
            jet1(hus)
            jet2(hus)
            mu=charac(hess(hus))
            if mu < 1.55:
                z2u(hus)
                s = distrib(hus)
                print(s)
            f.write(str(mu))
            f.write("\n")
            print j

def showstats(nbbins=50):
    with open('stats.csv','r') as f:
        X=numpy.loadtxt('stats.csv')
        plt.hist(X,bins=nbbins)
        plt.show()
            
def findLow():
    hus=buildHex()
    mumin=1.5
    while True:
        initialize(hus)
        well(hus)
        jet1(hus)
        jet2(hus)
        mu=charac(hess(hus))
        if mu < mumin:
            print "new minimum ! ", mu
            mumin=mu
            z2u(hus)
            s=distrib(hus)
            print(s)
            with open('lowkag.gml','w') as f:
                networkx.write_gml(hus,f,str)

def showLow():
    with open('lowkag.gml','r') as f:
        hus=networkx.read_gml(f)
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

