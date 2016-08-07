import numpy as np

def distance((x1,y1),(x2,y2)):
    """Calculates distance between two coordinates

    Arguments:
    (x1,y1) -- First x,y coordinate as tuple
    (x2,y2) -- Second x,y coordinate as tuple

    Returns:
    distance

    """

    distance = np.sqrt((x1-x2)**2+(y1-y2)**2)
    return distance

def find_clusters(coords,threshold=0.03,indices=None):
    """Finds clusters in list of coordinates

    Keyword arguments:
    coords -- python list of tuples for peak coordinates. e.g. [(x1,y1),(x2,y2),(x3,y3)] 
    threshold -- distance within which peaks are considered clustered
    indices -- list of coordinate names or indices

    Returns:
    clusters -- list of clustered coordinates [[(index, (x, y)), (index, (x, y))],[etc]]

    """

    if indices is None:
        indices = range(len(coords))

    elif len(indices) != len(coords):
        raise ValueError("Length of indices must be the same as the number of coordinates")

    coords = zip(indices,coords)
    clusters = []
    cluster_indices = []
    while coords:
        clu = []
        sear = []

        t = coords.pop(-1)

        sear.append(t)
        clu.append(t)

        while sear:

            t2 = sear.pop(-1)

            for j in coords:
                if distance(t2[1],j[1])<threshold:
                    clu.append(j)
                    sear.append(j)

                    coords.remove(j)
#            print sear
        clusters.append(clu)
#    print clusters

    return clusters

def makeFudaList(clusters,paramfile="param.fuda"):
    """ Creates list of peaks for input into fuda

    Keyword arguments:
    clusters -- list of peak clusters [[(index, (x, y)), (index, (x, y))],[etc]]. Ideally the output of find_clusters. 
    paramfile -- name of the param file you want to append to (default "param.fuda")

    Returns:
    output file with OVERLAP peaks in FuDA format

    """

    out = open(paramfile,"a")

    for cluster in clusters:
        opeaks = "OVERLAP_PEAKS=("
        if len(cluster) == 1:
            pass
        else:
            for num, peak in enumerate(cluster):
                if num == len(cluster)-1:
                    opeaks+="%s)"%str(peak[0])
                else:
                    opeaks+="%s;"%str(peak[0])
            out.write(opeaks+"\n")
            #print "Peak %d with %s"%(peak+1,opeaks)

    out.close()

if __name__=="__main__":
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_style("ticks")

    test_data = "/Users/jacobbrady/Documents/NMRDATA/Ddx4/14ftoa/160723_250mgml/uoft800/T1/30degrees/15/fuda_analysis/peaks.fuda"
    peaks = pd.read_table(test_data,names=["Number","F1","F2"])
    """ Scale peaks by gamma """
    scaled_F1 = peaks.F1/10.
    coords = zip(scaled_F1,peaks.F2)#[:50]

    indices = list(peaks.Number)#[:50]
    clusters = find_clusters(coords,0.03,indices)
    print clusters

    makeFudaList(clusters)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for peak in clusters:
        #print list(peak)
        coords = [[i[0],i[1][0],i[1][1]] for i in peak]
        print coords
        coords = np.vstack(coords)
        inds = coords[:,0]
        x = coords[:,1]*10. # Rescale
        y = coords[:,2]
        ax.plot(y,x,"o",label="%s"%str(inds))

    ax.invert_yaxis()
    ax.invert_xaxis()
    plt.legend(loc=0,ncol=4)
    plt.show()
    






