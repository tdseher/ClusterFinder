######################################################################
#                                                                    #
#                 FIND GENE CLUSTERS GIVEN A THRESHOLD               #
#                           Peter Cimermancic                        #
#                               April 2010                           #
#                                                                    #
######################################################################

import numpy as np

def get_clusters(file,X=2,TRSH=0.2):
    '''
    open HMM output file and return lists of clusters
    containing the same atributes as in the HMM output:
     - domain with clusters probability higher than TRSH
     - join clusters that are less than X nt appart
     - join clusters that have breaks in the middle of the gene
    '''
    #2265131651      2265149182      3       224     3       224     Q       0.113366392809
    # read in HMM output file
    data = open(file,'r')
    D = data.readlines()
    data.close()
    
    # read in
    Genes = []
    visitedGenes = []
    Pfams = []
    NAME = '_'.join(D[0].strip().split('\t')[2].split(' '))

    for n,d in enumerate(D):
        d = d.strip().split('\t')
        
        seqID = d[0]
        geneFXN = d[1]
        geneID = d[2]
        gstart = int(d[3])
        gstop = int(d[4])
        pstart = int(d[5]) + gstart
        pstop = int(d[6]) + gstart
        pfamID = d[7]
        Prob = float(d[8])
    
        if geneID not in visitedGenes:
            visitedGenes.append(geneID)
        
        if Prob >= TRSH:
            Genes.append( geneID )
    
    for n,d in enumerate(D):
        d = d.strip().split('\t')
        geneID = d[2]
    
        if geneID in Genes:
            Pfams.append( (visitedGenes.index(geneID),d) )
    
    CLUSTERS = []
    cluster = []
    last = 0
    
    for x in xrange(len(Pfams)-1):            
        if (Pfams[x][1] == Pfams[x+1][0] or Pfams[x+1][0]-Pfams[x][0] <= X) and Pfams[x][1][0] == Pfams[x+1][1][0]:
            cluster.append(Pfams[x][1])
            last = x
        
        else:
            if last == x-1 and Pfams[x][1][0] == Pfams[x-1][1][0]:
                cluster.append(Pfams[x][1])
            if len(cluster) > 0:
                CLUSTERS.append(cluster)
                cluster = []
        
        if x+2 == len(Pfams):
            #cluster.append(Pfams[x][1])
            cluster.append(Pfams[x+1][1])
            if len(cluster) > 0:
                CLUSTERS.append(cluster)
    print "Clustering completed!"
    return CLUSTERS,NAME

'''
import glob

files = glob.glob('NEW_NEWEST_DRAFT/*')
#files = ['Streptomyces_coelicolor_A3(2).out','Sorangium_cellulosum_So_ce_56.out','Saccharopolyspora_erythraea_NRRL_2338.out']
N = 0
filtering = 1

for fil in files:
    T,Name = get_clusters(fil,X=2,TRSH=0.2)
    
    for n,t in enumerate(T):
        if filtering == 1:
            Size, Score = int(t[-1][7]) - int(t[0][6]), np.mean(np.array([float(t[i][16]) for i in xrange(len(t))]))
            if Size < 2000 or Score < 0.4:
                continue
        N += 1
        for m,g in enumerate(t):
            name = Name + '_%i' % n
            print '\t'.join([name]+[str(i) for i in g])
print N
'''

