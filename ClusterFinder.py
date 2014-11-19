#!/usr/bin/env python

"""ClusterFinder - main script"""

# define imports
import sys
import os
PREFIX_PATH = os.path.dirname(os.path.realpath(sys.argv[0]))
PATH = os.path.join(PREFIX_PATH, "scr")

sys.path.append(PATH)
from Predict import *
from FindClusters import *
import re
import argparse
from operator import itemgetter

# Define global variables
__date__ = "2014"
__author__ = "Thaddeus D. Seher"
__twitter__ = "@tdseher"
__program__ = os.path.basename(sys.argv[0])
__citation__ = """Cimermancic P, Medema MH, Claesen J, Kurita K, Wieland
  Brown LC, et al. (2014) Insights into secondary metabolism
  from a global analysis of prokaryotic biosynthetic gene
  clusters. Cell 158: 412-421. doi: 10.1016/j.cell.2014.06.034"""
__description__ = """\
description:
  ClusterFinder - Predicting biosynthetic gene clusters in genomes

  ClusterFinder is a hidden Markov model-based probabilistic algorithm
  that aims to identify gene clusters of both known and unknown classes.
  ClusterFinder is based on a training set of 732 (minus 55) biosynthetic
  gene clusters (BGCs) plus some small molecule products. Also, Non-BGC
  regions were collected from 100 randomly selected genomes, defined as
  those regions without significant sequence similarity to the BGC state
  training set sequences.
  
  To scan a genome for BGCs, first you must map your protein sequence
  to Pfam domains.
  
  ClusterFinder assigns each domain a probability of being part of a gene
  cluster based on the frequencies at which these domains occur in the
  BGC and non-BGC training sets, and the identities of neighboring domains.
  
  Since ClusterFinder is based solely on Pfam domain frequencies, and
  nature uses distinct assemblages of the same enzyme superfamilies to
  construct unrelated natural product classes, ClusterFinder exhibits
  relatively little training set bias and is capable of identifying new
  classes of gene clusters effectively.
  
  The ClusterFinder prediction algorithm for BGC identification is a
  two-state Hidden Markov Model (HMM), with one hidden state
  corresponding to biosynthetic gene clusters (BGC state) and a second
  hidden state corresponding to the rest of the genome (non-BCG state).

input file format:
  COLUMN DESCRIPTION
       1 Protein sequence ID (treated as Locus Tag)
       2 Sequencing status (unused)
       3 Organism name (unused)
       4 Scaffold ID
       5 Organism ID (unused)
       6 Pfam domain name
       7 Gene start
       8 Gene end
       9 Strand (unused)
      10 Pfam template start (unused)
      11 Pfam template end (unused)
      12 Pfam start
      13 Pfam end
      14 PfamID
      15 Pfam E-score (unused)
      16 Enzyme ID (unused)
  
  If your input file format differs from the one above, please
  either modify the input file, or change the way this script
  parses the lines.

output file format:
  COLUMN DESCRIPTION
       1 Cluster
       2 Scaffold ID
       3 Pfam domain name
       4 Protein sequence ID (Locus Tag)
       5 Gene Start
       6 Gene End
       7 Pfam Start
       8 Pfam End
       9 PfamID
      10 Probability

output files:
  [output]/[organism].out
    tab-delimited file with select columns from [input], but
    with an additional column containing probability values
  [output]/[organism].clusters.out
    tab-delimited file with select columns from [input] and
    probabilities, but with an additional column specifying
    gene clusters for domains that have passed the filtering
    steps
  [output]/[method].[organism].png
    line plot generated when "--graph" option is enabled

requirements:
 * Python (2.X)
 * numpy
 * matplotlib (optional)
  
authors:
  Original program by Peter Cimermancic & Michael Fischbach.
  Modifications Copyright {__date__} {__author__} ({__twitter__}).

reference:
  {__citation__}

""".format(**locals())

__epilog__ = """\
example:
  $ python {__program__} example_input.txt example_org.out example_org.clusters.out

""".format(**locals())

def parse_arguments():
    """Creates the argument parser"""
    
    # --- Create the argument parser
    parser = argparse.ArgumentParser(
        description=__description__,
        epilog=__epilog__,
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # --- Add mandatory arguments
    parser.add_argument("input", help="input tab-delimited file")
    parser.add_argument("output", help="name of directory to store output files")
    parser.add_argument("organism", help="name of the organism")
    
    # --- Add optional arguments
    parser.add_argument('-g', '--graph', action='store_true', help='generate a figure of the model (default: False)')
    parser.add_argument('-m', '--method', type=str, choices=['HMM', 'viterbi'], default='HMM', help='Forward-Backward HMM or Viterbi alogirthm (default: HMM)')
    
    # --- parse the arguments
    args = parser.parse_args()
    
    return args

def main():
    # --- load arguments and parse them
    args = parse_arguments()
    
    
    # --- read in the HMM model data
    PATH_TO_FREQUENCIES = os.path.join(PREFIX_PATH, "freqdata") # link to where the frequencies are
    
    # Creates the output directory if it does not exist
    try: 
        os.makedirs(args.output)
    except OSError:
        if not os.path.isdir(args.output):
            raise
    
    out4 = open(os.path.join(PATH_TO_FREQUENCIES, 'TP_arr_A_latest.pkl'), 'rb')
    out5 = open(os.path.join(PATH_TO_FREQUENCIES, 'NewTS_all_B_reduced_6filter.pkl'), 'rb')
    out6 = open(os.path.join(PATH_TO_FREQUENCIES, 'SP_arr.pkl'), 'rb')
    out7 = open(os.path.join(PATH_TO_FREQUENCIES, 'NewTS_all_B_index.pkl'), 'rb')
    
    A=pickle.load(out4)
    B=pickle.load(out5)
    start_p=pickle.load(out6)
    index_dict=pickle.load(out7)
    
    out4.close()
    out5.close()
    out6.close()
    out7.close()
    
    
    # --- read in an input file
    data = open(args.input)
    D = data.readlines()
    data.close()
    
    DATA = []
    for d in D:
        d = d.strip().split('\t')
        
        ########## input file format parsing ##########
        ###### customize according to your input ######
        gene_id = d[0]
        gene,pfamID = d[5],d[13].replace('pfam','PF')
        genstart,genestop = int(d[6]),int(d[7])
        try: pfamstart,pfamstop = int(d[11]),int(d[12])
        except ValueError: pfamstart,pfamstop = genstart,genestop
        scaffold = d[3]
        ###############################################
        
        # - if gene has a Pfam annotation
        if pfamID != 'n/a':
            DATA.append([scaffold, gene, gene_id, genstart, genestop, pfamstart, pfamstop, pfamID])
        # - else
        else:
            DATA.append([scaffold, gene, gene_id, genstart, genestop, pfamstart, pfamstop, 'Q'])
    
    
    # --- sort Pfams as they appear in the genome
    PfamSorted = sorted(DATA, key=itemgetter(0,3,5))
    
    
    # --- filter out the repeats (optional)
    repeats = set([
        'PF07721', 'PF05593', 'PF07719', 'PF00515', 'PF00132', 'PF03130',
        'PF01839', 'PF01816', 'PF07720', 'PF00400', 'PF05594', 'PF07661',
        'PF02985', 'PF06049', 'PF08238', 'PF06696', 'PF00353', 'PF02412',
        'PF00023', 'PF02071', 'PF03991', 'PF01469', 'PF07676', 'PF00514',
        'PF00904', 'PF07634', 'PF02370', 'PF03335', 'PF01851', 'PF04728',
        'PF06715', 'PF03373', 'PF04680', 'PF00805', 'PF04508', 'PF07918',
        'PF01535', 'PF01011', 'PF05017', 'PF06671', 'PF00818', 'PF03406',
        'PF00399', 'PF09373', 'PF01744', 'PF01436', 'PF01239', 'PF05906',
        'PF03729', 'PF00404', 'PF04022', 'PF02363', 'PF02524', 'PF07981',
        'PF02095', 'PF00414', 'PF00560', 'PF05001', 'PF02162', 'PF01473',
        'PF05465', 'PF02493', 'PF03578', 'PF08043', 'PF06392', 'PF07142',
        'PF08309', 'PF02184'
    ])
    PfamSorted = [i for i in PfamSorted if i[-1] not in repeats]
    
    
    # - check that there are Pfam domain in the data
    pfpatt = re.compile(r"""PF+\d+\d+\d+\d+\d""")
    if len([re.match(pfpatt, p[-1]) for p in PfamSorted])==0:
        print 'There is no Pfam domain in you dataset!'
        exit()
    
    
    # --- run the HMM
    do_HMM([i[-1] for i in PfamSorted], A, B, index_dict, start_p, graph=args.graph, report=True, outs=PfamSorted, method=args.method, path=args.output, name=args.organism)
    
    
    # --- filter for biosynthetic domains
    data = open(os.path.join(PATH_TO_FREQUENCIES, 'biosynthetic_pfams.txt'), 'r')
    BioPfams = set([i.strip().split('\t')[0] for i in data.readlines()])
    data.close()
    
    
    # --- run clustering method
    Clusters, OrgName = get_clusters(os.path.join(args.output, args.organism+'.out'), X=2, TRSH=0.2)
    CLUSTERS = [] # store the new clusters
    for clstr in Clusters:
        pfams =  set([i[-4] for i in clstr])
        if len(pfams - BioPfams)>=3: # - requires at least three biosynthetic domains
            CLUSTERS.append(clstr)
    
    
    # --- write the data (in ClusterFinder format)
    output = open(os.path.join(args.output, args.organism+'.clusters.out'), 'w')
    for x,clstr in enumerate(CLUSTERS):
        for domain in clstr:
            output.write('\t'.join([str(j) for j in [args.organism+'_'+str(x)]+domain])+'\n')
    output.close()


if __name__ == '__main__':
    main()

