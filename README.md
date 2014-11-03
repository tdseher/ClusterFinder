# ClusterFinder

### usage:
ClusterFinder.py [-h] [-g] [-m {HMM,viterbi}] input output organism

### description:
ClusterFinder - Predicting biosynthetic gene clusters in genomes

### input file format:
COLUMN DESCRIPTION:
 1. GeneID
 2. Sequencing status
 3. Organism name
 4. Scaffold OID
 5. Organism OID
 6. Locus Tag
 7. Gene Start
 8. Gene End
 9. Strand
 10. Pfam Template Start
 11. Pfam Template End
 12. Pfam Start
 13. Pfam End
 14. PfamID
 15. Pfam E-score
 16. Enzyme ID

If your input file format differs from the one above, please either modify the input file, or change the way this script parses the lines.

### output files:
	[organism].out
		tab-delimited file with same columns as [input], but with an
		additional column containing probability values
	[organism].clusters.out
		tab-delimited file with input and probabilities, but only for
		domains of gene clusters that have passed the filtering steps
	[method].[organism].png

### requirements:
 * Python (2.X)
 * numpy
 * matplotlib (optional)

### authors:
Original program by Peter Cimermancic & Michael Fischbach.
Modifications Copyright 2014 Thaddeus D. Seher ([@tdseher])

### references:
Cimermancic P, Medema MH, Claesen J, Kurita K, Wieland Brown LC, et al. (2014) Insights into secondary metabolism from a global analysis of prokaryotic biosynthetic gene clusters. Cell 158: 412-421. doi: 10.1016/j.cell.2014.06.034 (http://www.sciencedirect.com/science/article/pii/S0092867414008265)

### positional arguments:
	input                 input tab-delimited file
	output                name of directory to store output files
	organism              name of the organism

### optional arguments:
	-h, --help            show this help message and exit
	-g, --graph           generate a figure of the model (default: False)
	-m {HMM,viterbi}, --method {HMM,viterbi}
	                      Forward-Backward HMM or Viterbi alogirthm (default: HMM)

###example:
	$ python ClusterFinder.py example_input.txt example_org.out example_org.clusters.out

[@tdseher]:http://twitter.com/tdseher