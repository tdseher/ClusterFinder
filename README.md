# ClusterFinder

## ClusterFinder.py

### usage:
ClusterFinder.py [-h] [-g] [-m {HMM,viterbi}] input output organism

### description:
ClusterFinder - Predicting biosynthetic gene clusters in genomes

ClusterFinder is a hidden Markov model-based probabilistic algorithm that aims to identify gene clusters of both known and unknown classes. ClusterFinder is based on a training set of 732 (minus 55) biosynthetic gene clusters (BGCs) plus some small molecule products. Also, Non-BGC regions were collected from 100 randomly selected genomes, defined as those regions without significant sequence similarity to the BGC state training set sequences.
  
To scan a genome for BGCs, first you must map your protein sequence to Pfam domains.
  
ClusterFinder assigns each domain a probability of being part of a gene cluster based on the frequencies at which these domains occur in the BGC and non-BGC training sets, and the identities of neighboring domains.
  
Since ClusterFinder is based solely on Pfam domain frequencies, and nature uses distinct assemblages of the same enzyme superfamilies to construct unrelated natural product classes, ClusterFinder exhibits relatively little training set bias and is capable of identifying new classes of gene clusters effectively.

The ClusterFinder prediction algorithm for BGC identification is a two-state Hidden Markov Model (HMM), with one hidden state corresponding to biosynthetic gene clusters (BGC state) and a second hidden state corresponding to the rest of the genome (non-BCG state).

### input file format:
COLUMN DESCRIPTION
 1. Protein sequence ID (treated as Locus Tag)
 2. Sequencing status (unused)
 3. Organism name (unused)
 4. Scaffold ID
 5. Organism ID (unused)
 6. Pfam domain name
 7. Gene start
 8. Gene end
 9. Strand (unused)
 10. Pfam template start (unused)
 11. Pfam template end (unused)
 12. Pfam start
 13. Pfam end
 14. Pfam ID
 15. Pfam E-score (unused)
 16. Enzyme ID (unused)

If your input file format differs from the one above, please either modify the input file, or change the way this script parses the lines.

### output file format:
COLUMN DESCRIPTION
 1. Cluster
 2. Scaffold ID
 3. Pfam domain name
 4. Protein sequence ID (Locus Tag)
 5. Gene Start
 6. Gene End
 7. Pfam Start
 8. Pfam End
 9. Pfam ID
 10. Probability

### output files:
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

### requirements:
 * [Python](http://www.python.org/)
 * [NumPy](http://www.numpy.org/)
 * [Matplotlib](http://www.matplotlib.org/) (optional)

### authors:
Original program by Peter Cimermancic & Michael Fischbach.
Modifications Copyright 2014 Thaddeus D. Seher ([@tdseher](http://www.twitter.com/tdseher)).

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
	$ cd example
	$ python ClusterFinder.py --graph example.input example.ClusterFinder Streptomyces_avermitilis_MA-4680

## ClusterFinder-prepare.py
### usage:
ClusterFinder-prepare.py [-h] [-s STATUS] [-o ORGANISM] [-c SCAFFOLD_ID] [-i ORGANISM_ID] \*.gff \*.pfam_scan

### description:
Combines a \*.gff (CDS annotations) file with a \*.pfam_scan (Pfam domains) file together to create the input necessary for ClusterFinder.

### requirements:
 * [Python](http://www.python.org/)
 * [HMMer 3.1](http://hmmer.janelia.org/)
 * [pfam_scan.pl](ftp://ftp.sanger.ac.uk/pub/databases/Pfam/Tools/PfamScan.tar.gz)
 * [Pfam-A.hmm](ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz)

### author:
Copyright 2014 Thaddeus D. Seher ([@tdseher](http://www.twitter.com/tdseher)).

### positional arguments:
	*.gff           Prodigal/RefSeq *.gff file
	*.pfam_scan     pfam_scan.pl *.pfam_scan file

### optional arguments:
	-h, --help      show this help message and exit
	-s STATUS, --status STATUS
	                sequencing satus (default: draft)
	-o ORGANISM, --organism ORGANISM
                    organism name (default: unknown)
	-c SCAFFOLD_ID, --scaffold_id SCAFFOLD_ID
	                scaffold id (default: unknown)
	-i ORGANISM_ID, --organism_id ORGANISM_ID
	                organism id (default: unknown)
	-e, --exclude_unmatched
	                do not include genes without Pfam domains (default: False)

### example:
	First, create a directory to store the Pfam databases
	$ mkdir Pfam-A+B.hmm
	$ cd Pfam-A+B.hmm
	
	Download the Pfam-A.hmm and Pfam-B.hmm databases for use with pfam_scan.pl
	$ wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
	$ wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
	$ wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-B.hmm.gz
	$ wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-B.hmm.dat.gz
	$ wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz
	
	Next, extract the database files
	$ gunzip *.gz
	
	Run hmmpress to create 4 additional files for each database
	$ hmmpress Pfam-A.hmm
	$ hmmpress Pfam-B.hmm
	
	Then obtain a sample genome or metagenome
	$ cd ..
	$ wget ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Streptomyces_avermitilis_MA_4680_uid57739/NC_003155.fna
	
	Run Prodigal on the genome FASTA and save its output as GFF and protein FASTA (include "-p meta" if scanning a metagenome)
	$ prodigal -i NC_003155.fna -a NC_003155.prodigal.faa -f gff > NC_003155.prodigal.gff 2> NC_003155.prodigal.err
	
	Run pfam_scan.pl on the protein output file (ClusterFinder does not use Pfam-B.hmm domains)
	$ pfam_scan.pl -outfile NC_003155.pfam_scan -fasta NC_003155.prodigal.faa -dir Pfam-A+B.hmm
	
	Use this program to convert to input required for ClusterFinder
	$ python ClusterFinder-prepare.py NC_003155.prodigal.gff NC_003155.pfam_scan --status finished --organism 'Streptomyces avermitilis MA-4680' --scaffold_id 'gi|148878541|dbj|BA000030.3|' --organism_id 227882 > NC_003155.prepare
	
	Run ClusterFinder
	$ python ClusterFinder.py NC_003155.prepare NC_003155.ClusterFinder Streptomyces_avermitilis_MA-4680
