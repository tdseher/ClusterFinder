#!/usr/bin/env python
import sys
import os
import re
import argparse

# Define global variables
__date__ = "2014"
__author__ = "Thaddeus D. Seher"
__twitter__ = "@tdseher"
__program__ = os.path.basename(sys.argv[0])
__description__ = """\
description:
  Combines a *.gff (CDS annotations) file with a *.pfam_scan (Pfam domains)
  file together to create the input necessary for ClusterFinder.

requirements:
  Python (http://www.python.org/)
  HMMer 3.1 (http://hmmer.janelia.org/)
  pfam_scan.pl (ftp://ftp.sanger.ac.uk/pub/databases/Pfam/Tools/PfamScan.tar.gz)
  Pfam-A.hmm (ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz)

author:
  Copyright {__date__} {__author__} ({__twitter__}).

""".format(**locals())

__epilog__ = """\
example:
  First, create a directory to store the Pfam databases
  $ mkdir Pfam-A+B.hmm
  $ cd Pfam-A+B.hmm
  
  Download the Pfam-A.hmm and Pfam-B.hmm databases for use with pfam_scan.pl
  $ wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
  $ wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
  $ wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-B.hmm.gz
  $ wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-B.hmm.dat.gz
  $ ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz
  
  Next, extract the database files
  $ gunzip *.gz
  
  Run hmmpress to create 4 additional files for each database
  $ hmmpress Pfam-A.hmm
  $ hmmpress Pfam-B.hmm
  
  Then obtain a sample genome or metagenome
  $ cd ..
  $ wget ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Streptomyces_avermitilis_MA_4680_uid57739/NC_003155.fna
  
  Run Prodigal on the genome FASTA and save its output as GFF and
  protein FASTA (include "-p meta" if scanning a metagenome)
  $ prodigal -i NC_003155.fna -a NC_003155.prodigal.faa -f gff
    > NC_003155.prodigal.gff 2> NC_003155.prodigal.err
  
  Run pfam_scan.pl on the protein output file (ClusterFinder does
  not use Pfam-B.hmm domains)
  $ pfam_scan.pl -outfile NC_003155.pfam_scan -fasta NC_003155.prodigal.faa
    -dir Pfam-A+B.hmm
  
  Use this program to convert to input required for ClusterFinder
  $ python {__program__} NC_003155.prodigal.gff NC_003155.pfam_scan
    --status finished --organism 'Streptomyces avermitilis MA-4680'
    --scaffold_id 'gi|148878541|dbj|BA000030.3|' --organism_id 227882
    > NC_003155.prepare
  
  Run ClusterFinder
  $ python ClusterFinder.py NC_003155.prepare NC_003155.ClusterFinder Streptomyces_avermitilis_MA-4680

""".format(**locals())

# target output (see README.md)
# GeneID	Sequencing status	Organism name	Scaffold ID	Organism ID	Locus Tag	Gene Start	Gene End	Strand	Pfam Template Start	Pfam Template End	Pfam Start	Pfam End	PfamID	Pfam E-score	Enzyme ID
# 637203094	Finished	Streptomyces avermitilis MA-4680	637000121	637000304	SAV2	2466	2699	+	26	71	3	74	pfam04851	36.1	n/a
# 637203095	Finished	Streptomyces avermitilis MA-4680	637000121	637000304	SAV3	2765	3310	+	2	80	40	138	pfam02945	68.2	n/a
# 637203098	Finished	Streptomyces avermitilis MA-4680	637000121	637000304	SAV6	4936	7563	+	2	183	1	210	pfam04851	56.5	n/a
# 637203098	Finished	Streptomyces avermitilis MA-4680	637000121	637000304	SAV6	4936	7563	+	11	77	366	448	pfam00271	21.1	n/a
# 637203098	Finished	Streptomyces avermitilis MA-4680	637000121	637000304	SAV6	4936	7563	+	3	68	565	640	pfam03457	59	n/a

def parse_arguments():
    """Creates the argument parser"""
    
    # Create the argument parser
    parser = argparse.ArgumentParser(
        description=__description__,
        epilog=__epilog__,
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # Add mandatory arguments
    parser.add_argument("gff", metavar="*.gff", help="Prodigal/RefSeq *.gff file")
    parser.add_argument("domains", metavar="*.pfam_scan", help="pfam_scan.pl *.pfam_scan file")
    
    # Add optional arguments
    parser.add_argument('-s', '--status', type=str, default="draft", help='sequencing satus (default: draft)')
    parser.add_argument('-o', '--organism', type=str, default="unknown", help='organism name (default: unknown)')
    parser.add_argument('-c', '--scaffold_id', type=str, default="unknown", help='scaffold id (default: unknown)')
    parser.add_argument('-i', '--organism_id', type=str, default="unknown", help='organism id (default: unknown)')
    parser.add_argument('-e', '--exclude_unmatched', action='store_true', help='do not include genes without Pfam domains (default: False)')
    
    # parse the arguments
    args = parser.parse_args()
    
    return args

def convert_pfam_id(input_id):
    # removes the decimal from PFAM id
    split_id = input_id.split(".", 1)
    return split_id[0]

def parse_gff(filename):
    output = []
    with open(filename, 'r') as flo:
        for line in flo:
            if not line.startswith('#'):
                # seqid	sorce/algorithm/software	type/feature	start	end	score	strand	phase	attributes/extra
                split_line = line.rstrip().split("\t")
                if (split_line[2] == 'CDS'):
                    if split_line[1].startswith("RefSeq"):
                        split_attributes = split_line[-1].split(";")
                        tag, uid = split_attributes[0].split("=")
                        cds_num = re.findall(r'(\d+)$', uid)[0]
                        
                        rs_tag, rs_uid = split_attributes[1].split("=")
                        output.append((split_line[0] + "_" + cds_num, split_line[3], split_line[4], split_line[6], rs_uid)) # contig_id, start, end, strand, RefSeqID
                    else:
                        split_attributes = split_line[-1].split(";")
                        tag, uid = split_attributes[0].split("=")
                        cds_num = re.findall(r'(\d+)$', uid)[0]
                        output.append((split_line[0] + "_" + cds_num, split_line[3], split_line[4], split_line[6])) # contig_id, start, end, strand
    return output

def parse_domains(filename):
    output = {}
    # open "domains" file
    with open(filename, 'r') as flo:
        for line in flo:
            line = line.rstrip()
            if (not line.startswith('#')) and (len(line) > 0):
                # hmmscan file
                #                                                                                          --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
                #   target name        accession   tlen query name                       accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
                #  ------------------- ---------- -----             -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
                # ResIII               PF04851.10   184 lcl|BA000030.3_prot_BAC67711.1_2 -             77   5.4e-09   36.1   0.0   1   1   7.5e-13   5.5e-09   36.1   0.0    26    71     6    53     3    74 0.89 Type III restriction enzyme, res subunit
                # DEAD                 PF00270.24   169 lcl|BA000030.3_prot_BAC67711.1_2 -             77   7.3e-07   28.8   0.0   1   1   1.1e-10     8e-07   28.7   0.0    19    78    10    67     4    76 0.87 DEAD/DEAH box helicase
                # Endonuclease_7       PF02945.10    82 lcl|BA000030.3_prot_BAC67712.1_3 -            181   3.6e-19   68.2   1.4   1   1   2.4e-23   3.6e-19   68.2   1.4     2    80    41   136    40   138 0.87 Recombination endonuclease VII
                # DUF2249              PF10006.4     69 lcl|BA000030.3_prot_BAC67713.1_4 -            168     0.093   12.3   0.0   1   1   1.4e-05      0.21   11.1   0.0    27    48   101   122    91   146 0.80 Uncharacterized conserved protein (DUF2249)
                # zf-CGNR              PF11706.3     44 lcl|BA000030.3_prot_BAC67714.1_5 -            312       6.7    6.3  16.5   1   1    0.0056        83    2.8  16.5     2    44   235   284   234   284 0.73 CGNR zinc finger
                # 00000000000000000    111111111         33333333333333333333333333                                     7                                                    15    16               19    20
                #split_line = re.split(r'\s+', line.rstrip(), 22)
                #try:
                #    # gene, pfam_t-start, pfam_t-end, pfam_start, pfam_end, pfam_id, e-score, enzyme
                #    output[split_line[3]].append((split_line[0], split_line[15], split_line[16], split_line[19], split_line[20], split_line[1], split_line[7], 'n/a'))
                #except KeyError:
                #    output[split_line[3]] = [(split_line[0], split_line[15], split_line[16], split_line[19], split_line[20], split_line[1], split_line[7], 'n/a')]
                
                # pfam_scan.pl file
                # <seq id>              <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc>   <hmm name>        <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan>
                # gi|29826542|ref|NP_821176.1|          6              53                3             74 PF04851.10  ResIII            Family          26        71          184        36.1   5.5e-09              1 CL0023
                # gi|29826543|ref|NP_821177.1|         41             136               40            138 PF02945.10  Endonuclease_7    Domain           2        80           82        68.2   3.6e-19              1 CL0263
                # gi|29826546|ref|NP_821180.1|          2             209                1            210 PF04851.10  ResIII            Family           2       183          184        56.5   2.9e-15              1 CL0023
                # gi|29826546|ref|NP_821180.1|        370             447              366            448 PF00271.26  Helicase_C        Family          11        77           78        21.3   0.00017              1 CL0023
                # gi|29826546|ref|NP_821180.1|        567             640              565            640 PF03457.9   HA                Domain           3        68           68        59.0   3.4e-16              1 No_clan
                # gi|29826546|ref|NP_821180.1|        641             715              641            716 PF03457.9   HA                Domain           1        67           68        48.7   5.9e-13              1 No_clan
                # gi|29826546|ref|NP_821180.1|        719             785              718            786 PF03457.9   HA                Domain           2        67           68        43.6   2.3e-11              1 No_clan
                # gi|29826546|ref|NP_821180.1|        803             872              800            873 PF03457.9   HA                Domain           5        67           68        49.8   2.6e-13              1 No_clan
                
                split_line = re.split(r'\s+', line)
                try:
                    # gene, pfam_t-start, pfam_t-end, pfam_start, pfam_end, pfam_id, e-score, enzyme
                    output[split_line[0]].append((split_line[6], split_line[8], split_line[9], split_line[3], split_line[4], split_line[5], split_line[11], 'n/a'))
                except KeyError:
                    output[split_line[0]] = [(split_line[6], split_line[8], split_line[9], split_line[3], split_line[4], split_line[5], split_line[11], 'n/a')]
    return output

if (__name__ == '__main__'):
    # load arguments and parse them
    args = parse_arguments()
    
    locations = parse_gff(args.gff)
    mappings = parse_domains(args.domains)
    
    for gene in locations:
        # First we try the concatenation between the scaffold ID and the CDS number
        try:
            for domain in mappings[gene[0]]:
                # GeneID	Sequencing status	Organism name	Scaffold ID	Organism ID	Locus Tag	Gene Start	Gene End	Strand	Pfam Template Start	Pfam Template End	Pfam Start	Pfam End	PfamID	Pfam E-score	Enzyme ID
                print "\t".join([gene[0], args.status, args.organism, args.scaffold_id, args.organism_id, domain[0], gene[1], gene[2], gene[3], domain[1], domain[2], domain[3], domain[4], convert_pfam_id(domain[5]), domain[6], domain[7]])
        # If that does not exist, then we try the RefSeqID
        except KeyError:
            try:
                # We loop through all the keys for the genes with Pfam domains
                for prot_id in mappings:
                    temp = prot_id.split("|")
                    if (gene[4] == temp[3]):
                        for domain in mappings[prot_id]:
                            print "\t".join([prot_id, args.status, args.organism, args.scaffold_id, args.organism_id, domain[0], gene[1], gene[2], gene[3], domain[1], domain[2], domain[3], domain[4], convert_pfam_id(domain[5]), domain[6], domain[7]])
                        # once we find the correct gene, we can stop looking
                        break
                # If we did not find that the gene had any Pfam domains, then we print the entry with several 'n/a' fields
                else:
                    if not args.exclude_unmatched:
                        print "\t".join([gene[0], args.status, args.organism, args.scaffold_id, args.organism_id, 'n/a', gene[1], gene[2], gene[3], 'n/a', 'n/a', 'n/a', 'n/a', 'n/a', 'n/a', 'n/a'])
                    
            except IndexError:
                if not args.exclude_unmatched:
                    print "\t".join([gene[0], args.status, args.organism, args.scaffold_id, args.organism_id, 'n/a', gene[1], gene[2], gene[3], 'n/a', 'n/a', 'n/a', 'n/a', 'n/a', 'n/a', 'n/a'])
    # program end
