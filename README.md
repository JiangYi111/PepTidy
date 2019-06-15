# PepTidy
PepTidy: an automated pipeline for annotating and mining peptidases

#www.peptidy.cn

#https://github.com/JiangYi111/PepTidy/

#Version 1.0.0

#This script is under the Artistic Licence

#https://opensource.org/licenses/artistic-license-2.0


#Author and Contact Information
Yi Jiang
yijiang@peptidy.cn

#Bug reporting
bugs@peptidy.cn

#INTRODUCTION
PepTidy is used to annotate and mine peptidases. The input file of PepTidy is protein sequences in FastA format1. PepTidy can detect homologues of peptidases by sequence pairwise alignments. PepTidy can also distinguish active peptidases from non-peptidease homologues by checking the completeness of important residues for activity of peptide bond hydrolysis, such as catalytic triad and metal binding sites, which is the most important advantage over annotation of peptideases merely relying on sequence similarity. Users can mine these putative peptidases by annotation informations, such as peptidase family, E.C. number, active site residue, metal binding site residue, peptidase name and so on.

#INSTALLATION on Linux or MacOSX

BLAST+ installation
$  cd ~
$  wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/ncbi-blast-2.8.1+-x64-linux.tar.gz
$  tar -xzvf ncbi-blast-2.8.1+-x64-linux.tar.gz
$  vim .bashrc
   add 'export PATH="$HOME/ncbi-blast-2.8.1+/bin:$PATH"' (without single quotes) to the last line.
$  source .bashrc

PepTidy installation
$  cd ~
$  wegt http://www.peptidy.cn/downloads/PepTidy_1.0.0.tar.gz
$  tar -xzvf PepTidy_1.0.0.tar.gz
$  makeblastdb -in ~/PepTidy_1.0.0/MEROPS_database/protease_lib -dbtype prot -out ~/PepTidy_1.0.0/MEROPS_database/protease_lib
$  makeblastdb -in ~/PepTidy_1.0.0/MEROPS_database/pepunit_lib -dbtype prot -out ~/PepTidy_1.0.0/MEROPS_database/pepunit_lib
$  perl ~/PepTidy_1.0.0/PepTidy.pl -help

If the help message of PepTidy.pl is shown, PepTidy is installed successfully. And users can run PepTidy to mine and annotate peptidases. If reporting ‘Can’t locate xxx.pm …’, please move on to next part “perl module installation”.

Perl module installation
$  sudo cpanm Getopt::Long
$  sudo cpanm Cwd
$  sudo cpanm File::Spec::Functions
$  sudo cpanm File::Basename
$  sudo cpanm Module::Load::Conditional

#USAGE
Please refer to 'www.peptidy.cn/downloads/PepTidy_Userguide.pdf'
