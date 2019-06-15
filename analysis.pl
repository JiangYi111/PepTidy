#! /usr/bin/perl

eval 'exec perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

#############################################################
#                                                           #
# analysis.pl                                               #
# analyzing PepTidy annotated resluts                       #
#                                                           #
# Authors: Yi Jiang                                         #
#                                                           #
# Contact   : yijiang@peptidy.cn                            #
# Bug report: bugs@peptidy.cn                               #
#                                                           #
# Release date: April 15th, 2019                            #
#                                                           #
# This script is under the Artistic Licence                 #
# https://opensource.org/licenses/artistic-license-2.0      #
#                                                           #
# Usage:                                                    #
# perl analysis.pl [option] --in                            #
#                                                           #
#############################################################

use Getopt::Long;
use Module::Load::Conditional qw(can_load check_install requires);
use Cwd;
use File::Spec::Functions qw(rel2abs);
use File::Basename qw(dirname basename);

$0 = rel2abs($0);
my $directory = dirname($0);

use strict;
use warnings;

my $usage = <<'ENDUSAGE';

analysis.pl is to analyze PepTidy annotated resluts.

SYNOPSIS

perl analysis.pl [option] --in

  --help                      To print this help message.
                              --help
  --in                        Path to directory containing Peptidy resluts files.
                              --in <String>
                              Example [--in /path/to/dir_with_Peptidy_resluts_files]
  --out                       Path to directory saving analysis results.
                              --out <String>
                              Example [--out /path/to/out_dir]
                              Default = '/path/to/dir_with_Peptidy_resluts_files/result_analysis'
  --statistics                To save statistical informations of annotated resluts to file 'statistics'.
                              --statistics                    
  --custom                    If resluts files are annotated by custom database, please set it.
                              --custom
  --db_build                  To build databases for PepTidy mining.
                              --db_build
  --db_name                   Name of database which will be builded. Start with letter. And only number,
                              letter and underline are vaild.
                              --db_name <String>
                              Example [--db_name name]
  Options to list amino acid of key sites and save to file 'actsites'.
  ------------------------------------------------------------------------------------------------------
  --active_list               To list amino acids of all active sites.
                              --active_list
  --metal_list                To list amino acids of all metal sites.
                              --metal_list
  --binding_list              To list amino acids of all binding sites.
                              --binding_list
  --active_metal_list         To list amino acids of all active and metal sites.
                              --active_metal_list
  and so on [--active_binding_list] 
            [--metal_binding_list] 
            [--active_metal_binding_list]

  Options to list E.C. number, peptidase family, protein name, query name and database corss-references.
  ------------------------------------------------------------------------------------------------------
  --EC_list                   To list all E.C. numbers and save to file 'EC_number'.
                              --EC_list
  --family_list               To list all peptidase families and save to file 'family'.
                              --family_list
  --name_list                 To list all annotated protein names and save to file 'name'.
                              --name_list
  --query_list                To list all annotated query names and save to file 'query'.
                              --query_list
  --annotation_list           To list all database corss-references from Swiss-Prot and save to file 'annotation'.
                              --annotation_list

  Options to search Peptidy annotated results by criteria. Matched results are saved to file 'filter'.
  ------------------------------------------------------------------------------------------------------
  --query                     To include the query.
                              --query <String>
                              Example [--query query1 --query query2 ... --query queryn]
  --de_query                  To exclude the query.
                              --de_query <String>
                              Example [--de_query query1 --de_query query2 ... --de_query queryn]
  --EC                        To include the EC number.
                              --EC <String>
                              Example [--EC EC1 --EC EC2 ... --EC ECn]
  --de_EC                     To exclude the EC number.
                              --de_EC <String>
                              Example [--de_EC EC1 --de_EC EC2 ... --de_EC ECn]
  --family                    To include the peptidase family.
                              --family <String>
                              Example [--family family1 --family family2 ... --family familyn]
                              If you want search all subfamilies, please use 'A00' or 'M00' or 'C00' or
                              'S00' or 'T00' or 'P00' or 'N00' or 'U00' or 'G00'.
                              such as [--family A11 --family C00 --family M1]
  --de_family                 To exclude the peptidase family.
                              --de_family <String>
                              Example [--de_family family1 --de_family family2 ... --de_family familyn]
  --name                      To include the annotated protein name.
                              --name <String>
                              Example [--name name1 --name name2 ... --name namen]
                              If blank character( ) in name, please change it to comma(,).
                              If underline(_) in name , please change it to comma(,).
                              (changing 'Carboxypeptidase Y' to 'Carboxypeptidase,Y')
                              (changing 'ADD_CBB' to 'ADD,CBB')
                              such as [--name extracellular,protease --name Prolyl,aminopeptidase
                                       --name carboxypeptidase --name KEX2 --name Pepsin,A ...] 
  --de_name                   To exclude the annotated protein name.
                              --de_name <String>
                              Example [--de_name name1 --de_name name2 ... --de_name namen]
  --active                    To include the amino acids of active sites.
                              --active <String>
                              Example [--active active1 --active active2 ... --active activen]
                              such as [--active HDS --active HC --active DD ...]
  --de_active                 To exclude the amino acids of active sites.
                              --de_active <String>
                              Example [--de_active active1 --de_active active2 ... --de_active activen]
  --metal                     To include the amino acids of metal sites.
                              --metal <String>
                              Example [--metal metal1 --metal metal2 ... --metal metaln]
  --de_metal                  To exclude the amino acids of metal sites.
                              --de_metal <String>
                              Example [--de_metal metal1 --de_metal metal2 ... --de_metal metaln]
  --binding                   To include the amino acids of binding sites.
                              --binding <String>
                              Example [--binding binding1 --binding binding2 ... --binding bindingn]
  --de_binding                To exclude the amino acids of binding sites.
                              --de_binding <String>
                              Example [--de_binding binding1 --de_binding binding2 ... --de_binding bindingn]
  --active_metal              To include the amino acids of active and metal sites.
                              --active_metal <String>
                              Example [--active_metal active_metal1 ... --active_metal active_metaln]
                              such as [--active_metal HE_HHD --active_metal HDC_HEHE ...]
  --de_active_metal           To exclude the amino acids of active and metal sites.
                              --de_active_metal <String>
                              Example [--de_active_metal active_metal1 ... --de_active_metal active_metaln]
  and so on [--active_binding]       [--de_active_binding] 
            [--metal_binding]        [--de_metal_binding]       
            [--active_metal_binding] [--de_active_metal_binding]

  --EC_family                 To include the EC number and the peptidase family.
                              --EC_family <String>
                              Example [--EC_family EC_family1 --EC_family EC_family2 ... --EC_family EC_familyn]
                              such as [--EC_family 3.4.11.5_S33 --EC_family 3.4.23.-_A11 --EC_family 3.4.-.5_S8]
  --de_EC_family              To exclude the EC number and the peptidase family.
                              Example [--de_EC_family EC_family1 ... --de_EC_family EC_familyn]
  and so on [--EC_name]                             [--de_EC_name] 
            [--EC_active]                           [--de_EC_active] 
            [--EC_metal]                            [--de_EC_metal]
            [--EC_binding]                          [--de_EC_binding] 
            [--EC_active_metal]                     [--de_EC_active_metal] 
            [--EC_active_binding]                   [--de_EC_active_binding] 
            [--EC_metal_binding]                    [--de_EC_metal_binding]
            [--EC_active_metal_binding]             [--de_EC_active_metal_binding] 
            [--family_name]                         [--de_family_name]
            [--family_active]                       [--de_family_active] 
            [--family_metal]                        [--de_family_metal] 
            [--family_binding]                      [--de_family_binding] 
            [--family_active_metal]                 [--de_family_active_metal] 
            [--family_active_binding]               [--de_family_active_binding] 
            [--family_metal_binding]                [--de_family_metal_binding] 
            [--family_active_metal_binding]         [--de_family_active_metal_binding] 
            [--name_active]                         [--de_name_active]
            [--name_metal]                          [--de_name_metal] 
            [--name_binding]                        [--de_name_binding] 
            [--name_active_metal]                   [--de_name_active_metal] 
            [--name_active_binding]                 [--de_name_active_binding] 
            [--name_metal_binding]                  [--de_name_metal_binding] 
            [--name_active_metal_binding]           [--de_name_active_metal_binding] 
            [--EC_family_active]                    [--de_EC_family_active] 
            [--EC_family_metal]                     [--de_EC_family_metal]
            [--EC_family_binding]                   [--de_EC_family_binding] 
            [--EC_family_active_metal]              [--de_EC_family_active_metal] 
            [--EC_family_active_binding]            [--de_EC_family_active_binding] 
            [--EC_family_metal_binding]             [--de_EC_family_metal_binding] 
            [--EC_family_active_metal_binding]      [--de_EC_family_active_metal_binding] 
            [--EC_name_active]                      [--de_EC_name_active] 
            [--EC_name_metal]                       [--de_EC_name_metal] 
            [--EC_name_binding]                     [--de_EC_name_binding] 
            [--EC_name_active_metal]                [--de_EC_name_active_metal] 
            [--EC_name_active_binding]              [--de_EC_name_active_binding] 
            [--EC_name_metal_binding]               [--de_EC_name_metal_binding] 
            [--EC_name_active_metal_binding]        [--de_EC_name_active_metal_binding] 
            [--family_name_active]                  [--de_family_name_active] 
            [--family_name_metal]                   [--de_family_name_metal] 
            [--family_name_binding]                 [--de_family_name_binding] 
            [--family_name_active_metal]            [--de_family_name_active_metal] 
            [--family_name_active_binding]          [--de_family_name_active_binding] 
            [--family_name_metal_binding]           [--de_family_name_metal_binding] 
            [--family_name_active_metal_binding]    [--de_family_name_active_metal_binding] 
            [--EC_family_name_active]               [--de_EC_family_name_active] 
            [--EC_family_name_metal]                [--de_EC_family_name_metal] 
            [--EC_family_name_binding]              [--de_EC_family_name_binding] 
            [--EC_family_name_active_metal]         [--de_EC_family_name_active_metal] 
            [--EC_family_name_active_binding]       [--de_EC_family_name_active_binding] 
            [--EC_family_name_metal_binding]        [--de_EC_family_name_metal_binding] 
            [--EC_family_name_active_metal_binding] [--de_EC_family_name_active_metal_binding] 

  --evalue                    To include scope of E-value.
                              --evalue <String>
                              Example [--evalue evalue_small#evalue_big] 
                              such as [--evalue 1e-100#1] including evalue from 1e-100 to 1
                              such as [--evalue 1e-5#1e-5] including evalue equals to 1e-5
                              Example [--evalue Bevalue_small]
                                      ('Bevalue_small' means big than evalue_small)
                              such as [--evalue B0]  including evalue big than 0 and including 0
                              such as [--evalue B1e-100] including evalue big than 1e-100 
                                                         and including 1e-100
                              Example [--evalue Sevalue_big]
                                      ('Sevalue_big' means small than evalue_big)
                              such as [--evalue S0.5]  including evalue small than 0.5 and including 0.5 
                              Example [--evalue same]
                              if '--evalue' is set to 'same', queries included are same with sequences 
                              in MEROPS database or custom database.
                              such as [--evalue same]
                              !!!Note: 1e-324 and numbers small than it are equal to 0.
  --all                       To print all iterms of matched queries.
                              --all

  Options to obtain informations from uniprot database. Results are saved to file 'filter'.
  ------------------------------------------------------------------------------------------------------
  --reference                 To obtain reference informations from uniprot database.
                              including PubMed id, DOI, paper title, journal name.
                              --reference
  --note                      To obtain comments or notes from uniprot database.
                              --note
  --annotation                To obtain database corss-references from uniprot database.
                              including KEGG, KO, GO, MEROPS, eggNOG, Pfam, OrthoDB, InterPro and so on.
                              --annotation <String>
                              Example [--annotation annotation1 ... --annotation annotationn]
                              such as [--annotation KEGG --annotation KO --annotation GO ...]
                              If you want obtain all database corss-references from uniprot database,
                              please use 'ALLADD' as the argument of '--annotation'.
                              such as [--annotation ALLADD]
  --active_site               To obtain comments of active sites from uniprot database.
                              --active_site
  --metal_site                To obtain comments of metal sites from uniprot database.
                              --metal_site
  --binding_site              To obtain comments of binding sites from uniprot database.
                              --binding_site


   !!!NOTE: The priority of '--de_XXX' are higher than '--XXX'.
   !!!NOTE: If resluts files are annotated by custom database, please set '--custom'.

perl analysis.pl [option] --in

DESCRIPTION

  Example:
    
    perl analysis.pl [option] --in

ENDUSAGE

my $help;                                  # to print usage
my $currentDir = cwd();                    # working superdirectory where programme is called from
my $workDir;                               # in the working directory results and temporary files are stored
my $statistics;                            # to save statistical informations of annotated resluts to file 'statistics'
my $output;                                # path to directory saving analysis results
my $input;                                 # path to Peptidy annotated resluts files
my $annotation_seq;                        # annotated query sequences file
my $annotation_query;                      # annotated results
my $annotation_uniprot;                    # annotated results 
my %query_name_hash;                       # query names hash
my %bais_hash;                             # bias query names hash;
my %query_result_all_hash;
my %query_same_hash;
my %query_family_hash;
my $evalue_max;
my $evalue_min;
my $custom;
my $db_build;                              # to build databases for PepTidy mining
my $db_name;                               # database name 
my $database_seq;
my $database_actsite;
my $database_AMB;
my $database_length;
my $database_name;
my $database_uniprot;
my %db_query_name_hash;
# 1e-324 and numbers small than it are equal to 0
my %uniprot_name_hash;
my $active_list;                           # to list amino acids of all active sites
my %active_list_hash;
my $metal_list;                            # to list amino acids of all metal sites
my %metal_list_hash;
my $binding_list;                          # to list amino acids of all binding sites
my %binding_list_hash;
my $active_metal_list;                     # to list amino acids of all active and metal sites
my %active_metal_list_hash;
my $active_binding_list;                   # to list amino acids of all active and binding sites
my %active_binding_list_hash;
my $metal_binding_list;                    # to list amino acids of all metal and binding sites
my %metal_binding_list_hash;
my $active_metal_binding_list;             # to list amino acids of all active and metal and binding sites
my %active_metal_binding_list_hash;
my $EC_list;                               # to list all E.C. numbers and save to file 'EC_number'
my %EC_list_hash;
my $family_list;                           # to list all peptidase families and save to file 'family'
my %family_list_hash;
my $name_list;                             # to list all annotated protein names and save to file 'name'
my %name_list_hash;
my $query_list;                            # to list all annotated query names and save to file 'query'
my $annotation_list;
my %annotation_list_hash;                  # to list all database corss-references from Swiss-Prot

my %query_name_search_all_hash;
my @query_array;                           # to include the query
my %query_hash;
my @de_query_array;                        # to exclude the query
my @EC_array;                              # to include the EC number
my %EC_hash;
my @de_EC_array;                           # to exclude the EC number
my %de_EC_hash;
my @family_array;                          # to include the peptidase family
my %family_hash;
my @de_family_array;                       # to exclude the peptidase family
my %de_family_hash;
my @name_array;                            # to include the annotated protein name
my %name_hash;
my @de_name_array;                         # to exclude the annotated protein name
my %de_name_hash;
my @active_array;                          # to include the amino acids of active sites
my %active_hash;
my @de_active_array;                       # to exclude the amino acids of active sites
my %de_active_hash;
my @metal_array;                           # to include the amino acids of metal sites
my %metal_hash;
my @de_metal_array;                        # to exclude the amino acids of metal sites
my %de_metal_hash;
my @binding_array;                         # to include the amino acids of binding sites
my %binding_hash;
my @de_binding_array;                      # to exclude the amino acids of binding sites
my %de_binding_hash;
my @active_metal_array;
my %active_metal_hash;
my @de_active_metal_array;
my %de_active_metal_hash;
my @active_binding_array;
my %active_binding_hash;
my @de_active_binding_array;
my %de_active_binding_hash;
my @metal_binding_array;
my %metal_binding_hash;
my @de_metal_binding_array;
my %de_metal_binding_hash;
my @active_metal_binding_array;
my %active_metal_binding_hash;
my @de_active_metal_binding_array;
my %de_active_metal_binding_hash;
my @EC_family_array;
my %EC_family_hash;
my @de_EC_family_array;
my %de_EC_family_hash;
my @EC_name_array;
my %EC_name_hash;
my @de_EC_name_array;
my %de_EC_name_hash;
my @EC_active_array;
my %EC_active_hash;
my @de_EC_active_array;
my %de_EC_active_hash;
my @EC_metal_array;
my %EC_metal_hash;
my @de_EC_metal_array;
my %de_EC_metal_hash;
my @EC_binding_array;
my %EC_binding_hash;
my @de_EC_binding_array;
my %de_EC_binding_hash;
my @EC_active_metal_array;
my %EC_active_metal_hash;
my @de_EC_active_metal_array;
my %de_EC_active_metal_hash;
my @EC_active_binding_array;
my %EC_active_binding_hash;
my @de_EC_active_binding_array;
my %de_EC_active_binding_hash;
my @EC_metal_binding_array;
my %EC_metal_binding_hash;
my @de_EC_metal_binding_array;
my %de_EC_metal_binding_hash;
my @EC_active_metal_binding_array;
my %EC_active_metal_binding_hash;
my @de_EC_active_metal_binding_array;
my %de_EC_active_metal_binding_hash;
my @family_name_array;
my %family_name_hash;
my @de_family_name_array;
my %de_family_name_hash;
my @family_active_array;
my %family_active_hash;
my @de_family_active_array;
my %de_family_active_hash;
my @family_metal_array;            
my %family_metal_hash;             
my @de_family_metal_array;         
my %de_family_metal_hash; 
my @family_binding_array;
my %family_binding_hash;
my @de_family_binding_array;
my %de_family_binding_hash;
my @family_active_metal_array;            
my %family_active_metal_hash;             
my @de_family_active_metal_array;         
my %de_family_active_metal_hash;          
my @family_active_binding_array;            
my %family_active_binding_hash;             
my @de_family_active_binding_array;         
my %de_family_active_binding_hash;
my @family_metal_binding_array;
my %family_metal_binding_hash;
my @de_family_metal_binding_array;
my %de_family_metal_binding_hash;
my @family_active_metal_binding_array;
my %family_active_metal_binding_hash;
my @de_family_active_metal_binding_array;
my %de_family_active_metal_binding_hash;
my @name_active_array;
my %name_active_hash;
my @de_name_active_array;
my %de_name_active_hash;
my @name_metal_array;
my %name_metal_hash;
my @de_name_metal_array;
my %de_name_metal_hash;
my @name_binding_array;
my %name_binding_hash;
my @de_name_binding_array;
my %de_name_binding_hash;
my @name_active_metal_array;
my %name_active_metal_hash;
my @de_name_active_metal_array;
my %de_name_active_metal_hash;
my @name_active_binding_array;
my %name_active_binding_hash;
my @de_name_active_binding_array;
my %de_name_active_binding_hash;
my @name_metal_binding_array;
my %name_metal_binding_hash;
my @de_name_metal_binding_array;
my %de_name_metal_binding_hash;
my @name_active_metal_binding_array;
my %name_active_metal_binding_hash;
my @de_name_active_metal_binding_array;
my %de_name_active_metal_binding_hash;
my @EC_family_active_array;
my %EC_family_active_hash;
my @de_EC_family_active_array;
my %de_EC_family_active_hash;
my @EC_family_metal_array;
my %EC_family_metal_hash;
my @de_EC_family_metal_array;
my %de_EC_family_metal_hash;
my @EC_family_binding_array;
my %EC_family_binding_hash;
my @de_EC_family_binding_array;
my %de_EC_family_binding_hash;
my @EC_family_active_metal_array;
my %EC_family_active_metal_hash;
my @de_EC_family_active_metal_array;
my %de_EC_family_active_metal_hash;
my @EC_family_active_binding_array;
my %EC_family_active_binding_hash;
my @de_EC_family_active_binding_array;
my %de_EC_family_active_binding_hash;
my @EC_family_metal_binding_array;
my %EC_family_metal_binding_hash;
my @de_EC_family_metal_binding_array;
my %de_EC_family_metal_binding_hash;
my @EC_family_active_metal_binding_array;
my %EC_family_active_metal_binding_hash;
my @de_EC_family_active_metal_binding_array;
my %de_EC_family_active_metal_binding_hash;
my @EC_name_active_array;
my %EC_name_active_hash;
my @de_EC_name_active_array;
my %de_EC_name_active_hash;
my @EC_name_metal_array;
my %EC_name_metal_hash;
my @de_EC_name_metal_array;
my %de_EC_name_metal_hash;
my @EC_name_binding_array;
my %EC_name_binding_hash;
my @de_EC_name_binding_array;
my %de_EC_name_binding_hash;
my @EC_name_active_metal_array;
my %EC_name_active_metal_hash;
my @de_EC_name_active_metal_array;
my %de_EC_name_active_metal_hash;
my @EC_name_active_binding_array;
my %EC_name_active_binding_hash;
my @de_EC_name_active_binding_array;
my %de_EC_name_active_binding_hash;
my @EC_name_metal_binding_array;
my %EC_name_metal_binding_hash;
my @de_EC_name_metal_binding_array;
my %de_EC_name_metal_binding_hash;
my @EC_name_active_metal_binding_array;
my %EC_name_active_metal_binding_hash;
my @de_EC_name_active_metal_binding_array;
my %de_EC_name_active_metal_binding_hash;
my @family_name_active_array;
my %family_name_active_hash;
my @de_family_name_active_array;
my %de_family_name_active_hash;
my @family_name_metal_array;
my %family_name_metal_hash;
my @de_family_name_metal_array;
my %de_family_name_metal_hash;
my @family_name_binding_array;
my %family_name_binding_hash;
my @de_family_name_binding_array;
my %de_family_name_binding_hash;
my @family_name_active_metal_array;
my %family_name_active_metal_hash;
my @de_family_name_active_metal_array;
my %de_family_name_active_metal_hash;
my @family_name_active_binding_array;
my %family_name_active_binding_hash;
my @de_family_name_active_binding_array;
my %de_family_name_active_binding_hash;
my @family_name_metal_binding_array;
my %family_name_metal_binding_hash;
my @de_family_name_metal_binding_array;
my %de_family_name_metal_binding_hash;
my @family_name_active_metal_binding_array;
my %family_name_active_metal_binding_hash;
my @de_family_name_active_metal_binding_array;
my %de_family_name_active_metal_binding_hash;
my @EC_family_name_active_array;
my %EC_family_name_active_hash;
my @de_EC_family_name_active_array;
my %de_EC_family_name_active_hash;
my @EC_family_name_metal_array;
my %EC_family_name_metal_hash;
my @de_EC_family_name_metal_array;
my %de_EC_family_name_metal_hash;
my @EC_family_name_binding_array;
my %EC_family_name_binding_hash;
my @de_EC_family_name_binding_array;
my %de_EC_family_name_binding_hash;
my @EC_family_name_active_metal_array;
my %EC_family_name_active_metal_hash;
my @de_EC_family_name_active_metal_array;
my %de_EC_family_name_active_metal_hash;
my @EC_family_name_active_binding_array;
my %EC_family_name_active_binding_hash;
my @de_EC_family_name_active_binding_array;
my %de_EC_family_name_active_binding_hash;
my @EC_family_name_metal_binding_array;
my %EC_family_name_metal_binding_hash;
my @de_EC_family_name_metal_binding_array;
my %de_EC_family_name_metal_binding_hash;
my @EC_family_name_active_metal_binding_array;
my %EC_family_name_active_metal_binding_hash;
my @de_EC_family_name_active_metal_binding_array;
my %de_EC_family_name_active_metal_binding_hash;
my %include_hash;
my %exclude_hash;


my $evalue;
my %query_name_evalue_hash;
my $ALL;
my $reference;
my $note;
my $active_site;
my $metal_site;
my $binding_site;
my @annotation_array;
my %annotation_hash;
my $annotation_all;

my $count_QN_max = 1;
my %uniprot_id_hash;
my %reference_hash;
my %note_hash;
my %annotation_cross_hash;
my %active_site_hash;
my %metal_site_hash;
my %binding_site_hash;


my %family_sum_hash = ('A' => 'aspartic peptidases family',
                       'S' => 'serine peptidases family',
                       'M' => 'metallo peptidases family',
                       'C' => 'cysteine peptidases family',
                       'N' => 'asparagine peptide lyases family',
                       'G' => 'glutamic peptidases family',
                       'U' => 'peptidases of unknown catalytic type',
                       'T' => 'threonine peptidases family',
                       'P' => 'mixed peptidases family');

if(@ARGV==0){
  print "$usage\n";
  exit(0);
}

GetOptions( 'help!'                         => \$help,
            'in=s'                          => \$input,
            'statistics!'                   => \$statistics,
            'active_list!'                  => \$active_list,
            'metal_list!'                   => \$metal_list,
            'binding_list!'                 => \$binding_list,
            'active_metal_list!'            => \$active_metal_list,
            'active_binding_list!'          => \$active_binding_list,
            'metal_binding_list!'           => \$metal_binding_list,
            'active_metal_binding_list!'    => \$active_metal_binding_list,
            'EC_list!'                      => \$EC_list,
            'family_list!'                  => \$family_list,
            'name_list!'                    => \$name_list,
            'query_list!'                   => \$query_list, 
            'annotation_list!'              => \$annotation_list,
            'query=s'                       => \@query_array,
            'de_query=s'                    => \@de_query_array,
            'EC=s'                          => \@EC_array,
            'de_EC=s'                       => \@de_EC_array,
            'family=s'                      => \@family_array,
            'de_family=s'                   => \@de_family_array,
            'name=s'                        => \@name_array,
            'de_name=s'                     => \@de_name_array,
            'active=s'                      => \@active_array,
            'de_active=s'                   => \@de_active_array,
            'metal=s'                       => \@metal_array,
            'de_metal=s'                    => \@de_metal_array,
            'binding=s'                     => \@binding_array,
            'de_binding=s'                  => \@de_binding_array,
            'active_metal=s'                => \@active_metal_array,               
            'de_active_metal=s'             => \@de_active_metal_array,
            'active_binding=s'              => \@active_binding_array,
            'de_active_binding=s'           => \@de_active_binding_array,
            'metal_binding=s'               => \@metal_binding_array,
            'de_metal_binding=s'            => \@de_metal_binding_array,
            'active_metal_binding=s'        => \@active_metal_binding_array,
            'de_active_metal_binding=s'     => \@de_active_metal_binding_array,
            'EC_family=s'                   => \@EC_family_array,
            'de_EC_family=s'                => \@de_EC_family_array,
            'EC_name=s'                     => \@EC_name_array,
            'de_EC_name=s'                  => \@de_EC_name_array,
            'EC_active=s'                   => \@EC_active_array,
            'de_EC_active=s'                => \@de_EC_active_array,
            'EC_metal=s'                    => \@EC_metal_array,
            'de_EC_metal=s'                 => \@de_EC_metal_array,
            'EC_binding=s'                  => \@EC_binding_array,
            'de_EC_binding=s'               => \@de_EC_binding_array,
            'EC_active_metal=s'             => \@EC_active_metal_array,
            'de_EC_active_metal=s'          => \@de_EC_active_metal_array,
            'EC_active_binding=s'           => \@EC_active_binding_array,
            'de_EC_active_binding=s'        => \@de_EC_active_binding_array,
            'EC_metal_binding=s'            => \@EC_metal_binding_array,
            'de_EC_metal_binding=s'         => \@de_EC_metal_binding_array,
            'EC_active_metal_binding=s'     => \@EC_active_metal_binding_array,
            'de_EC_active_metal_binding=s'  => \@de_EC_active_metal_binding_array,
            'family_name=s'                 => \@family_name_array,
            'de_family_name=s'              => \@de_family_name_array,
            'family_active=s'               => \@family_active_array,
            'de_family_active=s'            => \@de_family_active_array,
            'family_metal=s'                => \@family_metal_array,
            'de_family_metal=s'             => \@de_family_metal_array,
            'family_binding=s'              => \@family_binding_array,
            'de_family_binding=s'           => \@de_family_binding_array,
            'family_active_metal=s'         => \@family_active_metal_array,
            'de_family_active_metal=s'      => \@de_family_active_metal_array,
            'family_active_binding=s'       => \@family_active_binding_array,
            'de_family_active_binding=s'    => \@de_family_active_binding_array,
            'family_metal_binding=s'        => \@family_metal_binding_array,
            'de_family_metal_binding=s'     => \@de_family_metal_binding_array,
            'family_active_metal_binding=s' => \@family_active_metal_binding_array,
            'de_family_active_metal_binding=s' => \@de_family_active_metal_binding_array,
            'name_active=s'                 => \@name_active_array,
            'de_name_active=s'              => \@de_name_active_array,
            'name_metal=s'                  => \@name_metal_array,
            'de_name_metal=s'               => \@de_name_metal_array,
            'name_binding=s'                => \@name_binding_array,
            'de_name_binding=s'             => \@de_name_binding_array,
            'name_active_metal=s'           => \@name_active_metal_array,
            'de_name_active_metal=s'        => \@de_name_active_metal_array,
            'name_active_binding=s'         => \@name_active_binding_array,
            'de_name_active_binding=s'      => \@de_name_active_binding_array,
            'name_metal_binding=s'          => \@name_metal_binding_array,
            'de_name_metal_binding=s'       => \@de_name_metal_binding_array,
            'name_active_metal_binding=s'   => \@name_active_metal_binding_array,
            'de_name_active_metal_binding=s' => \@de_name_active_metal_binding_array,
            'EC_family_active=s'            => \@EC_family_active_array,
            'de_EC_family_active=s'         => \@de_EC_family_active_array,
            'EC_family_metal=s'             => \@EC_family_metal_array,
            'de_EC_family_metal=s'          => \@de_EC_family_metal_array,
            'EC_family_binding=s'           => \@EC_family_binding_array,
            'de_EC_family_binding=s'        => \@de_EC_family_binding_array,
            'EC_family_active_metal=s'      => \@EC_family_active_metal_array,
            'de_EC_family_active_metal=s'   => \@de_EC_family_active_metal_array,
            'EC_family_active_binding=s'    => \@EC_family_active_binding_array,
            'de_EC_family_active_binding=s' => \@de_EC_family_active_binding_array,
            'EC_family_metal_binding=s'     => \@EC_family_metal_binding_array,
            'de_EC_family_metal_binding=s'  => \@de_EC_family_metal_binding_array,
            'EC_family_active_metal_binding=s'    => \@EC_family_active_metal_binding_array,
            'de_EC_family_active_metal_binding=s' => \@de_EC_family_active_metal_binding_array,
            'EC_name_active=s'              => \@EC_name_active_array,
            'de_EC_name_active=s'           => \@de_EC_name_active_array,
            'EC_name_metal=s'               => \@EC_name_metal_array,
            'de_EC_name_metal=s'            => \@de_EC_name_metal_array,
            'EC_name_binding=s'             => \@EC_name_binding_array,
            'de_EC_name_binding=s'          => \@de_EC_name_binding_array,
            'EC_name_active_metal=s'        => \@EC_name_active_metal_array,
            'de_EC_name_active_metal=s'     => \@de_EC_name_active_metal_array,
            'EC_name_active_binding=s'      => \@EC_name_active_binding_array,
            'de_EC_name_active_binding=s'   => \@de_EC_name_active_binding_array,
            'EC_name_metal_binding=s'       => \@EC_name_metal_binding_array,
            'de_EC_name_metal_binding=s'    => \@de_EC_name_metal_binding_array,
            'EC_name_active_metal_binding=s'    => \@EC_name_active_metal_binding_array,
            'de_EC_name_active_metal_binding=s' => \@de_EC_name_active_metal_binding_array,
            'family_name_active=s'          => \@family_name_active_array,
            'de_family_name_active=s'       => \@de_family_name_active_array,
            'family_name_metal=s'           => \@family_name_metal_array,
            'de_family_name_metal=s'        => \@de_family_name_metal_array,
            'family_name_binding=s'         => \@family_name_binding_array,
            'de_family_name_binding=s'      => \@de_family_name_binding_array,
            'family_name_active_metal=s'    => \@family_name_active_metal_array,
            'de_family_name_active_metal=s' => \@de_family_name_active_metal_array,
            'family_name_active_binding=s'  => \@family_name_active_binding_array,
            'de_family_name_active_binding=s'       => \@de_family_name_active_binding_array,
            'family_name_metal_binding=s'   => \@family_name_metal_binding_array,
            'de_family_name_metal_binding=s'        => \@de_family_name_metal_binding_array,
            'family_name_active_metal_binding=s'    => \@family_name_active_metal_binding_array,
            'de_family_name_active_metal_binding=s' => \@de_family_name_active_metal_binding_array,
            'EC_family_name_active=s'       => \@EC_family_name_active_array,
            'de_EC_family_name_active=s'    => \@de_EC_family_name_active_array,
            'EC_family_name_metal=s'        => \@EC_family_name_metal_array,
            'de_EC_family_name_metal=s'     => \@de_EC_family_name_metal_array,
            'EC_family_name_binding=s'      => \@EC_family_name_binding_array,
            'de_EC_family_name_binding=s'   => \@de_EC_family_name_binding_array,
            'EC_family_name_active_metal=s' => \@EC_family_name_active_metal_array,
            'de_EC_family_name_active_metal=s'         => \@de_EC_family_name_active_metal_array,
            'EC_family_name_active_binding=s'          => \@EC_family_name_active_binding_array,
            'de_EC_family_name_active_binding=s'       => \@de_EC_family_name_active_binding_array,
            'EC_family_name_metal_binding=s'           => \@EC_family_name_metal_binding_array,
            'de_EC_family_name_metal_binding=s'        => \@de_EC_family_name_metal_binding_array,
            'EC_family_name_active_metal_binding=s'    => \@EC_family_name_active_metal_binding_array,
            'de_EC_family_name_active_metal_binding=s' => \@de_EC_family_name_active_metal_binding_array,
            'annotation=s'                  => \@annotation_array,
            'evalue=s'                      => \$evalue,
            'all!'                          => \$ALL,
            'reference!'                    => \$reference,
            'note!'                         => \$note,
            'active_site!'                  => \$active_site,
            'metal_site!'                   => \$metal_site,
            'binding_site!'                 => \$binding_site,
            'custom!'                       => \$custom,
            'db_build!'                     => \$db_build,
            'db_name=s'                     => \$db_name,
            'out=s'                         => \$output) or die();

if($help){
  print $usage;
  exit(0);
}


######################## make some regular checks ########################
if($input){
  $input = rel2abs($input);
  if(! -e $input){
    print STDERR "ERROR: directory \'$input\' does not exist. Please reset option '--in'.\n$usage\n";
    exit(1);
  }
  if(! -d $input){
    print STDERR "ERROR: directory \'$input\' does not exist. Please reset option '--in'.\n$usage\n";
    exit(1);
  }
}else{
  print STDERR "ERROR: Please set option '--in'.\n$usage\n";
  exit(1);
}


if($output){
  $output = rel2abs($output);
  if(! -e $output){
    print STDERR "ERROR: directory \'$output\' does not exist. Please reset option '--out'.\n$usage\n";
    exit(1);
  }
  if(! -d $output){
    print STDERR "ERROR: directory \'$output\' does not exist. Please reset option '--out'.\n$usage\n";
    exit(1);
  }
  if(! -w $output){
    print STDERR "ERROR: Output directory \'$output\' does not have the right to write. Please check.\n";
    exit(1);
  }
}else{
  $output = $currentDir;
  if(! -w $output){
    print STDERR "ERROR: Output directory \'$output\' does not have the right to write. Please check.\n";
    exit(1);
  }
}

my $next_switch = 1;
my $count_next;
while($next_switch == 1){
  $count_next++;
  my $output_tmp = "$output/analysis_output_dir$count_next";
  if(! -e $output_tmp){
    my $CMD_mkdir = "mkdir $output_tmp";
    system("$CMD_mkdir")== 0 or die("failed to execute command \'$CMD_mkdir\': $!\n");
    $output = $output_tmp;
    $next_switch = 2;
    last;
  }
}

if($next_switch == 2){
  $annotation_seq = "$output/tmp_analysis_annotation_seq";
  $annotation_query = "$output/tmp_analysis_annotation_query";
  $annotation_uniprot = "$output/tmp_analysis_annotation_uniprot";
  if(-e $annotation_seq || -e $annotation_query || -e $annotation_uniprot){
    print STDERR "FATAL ERROR: Please redownload \'analysis.pl\'.\n";
    exit(1);
  }
  my $CMD_cp = "cp $input/*_annotated_query_seq $annotation_seq";
  system("$CMD_cp")== 0 or die("failed to execute command \'$CMD_cp\': $!\n");
  $CMD_cp = "cp $input/*_annotation_results_query $annotation_query";
  system("$CMD_cp")== 0 or die("failed to execute command \'$CMD_cp\': $!\n");
  $CMD_cp = "cp $input/*_annotation_results_uniport $annotation_uniprot";
  system("$CMD_cp")== 0 or die("failed to execute command \'$CMD_cp\': $!\n");
}else{
  print STDERR "FATAL ERROR: Please redownload \'analysis.pl\'.\n";
  exit(1);
}


if(@annotation_array){
  my $count_tmp = scalar(@annotation_array);
  foreach my $anno_tmp(@annotation_array){
    if($count_tmp == 1){
      if($anno_tmp =~ /^ALLADD$/i){
        $annotation_all = 1
      }else{
        $annotation_hash{$anno_tmp} += 1;
      }
    }else{
      $annotation_hash{$anno_tmp} += 1;
    }
  }
}


my $evalue_small;
my $evalue_big;
my $evalue_same;
my $small_zero = 1;
my $big_zero = 1;
my $zero_key = 1;

#query_name_evalue_hash
if($evalue){
  if($evalue =~ /^([^\#]+)#([^\#]+)$/){
    $evalue_small = $1;
    $evalue_big = $2;
    unless($evalue_big){
      if($evalue_small){
        print STDERR "ERROR: The minimum E_value is 0. Please change it.\n$usage\n";
        exit(1);
      }
      unless($evalue_small){
        $zero_key = 3;
      }
    }
    if($evalue_big){
      unless($evalue_small){
        $small_zero = 2; 
      }
    }
  }elsif($evalue =~ /^B([^B]+)$/){
    $evalue_small = $1;
    unless($evalue_small){
      $zero_key = 4;
    }
  }elsif($evalue =~ /^S([^S]+)$/){
    $evalue_big = $1;
    unless($evalue_big){
      print STDERR "ERROR: The minimum E_value is 0. Please change it.\n$usage\n";
      exit(1);
    }
  }elsif($evalue =~ /^same$/i){
    $evalue_same = 1;
  }else{
    print STDERR "ERROR: '$evalue' is a illegal arugument of '--evalue'. Please change it.\n$usage\n";
    exit(1);
  }
  if($evalue_small){
    my $next_tmp = 1;
    if($evalue_small =~ /^[0]+$/){
      $small_zero = 2;
      $next_tmp = 2;
    }
    if($evalue_small =~ /^[0]+\.[0]*$/){
      $small_zero = 2;
      $next_tmp = 2;
    }
    if($evalue_small =~ /^0E0$/i){
      $small_zero = 2;
      $next_tmp = 2;
    }
    if($next_tmp == 1){
      if($evalue_small > 0){

      }elsif($evalue_small < 0){
        print STDERR "ERROR: $evalue_small is less than zero. The minimum E_value is 0.\n'$evalue' is a illegal arugument of '--evalue'. Please change it.\n$usage\n";
        exit(1);
      }else{
        print STDERR "ERROR: $evalue_small is not a number.\n'$evalue' is a illegal arugument of '--evalue'. Please change it.\n$usage\n";
        exit(1);
      }
    }
  }
  if($evalue_big){
    my $next_tmp = 1;
    if($evalue_big =~ /^[0]+$/){
      $big_zero = 2;
      $next_tmp = 2;
    }
    if($evalue_big =~ /^[0]+\.[0]*$/){
      $big_zero = 2;
      $next_tmp = 2;
    }
    if($evalue_big =~ /^0E0$/i){
      $big_zero = 2;
      $next_tmp = 2;
    }
    if($next_tmp == 1){
      if($evalue_big > 0){

      }elsif($evalue_big < 0){
        print STDERR "ERROR: $evalue_big is less than zero. The minimum E_value is 0.\n'$evalue' is a illegal arugument of '--evalue'. Please change it.\n$usage\n";
        exit(1);
      }else{
        print STDERR "ERROR: $evalue_big is not a number.\n'$evalue' is a illegal arugument of '--evalue'. Please change it.\n$usage\n";
        exit(1);
      }
    }
  }
  if($big_zero == 2){
    if($small_zero == 2){
      $zero_key = 3;
    }
    if($small_zero == 1){
      print STDERR "ERROR: The minimum E_value is 0. Please change it.\n$usage\n";
      exit(1);
    }
  }
  if($big_zero == 1 && $small_zero == 2){
    if($evalue_big){
      $evalue_small = 0;
    }
    unless($evalue_big){
      $zero_key = 4;
    }
  }
  if($evalue_small && $evalue_big && $big_zero == 1 && $small_zero == 1){
    if($evalue_small > $evalue_big){
      print STDERR "ERROR: $evalue_small is big than $evalue_big.\n'$evalue' is a illegal arugument of '--evalue'. Please change it.\n$usage\n";
      exit(1);
    }
  }
}

if($db_build){
  unless($db_name){
    print STDERR "ERROR: If you want build databases for PepTidy mining, please also set '--db_name'.\n$usage\n";
    exit(1);
  }
  if($db_name){
    unless($db_name =~ /^[A-Za-z][\w]*$/){
      print STDERR "ERROR: Argument of '--db_name' should start with letter. And only letter, number and underline are vaild.\n'--db_name $db_name' Please it.\n$usage\n";
      exit(1);
    }
    my $output_db = "$output/$db_name" . '_database';
    my $cmd_mkdir = "mkdir $output_db";
    system("$cmd_mkdir") == 0 or die "failed to execute command \'$cmd_mkdir\': $!\n";

    $database_seq = "$output_db/$db_name". '_database_seq';
    $database_actsite = "$output_db/$db_name" . '_seq_actsites_all';
    $database_AMB = "$output_db/$db_name" . '_seq_actsites_AMB_sum';
    $database_length = "$output_db/$db_name" . '_seq_length';
    $database_name = "$output_db/$db_name" . '_seq_name';
    $database_uniprot = "$output_db/$db_name" . '_seq_uniport_annotation_sum';

    if(-e $database_seq || -e $database_actsite || -e $database_AMB || -e $database_length || -e $database_name || -e $database_uniprot){
      print STDERR "FATAL ERROR: Please redownload \'analysis.pl\'.\n";
      exit(1);
    }
  }
}
unless($db_build){
  if($db_name){
    print STDERR "ERROR: If you want build databases '$db_name' for PepTidy mining, please also set '--db_build'.\n$usage\n";
    exit(1);
  }
}

check_seq_format($annotation_seq);
my $count_sum = scalar(keys %query_name_hash);

check_query_format($annotation_query);
$count_sum = scalar(keys %uniprot_name_hash);

check_uniprot_format($annotation_uniprot);




sub check_uniprot_format{
  print STDOUT "Check uniprotformat ...";
  my $in = $_[0];
  my %uniprot_name_tmp_hash;
  my $start_key = 'NA';
  open IN, $in;
  while(<IN>){
    chomp;
    if(/^ID   .+/){
      if($start_key eq 'NA'){
        $start_key = 1;
      }elsif($start_key == 5){
        $start_key = 1;
      }else{
        print STDERR "ERROR: File '$input/*_annotation_results_uniprot' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
        exit(1);
      }
    }
    if(/^AC   (.+)$/){
      my $name_tmp = $1;
      my @name_array = ($name_tmp =~ /\w+/g);
      foreach my $name (@name_array){
        if($uniprot_name_hash{$name}){
          $uniprot_name_tmp_hash{$name} += 1;
        }
      }
    }
    if(/^OX   (.+)$/){
      if($start_key == 1){
        $start_key = 2;
      }else{
        print STDERR "ERROR: File '$input/*_annotation_results_uniprot' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
        exit(1);
      }
    }
    if(/^PE   (.+)$/){
      if($start_key == 2){
        $start_key = 3;
      }else{
        print STDERR "ERROR: File '$input/*_annotation_results_uniprot' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
        exit(1);
      }
    }
    if(/^SQ   (.+)$/){
      if($start_key == 3){
        $start_key = 4;
      }else{
        print STDERR "ERROR: File '$input/*_annotation_results_uniprot' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
        exit(1);
      }
    }
    if(/^\/\//){
      if($start_key == 4){
        $start_key = 5;
      }else{
        print STDERR "ERROR: File '$input/*_annotation_results_uniprot' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
        exit(1);
      }
    }
  }
  close IN;  
  
  my $count_check = scalar(keys %uniprot_name_tmp_hash);
  unless($count_sum == $count_check){
    print STDERR "ERROR: File '$input/*_annotation_results_uniprot' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
    exit(1);
  }
}

sub check_query_format{
  print STDOUT "Check query format ...\n";
  my $in = $_[0];
  my $query_name;
  my %query_name_tmp_hash;
  my $start_key = 'NA';
  my $QN = 1;
  my $AT = 1;
  my %AT_hash;
  my $DO = 1;
  my $FM = 1;
  my $MD = 1;
  my $EV = 1;
  my $SC = 1;
  my $UP = 1;
  my $DE = 1;
  my $end = 1;
  my $loop = 1;
  open IN, $in;
  while(<IN>){
    chomp;
    if(/^QN   (.+)$/){
      if($QN == 1 && $loop == 1){
        $query_name = $1;
        $query_name_tmp_hash{$query_name} += 1;
        unless($query_name_hash{$query_name}){
          print STDERR "ERROR1: File '$input/*_annotation_results_query' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
          exit(1);
        }
        $start_key = 1;
        $AT = 1;
      }
      if($QN == 1 && $loop > 1 && $end == 2 && $start_key == 10){
        $query_name = $1;
        $query_name_tmp_hash{$query_name} += 1;
        $start_key = 1;
        $AT = 1;
      }
    }
    if(/^AT   ([^\:]+): (.+)$/){
      my $AT_tag = $1;
      my $sites = $2;
      my $query_AT = $query_name . $AT_tag;
      if($query_name_tmp_hash{$query_name} == 1){
        $AT_hash{$query_AT} = $sites;
      }
      if($query_name_tmp_hash{$query_name} > 1){
        unless($AT_hash{$query_AT}){
          $bais_hash{$query_name} += 1;
        }
        if($AT_hash{$query_AT}){
          unless($sites eq $AT_hash{$query_AT}){
            $bais_hash{$query_name} += 1;
          }
        }
      }
      if($start_key == 1 && $AT == 1){
        $start_key = 2;
        $AT = 2;
        $DO = 1;
      }elsif($start_key == 2 && $AT == 2){
        
      }elsif($start_key == 9 && $AT == 2 && $DE == 2){
        $start_key = 2;
        unless($custom){
          $DO = 1;
        }
        $bais_hash{$query_name} += 1;
        if($custom){
          $MD = 1;
        }
      }else{
        print STDERR "ERROR2: File '$input/*_annotation_results_query' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
        exit(1);
      }
    }
    unless($custom){
    if(/^DO   Domain: .+/){
      if($start_key == 2 && $DO == 1){
        $start_key = 3;
        $DO = 2;
        $FM = 1;
      }else{
        print STDERR "ERROR3: File '$input/*_annotation_results_query' is incompatible.\n Please run PepTidy angain and add option '-merops' to obtain compatible annotation results.\n";
        exit(1);
      }
    }
    if(/^FM   Family: (\w).+/){
      my $family = $1;
      $query_family_hash{$family} += 1;
      if($start_key == 3 && $FM == 1){
        $start_key = 4;
        $FM = 2;
        $MD = 1;
      }else{
        print STDERR "ERROR4: File '$input/*_annotation_results_query' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
        exit(1);
      }
    }
    }
    if(/^MD   MEROPS_id: .+/){
      if($custom){
        unless($start_key == 2){
          print STDERR "ERROR5: File '$input/*_annotation_results_query' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
          exit(1);
        }
        if($start_key == 2){
          $start_key = 4;
        }
      }
      if($start_key == 4 && $MD == 1){
        $start_key = 5;
        $MD = 2;
        $EV = 1;
      }else{
        print STDERR "ERROR6: File '$input/*_annotation_results_query' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
        exit(1);
      }
    }
    if(/^EV   E-Value: (.+)/){
      my $e_value = $1;
      unless($e_value =~ /same/){
        if($evalue_max){
          if($e_value > $evalue_max){
            $evalue_max = $e_value;
          }
          if($e_value < $evalue_min){
            $evalue_min = $e_value ;
          }
        }else{
          $evalue_max = $e_value;
          $evalue_min = $e_value;
        }
      }

      if($e_value eq 'same'){
        $query_same_hash{$query_name} += 1;
        if($evalue_same){
          $query_name_evalue_hash{$query_name} += 1;
        }
      }elsif($zero_key == 3){
        if($e_value == 0){
          $query_name_evalue_hash{$query_name} += 1;
        }
      }elsif($zero_key == 4){
        $query_name_evalue_hash{$query_name} += 1;
      }else{
        if($evalue_small){
          if($e_value >= $evalue_small){
            if($evalue_big){
              if($e_value <= $evalue_big){
                $query_name_evalue_hash{$query_name} += 1;
              }
            }
            unless($evalue_big){
              $query_name_evalue_hash{$query_name} += 1;
            }
          }
        }
        unless($evalue_small){
          if($evalue_big){
            if($e_value <= $evalue_big){
              $query_name_evalue_hash{$query_name} += 1;
            }
          }
        }
      }
      if($start_key == 5 && $EV == 1){
        $start_key = 6;
        $EV = 2;
        $SC = 1;
      }else{
        print STDERR "ERROR7: File '$input/*_annotation_results_query' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
        exit(1);
      }
    }
    if(/^SC   Score: .+/){
      if($start_key == 6 && $SC == 1){
        $start_key = 7;
        $SC = 2;
        $UP = 1;
      }else{
        print STDERR "ERROR8: File '$input/*_annotation_results_query' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
        exit(1);
      }
    }
    if(/^UP   Uniport_id: (.+)$/){
      my $uniprot_id = $1;
      $uniprot_name_hash{$uniprot_id} += 1;
      if($start_key == 7 && $UP == 1){
        $start_key = 8;
        $UP = 2;
        $DE = 1;
      }else{
        print STDERR "ERROR9: File '$input/*_annotation_results_query' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
        exit(1);
      }
    }
    if(/^DE   .+/){
      if($start_key == 8 && $DE == 1){
        $start_key = 9;
        $DE = 2;
        $end = 1;
      }elsif($start_key == 9 && $DE == 2 && $end == 1){
        
      }else{
        print STDERR "ERROR10: File '$input/*_annotation_results_query' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
        exit(1);
      }
    }
    if(/^end$/){
      if($start_key == 9 && $end == 1){
        $start_key = 10;
        $QN = 1;
        $end = 2;
        $loop++;
      }else{
        print "$query_name\n"; ###add
        print STDERR "ERROR11: File '$input/*_annotation_results_query' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
        exit(1);
      }
    }
  }
  close IN;

  my $count_check = scalar(keys %query_name_tmp_hash);
  unless($count_check == $count_sum){
    print STDERR "ERROR12: File '$input/*_annotation_results_query' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
    exit(1);
  }
}

sub check_seq_format{
  print STDOUT "Check seq format ...\n";
  my $in = $_[0];
  my $count_line;
  my $count_start;
  my $count_add = 1;
  open IN, $in;
  while(<IN>){
    if(/^\n$/){
      $count_add++;
      next;
    }
    chomp;
    $count_line++;
    if($count_line == 1){
      unless(/^>/){
        print STDERR "ERROR: Error at $count_line line of $in. The first line must start with \'>\'. Please use fasta format.\n";
        exit(1);
      }
    }

    if(/^>(.+)$/){
      if($count_line > 1){
        if($count_line == ($count_start + 1)){
          my $count_start_tmp;
          $count_start_tmp = $count_start + $count_add -1;
          print STDERR "ERROR: Error at $count_start_tmp line of $in. No sequences. Please use fasta format.\n";
          exit(1);
        }
      }
      $count_start = $count_line;
      my $name_tmp = $1;
      $query_name_hash{$name_tmp} += 1;
      if($query_name_hash{$name_tmp} > 1){
        print STDERR "ERROR: query name: \'$name_tmp\' repeat in $in. Please make sure no repeat names in query file.\n";
        exit(1);
      }
    }

    unless(/^>/){
      my $tmp_seq = $_;
      my $last_char = substr($tmp_seq, -1);
      if($last_char =~ /\r/){
        chop($tmp_seq);
      }
      my $len = length $tmp_seq;
      my @ERROR_chr = ($tmp_seq =~ /[^ACDEFGHIKLMNPQRSTVWYX]/ig);
      if(@ERROR_chr){
        my $count_line_tmp;
        $count_line_tmp = $count_line + $count_add -1;
        print STDOUT "WARNING: Unexpected string \'@ERROR_chr\' at $count_line_tmp line of $in.\n";
      }
    }
  }
  close IN;
}

#check whether files exist
sub Check_peptidy_files{
  my $file = $_[0];
  unless($file){
    print STDERR "FATAL ERROR: Please redownload \'analysis.pl\'.\n";
    exit(1);
  }
  if(! -e $file || -d $file){
    print STDERR "ERROR: No file \'$file\'. Please check.\n";
    exit(1);
  }
}

print STDOUT "All input files and settings are ok.\n\n";
print STDOUT '-'x60;
print STDOUT "\n\nResults Summery\n\n";

############################ Main functon ##################################

if($statistics){
  statistics($annotation_query, 'statistics');
}else{
  statistics($annotation_query, 'NA');
}

if($query_list){
  my $out = "$output/query";
  open OUT, '>>', $out;
  my @tmp_array = sort (keys %query_name_hash);
  $count_sum = scalar @tmp_array;
  print OUT "QUERY SUM: $count_sum\n";
  foreach my $tmp (@tmp_array){
    print OUT "$tmp\n";
  }
  close OUT;
}

if($active_list || $metal_list  || $binding_list  || $active_metal_list  || $active_binding_list  || $metal_binding_list  || $active_metal_binding_list  || $EC_list  || $family_list  || $name_list || $db_build){
  read_query($annotation_query);
  if($active_list){
    if(%active_list_hash){
      my $out = "$output/actsites";
      open OUT, '>>', $out;
      my @tmp_array = sort (keys %active_list_hash);
      $count_sum = scalar @tmp_array;
      print OUT "\nActive sites SUM: $count_sum\n";
      foreach my $tmp (@tmp_array){
        print OUT "$tmp\n";
      }
      close OUT;
    }else{
      print STDOUT "No Active sites.\n";
    }
  }

  if($metal_list){
    if(%metal_list_hash){
      my $out = "$output/actsites";
      open OUT, '>>', $out;
      my @tmp_array = sort (keys %metal_list_hash);
      $count_sum = scalar @tmp_array;
      print OUT "\nMetal sites SUM: $count_sum\n";
      foreach my $tmp (@tmp_array){
        print OUT "$tmp\n";
      }
      close OUT;
    }else{
      print STDOUT "No Metal sites.\n";
    }
  }

  if($binding_list){
    if(%binding_list_hash){
      my $out = "$output/actsites";
      open OUT, '>>', $out;
      my @tmp_array = sort (keys %binding_list_hash);
      $count_sum = scalar @tmp_array;
      print OUT "\nBinding sites SUM: $count_sum\n";
      foreach my $tmp (@tmp_array){
        print OUT "$tmp\n";
      }
      close OUT;
    }else{
      print STDOUT "No Binding sites.\n";
    }
  }

  if($active_metal_list){
    if(%active_metal_list_hash){
      my $out = "$output/actsites";
      open OUT, '>>', $out;
      my @tmp_array = sort (keys %active_metal_list_hash);
      $count_sum = scalar @tmp_array;
      print OUT "\nActive and Metal sites SUM: $count_sum\n";
      foreach my $tmp (@tmp_array){
        print OUT "$tmp\n";
      }
      close OUT;
    }else{
      print STDOUT "No Active and Metal sites.\n";
    }
  }

  if($active_binding_list){
    if(%active_binding_list_hash){
      my $out = "$output/actsites";
      open OUT, '>>', $out;
      my @tmp_array = sort (keys %active_binding_list_hash);
      $count_sum = scalar @tmp_array;
      print OUT "\nActive and Binding sites SUM: $count_sum\n";
      foreach my $tmp (@tmp_array){
        print OUT "$tmp\n";
      }
      close OUT;
    }else{
      print STDOUT "No Active and Binding sites.\n";
    }
  }

  if($metal_binding_list){
    if(%metal_binding_list_hash){
      my $out = "$output/actsites";
      open OUT, '>>', $out;
      my @tmp_array = sort (keys %metal_binding_list_hash);
      $count_sum = scalar @tmp_array;
      print OUT "\nMetal and Binding sites SUM: $count_sum\n";
      foreach my $tmp (@tmp_array){
        print OUT "$tmp\n";
      }
      close OUT;
    }else{
      print STDOUT "No Metal and Binding sites.\n";
    }
  }

  if($active_metal_binding_list){
    if(%active_metal_binding_list_hash){
      my $out = "$output/actsites";
      open OUT, '>>', $out;
      my @tmp_array = sort (keys %active_metal_binding_list_hash);
      $count_sum = scalar @tmp_array;
      print OUT "\nActive and Metal and Binding sites SUM: $count_sum\n";
      foreach my $tmp (@tmp_array){
        print OUT "$tmp\n";
      }
      close OUT;
    }else{
      print STDOUT "No Active and Metal and Binding sites.\n";
    }
  }

  if($EC_list){
    if(%EC_list_hash){
      my $out = "$output/EC_number";
      open OUT, '>>', $out;
      my @tmp_array = sort (keys %EC_list_hash);
      $count_sum = scalar @tmp_array;
      print OUT "E.C. SUM: $count_sum\n";
      foreach my $tmp (@tmp_array){
        print OUT "$tmp\n";
      }
      close OUT;
    }else{
      print STDOUT "No E.C.\n";
    }
  }

  unless($custom){
    if($family_list){
      if(%family_list_hash){
        my $out = "$output/family";
        open OUT, '>>', $out;
        my @tmp_array = sort (keys %family_list_hash);
        $count_sum = scalar @tmp_array;
        print OUT "FAMILY SUM: $count_sum\n";
        foreach my $tmp (@tmp_array){
          print OUT "$tmp\n";
        }
        close OUT;
      }else{
        print STDOUT "No Family.\n";
      }
    }
  }

  if($name_list){
    if(%name_list_hash){
      my $out = "$output/name";
      open OUT, '>>', $out;
      my @tmp_array = sort (keys %name_list_hash);
      $count_sum = scalar @tmp_array;
      print OUT "NAME SUM: $count_sum\n";
      foreach my $tmp (@tmp_array){
        print OUT "$tmp\n";
      }
      close OUT;
    }else{
      print STDOUT "No Name.\n";
    }
  }

}else{
  if($ALL || $reference || @annotation_array || $note || $active_site || $metal_site || $binding_site || $annotation_list){
    read_query($annotation_query);
    read_uniprot($annotation_uniprot);
    if($annotation_list){
      if(%annotation_list_hash){
        my $out = "$output/annotation";
        open OUT, '>>', $out;
        my @tmp_array = sort (keys %annotation_list_hash);
        $count_sum = scalar @tmp_array;
        print OUT "NAME SUM: $count_sum\n";
        foreach my $tmp (@tmp_array){
          print OUT "$tmp\n";
        }
        close OUT;
      }else{
        print STDOUT "No Annotation.\n";
      }
    }
  }
}

my %query_search_include_name_hash;
my %query_search_all_hash;
my %query_search_exclude_name_hash;
my $tmp_all;

#$next_switch

search_query($annotation_query);
if(@de_query_array){
  foreach my $tmp (@de_query_array){
    $exclude_hash{"de_query: $tmp"} += 1;
    if($query_name_hash{$tmp}){
      $query_search_exclude_name_hash{$tmp} += 1;
    }else{
      foreach my $query_name (keys %query_name_search_all_hash){
        my $name_tmp_all = $query_name;
        if($name_tmp_all =~ /$tmp/i){
          $query_search_exclude_name_hash{$query_name} += 1;
        }
      }
    }
  }
}

if(@query_array){
  foreach my $tmp (@query_array){
    if($query_name_hash{$tmp}){
      $query_hash{$tmp} += 1;
    }else{
      foreach my $query_name (keys %query_name_search_all_hash){
        my $name_tmp_all = $query_name;
        if($name_tmp_all =~ /$tmp/i){
          $query_hash{$query_name} += 1;
        }
      }
    }
  }
}

foreach my $query_name (keys %query_search_all_hash){
  $tmp_all = $query_search_all_hash{$query_name};
  if(@de_EC_array){
    foreach my $tmp (@de_EC_array){
      $exclude_hash{"de_EC: $tmp"} += 1;
      my $tmp_search = '#E#' . $tmp . '_';
      if($tmp_all =~ /$tmp_search/i){
        $query_search_exclude_name_hash{$query_name} += 1;
        next;
      }
    }
  }
  unless($custom){
    if(@de_family_array){
      foreach my $tmp (@de_family_array){
        if($tmp =~ /^[ASMCGPUNT]00$/){
          $tmp =~ s/0/\\d/g;
        }
        if($tmp =~ /^([ASMCGPUNT])([1-9])$/){
          $tmp = $1 . '0' . $2;
        }
        if($tmp =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $tmp = $1 . $2;
        }
        $exclude_hash{"de_family: $tmp"} += 1;
        my $tmp_search = '#F#' . $tmp . '_';
        if($tmp_all =~ /$tmp_search/i){
          $query_search_exclude_name_hash{$query_name} += 1;
          next;
        }
      }
    }
  }
  if(@de_name_array){
    foreach my $tmp (@de_name_array){
      $exclude_hash{"de_name: $tmp"} += 1;
      if($tmp_all =~ /#N#[^#]*$tmp/i){
        $query_search_exclude_name_hash{$query_name} += 1;
        next;
      }
    }
  }
  if(@de_active_array){
    foreach my $tmp (@de_active_array){
      $exclude_hash{"de_active: $tmp"} += 1;
      my $tmp_search = '#A#' . $tmp . '_';
      if($tmp_all =~ /$tmp_search/i){
        $query_search_exclude_name_hash{$query_name} += 1;
        next;
      }
    }
  }
  if(@de_metal_array){
    foreach my $tmp (@de_metal_array){
      $exclude_hash{"de_metal: $tmp"} += 1;
      my $tmp_search = '#M#' . $tmp . '_';
      if($tmp_all =~ /$tmp_search/i){
        $query_search_exclude_name_hash{$query_name} += 1;
        next;
      }
    }
  }
  if(@de_binding_array){
    foreach my $tmp (@de_binding_array){
      $exclude_hash{"de_binding: $tmp"} += 1;
      my $tmp_search = '#B#' . $tmp . '_';
      if($tmp_all =~ /$tmp_search/i){
        $query_search_exclude_name_hash{$query_name} += 1;
        next;
      }
    }
  }
  if(@de_active_metal_array){
    $next_switch =1;
    foreach my $tmp (@de_active_metal_array){
      $exclude_hash{"de_active_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $active = $1;
        my $metal = $2;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_active_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_active_binding_array){
      $exclude_hash{"de_active_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $active = $1;
        my $binding = $2;
        match_active($active);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_metal_binding_array){
      $exclude_hash{"de_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $metal = $1;
        my $binding = $2;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_active_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_active_metal_binding_array){
      $exclude_hash{"de_active_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $active = $1;
        my $metal = $2;
        my $binding = $3;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal); 
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    } 
  }
  if(@de_EC_name_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_name_array){
      $exclude_hash{"de_EC_name: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $name = $2;
        match_EC($EC);
        next if $next_switch == 1;
        match_name($name);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_active_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_active_array){
      $exclude_hash{"de_EC_active: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $active = $2;
        match_EC($EC);
        next if $next_switch == 1;
        match_active($active);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_metal_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_metal_array){
      $exclude_hash{"de_EC_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $metal = $2;
        match_EC($EC);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_binding_array){
      $exclude_hash{"de_EC_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $binding = $2;
        match_EC($EC);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_active_metal_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_active_metal_array){
      $exclude_hash{"de_EC_active_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $active = $2;
        my $metal = $3;
        match_EC($EC);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_active_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_active_binding_array){
      $exclude_hash{"de_EC_active_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $active = $2;
        my $binding = $3;
        match_EC($EC);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_metal_binding_array){
      $exclude_hash{"de_EC_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $metal = $2;
        my $binding = $3;
        match_EC($EC);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_active_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_active_metal_binding_array){
      $exclude_hash{"de_EC_active_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1; 
        my $active = $2;
        my $metal = $3;
        my $binding = $4;
        match_EC($EC);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal); 
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  unless($custom){
  if(@de_family_name_array){
    $next_switch =1;
    foreach my $tmp (@de_family_name_array){
      $exclude_hash{"de_family_name: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $name = $2;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }  
  if(@de_family_active_array){
    $next_switch =1;
    foreach my $tmp (@de_family_active_array){
      $exclude_hash{"de_family_active: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $active = $2;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_active($active);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_family_metal_array){
    $next_switch =1;
    foreach my $tmp (@de_family_metal_array){
      $exclude_hash{"de_family_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $metal = $2;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_family_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_family_binding_array){
      $exclude_hash{"de_family_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $binding = $2;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_family_active_metal_array){
    $next_switch =1;
    foreach my $tmp (@de_family_active_metal_array){
      $exclude_hash{"de_family_active_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $active = $2;
        my $metal = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_family_active_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_family_active_binding_array){
      $exclude_hash{"de_family_active_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $active = $2;
        my $binding = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_family_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_family_metal_binding_array){
      $exclude_hash{"de_family_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $metal = $2;
        my $binding = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_family_active_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_family_active_metal_binding_array){
      $exclude_hash{"de_family_active_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $active = $2;
        my $metal = $3;
        my $binding = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  }
  if(@de_name_active_array){
    $next_switch =1;
    foreach my $tmp (@de_name_active_array){
      $exclude_hash{"de_name_active: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $name = $1;
        my $active = $2;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_name_metal_array){
    $next_switch =1;
    foreach my $tmp (@de_name_metal_array){
      $exclude_hash{"de_name_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $name = $1;
        my $metal = $2;
        match_name($name);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_name_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_name_binding_array){
      $exclude_hash{"de_name_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $name = $1;
        my $binding = $2;
        match_name($name);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_name_active_metal_array){
    $next_switch =1;
    foreach my $tmp (@de_name_active_metal_array){
      $exclude_hash{"de_name_active_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $name = $1;
        my $active = $2;
        my $metal = $3;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_name_active_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_name_active_binding_array){
      $exclude_hash{"de_name_active_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $name = $1;
        my $active = $2;
        my $binding = $3;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_name_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_name_metal_binding_array){
      $exclude_hash{"de_name_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $name = $1;
        my $metal = $2;
        my $binding = $3;
        match_name($name);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_name_active_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_name_active_metal_binding_array){
      $exclude_hash{"de_name_active_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $name = $1;
        my $active = $2;
        my $metal = $3;
        my $binding = $4;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  unless($custom){
  if(@de_EC_family_active_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_family_active_array){
      $exclude_hash{"de_EC_family_active: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $active = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_active($active);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_family_metal_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_family_metal_array){
      $exclude_hash{"de_EC_family_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $metal = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_family_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_family_binding_array){
      $exclude_hash{"de_EC_family_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $binding = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_family_active_metal_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_family_active_metal_array){
      $exclude_hash{"de_EC_family_active_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $active = $3;
        my $metal = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_active($active);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_family_active_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_family_active_binding_array){
      $exclude_hash{"de_EC_family_active_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $active = $3;
        my $binding = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_family_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_family_metal_binding_array){
      $exclude_hash{"de_EC_family_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $metal = $3;
        my $binding = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_family_active_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_family_active_metal_binding_array){
      $exclude_hash{"de_EC_family_active_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $active = $3;
        my $metal = $4;
        my $binding = $5;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal); 
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  }
  if(@de_EC_name_active_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_name_active_array){
      $exclude_hash{"de_EC_name_active: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $name = $2;
        my $active = $3;
        match_EC($EC);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_name_metal_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_name_metal_array){
      $exclude_hash{"de_EC_name_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $name = $2;
        my $metal = $3;
        match_EC($EC);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_name_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_name_binding_array){
      $exclude_hash{"de_EC_name_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $name = $2;
        my $binding = $3;
        match_EC($EC);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_name_active_metal_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_name_active_metal_array){
      $exclude_hash{"de_EC_name_active_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $name = $2;
        my $active = $3;
        my $metal = $4;
        match_EC($EC);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_name_active_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_name_active_binding_array){
      $exclude_hash{"de_EC_name_active_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $name = $2;
        my $active = $3;
        my $binding = $4;
        match_EC($EC);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_name_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_name_metal_binding_array){
      $exclude_hash{"de_EC_name_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $name = $2;
        my $metal = $3;
        my $binding = $4;
        match_EC($EC);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_name_active_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_name_active_metal_binding_array){
      $exclude_hash{"de_EC_name_active_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $name = $2;
        my $active = $3;
        my $metal = $4;
        my $binding = $5;
        match_EC($EC);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  unless($custom){
  if(@de_family_name_active_array){
    $next_switch =1;
    foreach my $tmp (@de_family_name_active_array){
      $exclude_hash{"de_family_name_active: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $name = $2;
        my $active = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_family_name_metal_array){
    $next_switch =1;
    foreach my $tmp (@de_family_name_metal_array){
      $exclude_hash{"de_family_name_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $name = $2;
        my $metal = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_family_name_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_family_name_binding_array){
      $exclude_hash{"de_family_name_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $name = $2;
        my $binding = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_family_name_active_metal_array){
    $next_switch =1;
    foreach my $tmp (@de_family_name_active_metal_array){
      $exclude_hash{"de_family_name_active_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $name = $2;
        my $active = $3;
        my $metal = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_family_name_active_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_family_name_active_binding_array){
      $exclude_hash{"de_family_name_active_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $name = $2;
        my $active = $3;
        my $binding = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_family_name_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_family_name_metal_binding_array){
      $exclude_hash{"de_family_name_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $name = $2;
        my $metal = $3;
        my $binding = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_family_name_active_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_family_name_active_metal_binding_array){
      $exclude_hash{"de_family_name_active_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $name = $2;
        my $active = $3;
        my $metal = $4;
        my $binding = $5;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_family_name_active_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_family_name_active_array){
      $exclude_hash{"de_EC_family_name_active: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $name = $3;
        my $active = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_family_name_metal_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_family_name_metal_array){
      $exclude_hash{"de_EC_family_name_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1; 
        my $family = $2;
        my $name = $3;
        my $metal = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_family_name_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_family_name_binding_array){
      $exclude_hash{"de_EC_family_name_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $name = $3;
        my $binding = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_family_name_active_metal_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_family_name_active_metal_array){
      $exclude_hash{"de_EC_family_name_active_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $name = $3;
        my $active = $4;
        my $metal = $5;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_family_name_active_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_family_name_active_binding_array){
      $exclude_hash{"de_EC_family_name_active_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $name = $3;
        my $active = $4;
        my $binding = $5;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_family_name_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_family_name_metal_binding_array){
      $exclude_hash{"de_EC_family_name_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $name = $3;
        my $metal = $4;
        my $binding = $5;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@de_EC_family_name_active_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@de_EC_family_name_active_metal_binding_array){
      $exclude_hash{"de_EC_family_name_active_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $name = $3;
        my $active = $4;
        my $metal = $5;
        my $binding = $6;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal); 
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_exclude_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  }

##############match#############

  if(@EC_array){
    foreach my $tmp (@EC_array){
      $include_hash{"EC: $tmp"} += 1;
      my $tmp_search = '#E#' . $tmp . '_';
      if($tmp_all =~ /$tmp_search/i){
        $query_search_include_name_hash{$query_name} += 1;
        next;
      }
    }
  }
  unless($custom){
  if(@family_array){
    foreach my $tmp (@family_array){
      if($tmp =~ /^[ASMCGPUNT]00$/){
        $tmp =~ s/0/\\d/g;
      }
      if($tmp =~ /^([ASMCGPUNT])([1-9])$/){
        $tmp = $1 . '0' . $2;
      }
      if($tmp =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
        $tmp = $1 . $2;
      }
      $include_hash{"family: $tmp"} += 1;
      my $tmp_search = '#F#' . $tmp . '_';
      if($tmp_all =~ /$tmp_search/i){
        $query_search_include_name_hash{$query_name} += 1;
        next;
      }
    }
  }
  }
  if(@name_array){
    foreach my $tmp (@name_array){
      $include_hash{"name: $tmp"} += 1;
      if($tmp_all =~ /#N#[^#]*$tmp/i){
        $query_search_include_name_hash{$query_name} += 1;
        next;
      }
    }
  }
  if(@active_array){
    foreach my $tmp (@active_array){
      $include_hash{"active: $tmp"} += 1;
      my $tmp_search = '#A#' . $tmp . '_';
      if($tmp_all =~ /$tmp_search/i){
        $query_search_include_name_hash{$query_name} += 1;
        next;
      }
    }
  }
  if(@metal_array){
    foreach my $tmp (@metal_array){
      $include_hash{"metal: $tmp"} += 1;
      my $tmp_search = '#M#' . $tmp . '_';
      if($tmp_all =~ /$tmp_search/i){
        $query_search_include_name_hash{$query_name} += 1;
        next;
      }
    }
  }
  if(@binding_array){
    foreach my $tmp (@binding_array){
      $include_hash{"binding: $tmp"} += 1;
      my $tmp_search = '#B#' . $tmp . '_';
      if($tmp_all =~ /$tmp_search/i){
        $query_search_include_name_hash{$query_name} += 1;
        next;
      }
    }
  }
  if(@active_metal_array){
    $next_switch =1;
    foreach my $tmp (@active_metal_array){
      $include_hash{"active_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $active = $1;
        my $metal = $2;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@active_binding_array){
    $next_switch =1;
    foreach my $tmp (@active_binding_array){
      $include_hash{"active_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $active = $1;
        my $binding = $2;
        match_active($active);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@metal_binding_array){
      $include_hash{"metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $metal = $1;
        my $binding = $2;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@active_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@active_metal_binding_array){
      $include_hash{"active_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $active = $1;
        my $metal = $2;
        my $binding = $3;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal); 
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    } 
  }
  if(@EC_name_array){
    $next_switch =1;
    foreach my $tmp (@EC_name_array){
      $include_hash{"EC_name: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $name = $2;
        match_EC($EC);
        next if $next_switch == 1;
        match_name($name);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_active_array){
    $next_switch =1;
    foreach my $tmp (@EC_active_array){
      $include_hash{"EC_active: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $active = $2;
        match_EC($EC);
        next if $next_switch == 1;
        match_active($active);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_metal_array){
    $next_switch =1;
    foreach my $tmp (@EC_metal_array){
      $include_hash{"EC_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $metal = $2;
        match_EC($EC);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_binding_array){
    $next_switch =1;
    foreach my $tmp (@EC_binding_array){
      $include_hash{"EC_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $binding = $2;
        match_EC($EC);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_active_metal_array){
    $next_switch =1;
    foreach my $tmp (@EC_active_metal_array){
      $include_hash{"EC_active_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $active = $2;
        my $metal = $3;
        match_EC($EC);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_active_binding_array){
    $next_switch =1;
    foreach my $tmp (@EC_active_binding_array){
      $include_hash{"EC_active_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $active = $2;
        my $binding = $3;
        match_EC($EC);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@EC_metal_binding_array){
      $include_hash{"EC_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $metal = $2;
        my $binding = $3;
        match_EC($EC);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_active_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@EC_active_metal_binding_array){
      $include_hash{"EC_active_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1; 
        my $active = $2;
        my $metal = $3;
        my $binding = $4;
        match_EC($EC);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal); 
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  unless($custom){
  if(@family_name_array){
    $next_switch =1;
    foreach my $tmp (@family_name_array){
      $include_hash{"family_name: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $name = $2;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@family_active_array){
    $next_switch =1;
    foreach my $tmp (@family_active_array){
      $include_hash{"family_active: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $active = $2;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_active($active);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@family_metal_array){
    $next_switch =1;
    foreach my $tmp (@family_metal_array){
      $include_hash{"family_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $metal = $2;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@family_binding_array){
    $next_switch =1;
    foreach my $tmp (@family_binding_array){
      $include_hash{"family_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $binding = $2;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@family_active_metal_array){
    $next_switch =1;
    foreach my $tmp (@family_active_metal_array){
      $include_hash{"family_active_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $active = $2;
        my $metal = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@family_active_binding_array){
    $next_switch =1;
    foreach my $tmp (@family_active_binding_array){
      $include_hash{"family_active_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $active = $2;
        my $binding = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@family_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@family_metal_binding_array){
      $include_hash{"family_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $metal = $2;
        my $binding = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@family_active_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@family_active_metal_binding_array){
      $include_hash{"family_active_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $active = $2;
        my $metal = $3;
        my $binding = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  }
  if(@name_active_array){
    $next_switch =1;
    foreach my $tmp (@name_active_array){
      $include_hash{"name_active: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $name = $1;
        my $active = $2;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@name_metal_array){
    $next_switch =1;
    foreach my $tmp (@name_metal_array){
      $include_hash{"name_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $name = $1;
        my $metal = $2;
        match_name($name);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@name_binding_array){
    $next_switch =1;
    foreach my $tmp (@name_binding_array){
      $include_hash{"name_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)/){
        my $name = $1;
        my $binding = $2;
        match_name($name);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@name_active_metal_array){
    $next_switch =1;
    foreach my $tmp (@name_active_metal_array){
      $include_hash{"name_active_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $name = $1;
        my $active = $2;
        my $metal = $3;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@name_active_binding_array){
    $next_switch =1;
    foreach my $tmp (@name_active_binding_array){
      $include_hash{"name_active_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $name = $1;
        my $active = $2;
        my $binding = $3;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@name_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@name_metal_binding_array){
      $include_hash{"name_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $name = $1;
        my $metal = $2;
        my $binding = $3;
        match_name($name);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@name_active_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@name_active_metal_binding_array){
      $include_hash{"name_active_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $name = $1;
        my $active = $2;
        my $metal = $3;
        my $binding = $4;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  unless($custom){
  if(@EC_family_active_array){
    $next_switch =1;
    foreach my $tmp (@EC_family_active_array){
      $include_hash{"EC_family_active: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $active = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_active($active);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_family_metal_array){
    $next_switch =1;
    foreach my $tmp (@EC_family_metal_array){
      $include_hash{"EC_family_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $metal = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_family_binding_array){
    $next_switch =1;
    foreach my $tmp (@EC_family_binding_array){
      $include_hash{"EC_family_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $binding = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_family_active_metal_array){
    $next_switch =1;
    foreach my $tmp (@EC_family_active_metal_array){
      $include_hash{"EC_family_active_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $active = $3;
        my $metal = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_active($active);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_family_active_binding_array){
    $next_switch =1;
    foreach my $tmp (@EC_family_active_binding_array){
      $include_hash{"EC_family_active_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $active = $3;
        my $binding = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_family_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@EC_family_metal_binding_array){
      $include_hash{"EC_family_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $metal = $3;
        my $binding = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_family_active_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@EC_family_active_metal_binding_array){
      $include_hash{"EC_family_active_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $active = $3;
        my $metal = $4;
        my $binding = $5;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal); 
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  }
  if(@EC_name_active_array){
    $next_switch =1;
    foreach my $tmp (@EC_name_active_array){
      $include_hash{"EC_name_active: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $name = $2;
        my $active = $3;
        match_EC($EC);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_name_metal_array){
    $next_switch =1;
    foreach my $tmp (@EC_name_metal_array){
      $include_hash{"EC_name_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $name = $2;
        my $metal = $3;
        match_EC($EC);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_name_binding_array){
    $next_switch =1;
    foreach my $tmp (@EC_name_binding_array){
      $include_hash{"EC_name_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $name = $2;
        my $binding = $3;
        match_EC($EC);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_name_active_metal_array){
    $next_switch =1;
    foreach my $tmp (@EC_name_active_metal_array){
      $include_hash{"EC_name_active_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $name = $2;
        my $active = $3;
        my $metal = $4;
        match_EC($EC);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_name_active_binding_array){
    $next_switch =1;
    foreach my $tmp (@EC_name_active_binding_array){
      $include_hash{"EC_name_active_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $name = $2;
        my $active = $3;
        my $binding = $4;
        match_EC($EC);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_name_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@EC_name_metal_binding_array){
      $include_hash{"EC_name_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $name = $2;
        my $metal = $3;
        my $binding = $4;
        match_EC($EC);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_name_active_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@EC_name_active_metal_binding_array){
      $include_hash{"EC_name_active_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $name = $2;
        my $active = $3;
        my $metal = $4;
        my $binding = $5;
        match_EC($EC);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  unless($custom){
  if(@family_name_active_array){
    $next_switch =1;
    foreach my $tmp (@family_name_active_array){
      $include_hash{"family_name_active: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $name = $2;
        my $active = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@family_name_metal_array){
    $next_switch =1;
    foreach my $tmp (@family_name_metal_array){
      $include_hash{"family_name_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $name = $2;
        my $metal = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@family_name_binding_array){
    $next_switch =1;
    foreach my $tmp (@family_name_binding_array){
      $include_hash{"family_name_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $name = $2;
        my $binding = $3;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@family_name_active_metal_array){
    $next_switch =1;
    foreach my $tmp (@family_name_active_metal_array){
      $include_hash{"family_name_active_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $name = $2;
        my $active = $3;
        my $metal = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@family_name_active_binding_array){
    $next_switch =1;
    foreach my $tmp (@family_name_active_binding_array){
      $include_hash{"family_name_active_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $name = $2;
        my $active = $3;
        my $binding = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@family_name_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@family_name_metal_binding_array){
      $include_hash{"family_name_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $name = $2;
        my $metal = $3;
        my $binding = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@family_name_active_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@family_name_active_metal_binding_array){
      $include_hash{"family_name_active_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $family = $1;
        my $name = $2;
        my $active = $3;
        my $metal = $4;
        my $binding = $5;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_family_name_active_array){
    $next_switch =1;
    foreach my $tmp (@EC_family_name_active_array){
      $include_hash{"EC_family_name_active: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $name = $3;
        my $active = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_family_name_metal_array){
    $next_switch =1;
    foreach my $tmp (@EC_family_name_metal_array){
      $include_hash{"EC_family_name_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1; 
        my $family = $2;
        my $name = $3;
        my $metal = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_family_name_binding_array){
    $next_switch =1;
    foreach my $tmp (@EC_family_name_binding_array){
      $include_hash{"EC_family_name_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $name = $3;
        my $binding = $4;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_family_name_active_metal_array){
    $next_switch =1;
    foreach my $tmp (@EC_family_name_active_metal_array){
      $include_hash{"EC_family_name_active_metal: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $name = $3;
        my $active = $4;
        my $metal = $5;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_family_name_active_binding_array){
    $next_switch =1;
    foreach my $tmp (@EC_family_name_active_binding_array){
      $include_hash{"EC_family_name_active_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $name = $3;
        my $active = $4;
        my $binding = $5;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_family_name_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@EC_family_name_metal_binding_array){
      $include_hash{"EC_family_name_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $name = $3;
        my $metal = $4;
        my $binding = $5;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_metal($metal);
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  if(@EC_family_name_active_metal_binding_array){
    $next_switch =1;
    foreach my $tmp (@EC_family_name_active_metal_binding_array){
      $include_hash{"EC_family_name_active_metal_binding: $tmp"} += 1;
      if($tmp =~ /^([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)_([^\_]+)/){
        my $EC = $1;
        my $family = $2;
        my $name = $3;
        my $active = $4;
        my $metal = $5;
        my $binding = $6;
        if($family =~ /^[ASMCGPUNT]00$/){
          $family =~ s/0/\\d/g;
        }
        if($family =~ /^([ASMCGPUNT])([1-9])$/){
          $family = $1 . '0' . $2;
        }
        if($family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
          $family = $1 . $2;
        }
        match_EC($EC);
        next if $next_switch == 1;
        match_family($family);
        next if $next_switch == 1;
        match_name($name);
        next if $next_switch == 1;
        match_active($active);
        next if $next_switch == 1;
        match_metal($metal); 
        next if $next_switch == 1;
        match_binding($binding);
      }
      next if $next_switch == 1;
      $query_search_include_name_hash{$query_name} += 1 if $next_switch == 2;
      next if $next_switch == 2;
    }
  }
  }

}


my $print_last = 1;
my $output_filter = "$output/filter";
unless(%exclude_hash){
  unless(%include_hash){
    unless($evalue){
      if(@query_array){
        if(%query_hash){
          open OUT, '>>', $output_filter;
          $count_sum = scalar(keys %query_hash);
          print OUT "'-query' with '@query_array'\n\n";
          print OUT "Matches SUM: $count_sum\n\n";
          print OUT "-----------------------------------------------------------------\n\n";
          foreach my $query_name (keys %query_hash){
            write_query($query_name);
            $db_query_name_hash{$query_name} += 1;
          }
          if($db_build && %db_query_name_hash){
            database_build();
          }
          close OUT;
          write_output_file();
          print STDOUT "'--query' result file: '$output_filter'.\n";
        }else{
          write_output_file();
          print STDOUT "After excuting '-query' with '@query_array', no annotated query left.\n";
        }
        exit(0);
      }
    }

    if($evalue){
      unless(%query_name_evalue_hash){
        write_output_file();
        print STDOUT "After excuting '--evalue $evalue', no annotated query left.\nPlease set the argument of '--evalue' in scope from $evalue_min to $evalue_max.\n";
        exit(0);
      }
      open OUT, '>>', $output_filter;
      $count_sum = scalar(keys %query_name_evalue_hash);
      print OUT "evalue: $evalue\nevalue matching SUM: $count_sum\n\n";
      print OUT "'-query' with '@query_array'\n\n" if @query_array;
      print OUT "-----------------------------------------------------------------\n\n";
      $count_sum = 0;
      foreach my $query_name (keys %query_name_evalue_hash){
        if(@query_array){
          next unless $query_hash{$query_name};
        }
        write_query($query_name);
        $db_query_name_hash{$query_name} += 1;
        $count_sum++;
        $print_last = 2;
      }

      if($db_build && %db_query_name_hash){
        database_build();
      }

      if($print_last == 1){
        write_output_file();
        print STDOUT "No annotated query left.\nNo matchs when set '-query' with '@query_array'.\n";
        print OUT "No annotated query left.\nNo matchs when set '-query' with '@query_array'.\n";
        close OUT;
        exit(0);
      }

      print OUT "\nfinal matching SUM: $count_sum\n";
      close OUT;

      write_output_file();
      print STDOUT "'--evalue' result file: '$output_filter'.\n";
      exit(0);
    }
    exit(0);
  }
  if(%include_hash){
    unless(%query_search_include_name_hash){
      write_output_file();
      print STDOUT "After excuting:\n";
      foreach my $criteria (sort(keys %include_hash)){
        print STDOUT "  --$criteria\n";
      }
      print STDOUT "No annotated query left.\n";
      exit(0);
    }
    if($evalue){
      unless(%query_name_evalue_hash){
        write_output_file();
        print STDOUT "After excuting '--evalue $evalue', no annotated query left.\nPlease set the argument of '--evalue' in scope of $evalue_min-$evalue_max.\n";
        exit(0);
      }
      open OUT, '>>', $output_filter;
      $count_sum = scalar(keys %query_name_evalue_hash);
      print OUT "evalue: $evalue\nevalue matching SUM: $count_sum\n\n";
      foreach my $criteria (sort(keys %include_hash)){
        print OUT "$criteria\n";
      }
      $count_sum = scalar(keys %query_search_include_name_hash);
      print OUT "search_criteria matching SUM: $count_sum\n\n";
      print OUT "'-query' with @query_array\n\n" if @query_array;
      print OUT "-----------------------------------------------------------------\n\n";
      $count_sum = 0;
      foreach my $query_name (keys %query_name_evalue_hash){
        if($query_search_include_name_hash{$query_name}){
          if(@query_array){
            next unless $query_hash{$query_name};
          }
          write_query($query_name);
          $count_sum++;
          $print_last = 2;
          $db_query_name_hash{$query_name} += 1;
        }
      }

      if($db_build && %db_query_name_hash){
        database_build();
      }

      if($print_last == 1){
        write_output_file();
        print STDOUT "No annotated query left.\n";
        print OUT "No annotated query left.\n";
        close OUT;
        exit(0);
      }
      print OUT "\nfinal matching SUM: $count_sum\n";
      close OUT;
      write_output_file();
      print STDOUT "'--evalue' result file: '$output_filter'.\n";
      exit(0);
    }
    open OUT, '>>', $output_filter;
    print OUT "After excuting:\n";
    foreach my $criteria (sort(keys %include_hash)){
      print OUT "  --$criteria\n";
    }
    $count_sum = scalar(keys %query_search_include_name_hash);
    print OUT "Search matching SUM: $count_sum\n\n";
    print OUT "'-query' with '@query_array'\n\n" if @query_array;
    print OUT "-----------------------------------------------------------------\n\n";
    $count_sum = 0;
    foreach my $query_name (keys %query_search_include_name_hash){
      if(@query_array){     
        next unless $query_hash{$query_name};
      }
      write_query($query_name);
      $db_query_name_hash{$query_name} += 1;
      $count_sum++; 
      $print_last = 2;
    }

    if($db_build && %db_query_name_hash){
      database_build();
    }

    if($print_last == 1){
      write_output_file();
      print STDOUT "No annotated query left.\nNo matchs when set '-query' with '@query_array'.\n";
      print OUT "No annotated query left.\nNo matchs when set '-query' with '@query_array'.\n";
      close OUT;
      exit(0);
    }

    print OUT "\nfinal matching SUM: $count_sum\n";
    close OUT;

    write_output_file();
    print STDOUT "search_criteria matching result file: '$output_filter'.\n";
    exit(0);
  }
  exit(0);
}

if(%exclude_hash){
  my $count_ex;
  my $count_left;
  my $count_all = scalar(keys %query_name_hash);
  if(%query_search_exclude_name_hash){
    $count_ex = scalar(keys %query_search_exclude_name_hash);
    if($count_ex == $count_all){
      write_output_file();
      print STDOUT "After excuting:\n";
      foreach my $criteria (sort(keys %include_hash)){
        print STDOUT "  --$criteria\n";
      }
      print STDOUT "No annotated query left.\n";
      exit(0);
    }
    if($count_ex > $count_all){
      print STDERR "FATAL ERROR: Please redownload \'analysis.pl\'.\n"; 
      exit(1);
    }
    $count_left = $count_all - $count_ex;
  }
  unless(%query_search_exclude_name_hash){
    $count_left = $count_all;
  }
  unless(%include_hash){
    if($evalue){
      unless(%query_name_evalue_hash){
        write_output_file();
        print STDOUT "After excuting '--evalue $evalue', no annotated query left.\nPlease set the argument of '--evalue' in scope of $evalue_min-$evalue_max.\n";
        exit(0);
      }
      open OUT, '>>', $output_filter;
      $count_sum = scalar(keys %query_name_evalue_hash);
      print OUT "evalue: $evalue\nevalue matching SUM: $count_sum\n\n";
      foreach my $criteria (sort(keys %exclude_hash)){
        print OUT "$criteria\n";
      }
      print OUT "de_search_criteria matching left SUM: $count_left\n\n";
      print OUT "'-query' with '@query_array'\n\n" if @query_array;
      print OUT "-----------------------------------------------------------------\n\n";
      $count_sum = 0;
      foreach my $query_name (keys %query_name_evalue_hash){
        unless($query_search_exclude_name_hash{$query_name}){
          if(@query_array){
            next unless $query_hash{$query_name};
          }
          write_query($query_name);
          $count_sum++;
          $print_last = 2;
          $db_query_name_hash{$query_name} += 1;
        }
      }

      if($db_build && %db_query_name_hash){
        database_build();
      }

      if($print_last == 1){
        write_output_file();
        print STDOUT "No annotated query left.\n";
        print OUT "No annotated query left.\n";
        close OUT;
        exit(0);
      }
      print OUT "\nfinal matching SUM: $count_sum\n";
      write_output_file();
      print STDOUT "'--evalue' result file: '$output_filter'.\n";
      close OUT;
      exit(0);
    }
    unless($evalue){
      open OUT, '>>', $output_filter;
      foreach my $criteria (sort(keys %exclude_hash)){
        print OUT "$criteria\n";
      }
      print OUT "\nde_search_criteria matching left SUM: $count_left\n";
      print OUT "'-query' with '@query_array'\n\n" if @query_array;
      print OUT "-----------------------------------------------------------------\n\n";
      $count_sum = 0;
      foreach my $query_name (keys %query_name_hash){
        unless($query_search_exclude_name_hash{$query_name}){
          if(@query_array){
            next unless $query_hash{$query_name};
          }
          write_query($query_name);
          $db_query_name_hash{$query_name} += 1;
          $count_sum++;
          $print_last = 2;
        }
      }

      if($db_build && %db_query_name_hash){
        database_build();
      }

      if($print_last == 1){
        write_output_file();
        print STDOUT "No annotated query left.\nNo matchs when set '-query' with '@query_array'.\n";
        print OUT "No annotated query left.\nNo matchs when set '-query' with '@query_array'.\n";
        close OUT;
        exit(0);
      }

      print OUT "\nfinal matching SUM: $count_sum\n";
      close OUT;

      write_output_file();
      print STDOUT "'--evalue' result file: '$output_filter'.\n";
      exit(0);
    }
    exit(0);
  }
  if(%include_hash){
    unless(%query_search_include_name_hash){
      write_output_file();
      print STDOUT "After excuting:\n";
      foreach my $criteria (sort(keys %include_hash)){
        print STDOUT "  --$criteria\n";
      }
      print STDOUT "No annotated query left.\n";
      exit(0);
    }
    if($evalue){
      unless(%query_name_evalue_hash){
        write_output_file();
        print STDOUT "After excuting '--evalue $evalue', no annotated query left.\nPlease set the argument of '--evalue' in scope of $evalue_min-$evalue_max.\n";
        exit(0);
      }
      open OUT, '>>', $output_filter;
      $count_sum = scalar(keys %query_name_evalue_hash);
      print OUT "evalue: $evalue\nevalue matching SUM: $count_sum\n\n";
      foreach my $criteria (sort(keys %include_hash)){
        print OUT "$criteria\n";
      }
      $count_sum = scalar(keys %query_search_include_name_hash);
      print OUT "search_criteria matching SUM: $count_sum\n\n";
      foreach my $criteria (sort(keys %exclude_hash)){
        print OUT "$criteria\n";
      }
      print OUT "de_search_criteria matching left SUM: $count_left\n\n";
      print OUT "'-query' with '@query_array'\n\n" if @query_array;
      print OUT "-----------------------------------------------------------------\n\n";
      $count_sum = 0;
      foreach my $query_name (keys %query_name_evalue_hash){
        unless($query_search_exclude_name_hash{$query_name}){
          if($query_search_include_name_hash{$query_name}){
            if(@query_array){
              next unless $query_hash{$query_name};
            }
            write_query($query_name);
            $count_sum++;
            $print_last = 2;
            $db_query_name_hash{$query_name} += 1;
          }
        }
      }

      if($db_build && %db_query_name_hash){
        database_build();
      }

      if($print_last == 1){
        write_output_file();
        print STDOUT "No annotated query left.\n";
        print OUT "No annotated query left.\n";
        close OUT;
        exit(0);
      }
      write_output_file();
      print STDOUT "'--evalue' result file: '$output_filter'.\n";
      print OUT "\nfinal matching SUM: $count_sum\n";
      close OUT;
      exit(0);

    }
    unless($evalue){
      open OUT, '>>', $output_filter;
      foreach my $criteria (sort(keys %include_hash)){
        print OUT "$criteria\n";
      }
      $count_sum = scalar(keys %query_search_include_name_hash);
      print OUT "search_criteria matching SUM: $count_sum\n\n";
      foreach my $criteria (sort(keys %exclude_hash)){
        print OUT "$criteria\n";
      }
      print OUT "de_search_criteria matching left SUM: $count_left\n\n";
      print OUT "'-query' with '@query_array'\n\n" if @query_array;
      print OUT "-----------------------------------------------------------------\n\n";
      $count_sum = 0;
      foreach my $query_name (keys %query_search_include_name_hash){
        unless($query_search_exclude_name_hash{$query_name}){
          if(@query_array){
            next unless $query_hash{$query_name};
          }
          write_query($query_name);
          $count_sum++;
          $print_last = 2;
          $db_query_name_hash{$query_name} += 1;
        }
      }
      
      if($db_build && %db_query_name_hash){
        database_build();
      }

      if($print_last == 1){
        write_output_file();
        print STDOUT "No annotated query left.\n";
        print OUT "No annotated query left.\n";
        close OUT;
        exit(0);
      }
      write_output_file();
      print STDOUT "'--evalue' result file: '$output_filter'.\n";
      print OUT "\nfinal matching SUM: $count_sum\n";
      close OUT;
      exit(0);
    }
  }
  exit(0);
}

write_output_file();
exit(0);






sub database_build{


  open OUT1, '>>', $database_seq;
  open OUT2, '>>', $database_actsite;
  open OUT3, '>>', $database_AMB;
  open OUT4, '>>', $database_length;
  open OUT5, '>>', $database_name;
  open OUT6, '>>', $database_uniprot;

  my $count_id_db;
  my %db_query_count_hash;
  my %db_query_seq_hash;
  my $db_seq_write_key = 1;
  my $name_all;

  open IN, $annotation_seq;
  while(<IN>){
    chomp;
    $_ =~ s/\r//g;
    my $tmp = $_;
    if(/^>(.+)$/){
      $name_all = $1;
      $db_seq_write_key = 1;
      if($db_query_name_hash{$name_all}){
        $count_id_db++;
        print OUT5 "MER$count_id_db#$name_all\n";
        $db_query_count_hash{$name_all} = $count_id_db;
        print OUT1 ">MER$count_id_db \n";
        $db_seq_write_key = 2;
        next;
      }
    }
    if($db_seq_write_key == 2){
      print OUT1 "$_\n";
      my $tmp = $_;
      $tmp =~ s/\r//g;
      my $db_seq_length = length($tmp);
      $db_query_seq_hash{$count_id_db} += $db_seq_length;
    }
  }
  close IN;
  close OUT1;
  close OUT5;
  $name_all = '';
  my $cmd_makeblastdb_db = "makeblastdb -in $database_seq -dbtype prot -out $database_seq";
  system("$cmd_makeblastdb_db")== 0 or die("failed to execute command \'$cmd_makeblastdb_db\': $!\n");

  my @db_seq_length_array = sort{$a <=> $b}(keys %db_query_seq_hash);
  foreach my $seq_id_tmp (@db_seq_length_array){
    print OUT4 "$seq_id_tmp \{peptidase unit: 1\-$db_query_seq_hash{$seq_id_tmp}\}\n";
  }
  close OUT4;

  open IN, $annotation_query;
  my $db_AT_count = 1;
  my %db_AT_hash;
  my %db_AT_AMB_hash;
  
  while(<IN>){
    chomp;
    $_ =~ s/\r//g;
    if(/^QN   (.+)$/){
      $name_all = $1;
      $db_AT_count = 1;
      %db_AT_hash = ();
      %db_AT_AMB_hash = ();
    }
    if($db_query_count_hash{$name_all}){
      if(/^AT   ACT_SITE: (.+)$/){
        if($db_AT_count == 1){
          my $tmp = $1;
          my @db_AT_array = ($tmp =~ /[A-Za-z]\d+/g);
          foreach my $AT_tmp (@db_AT_array){
            if($AT_tmp =~ /^([A-Za-z])(\d+)$/){
              my $AT_AA = $1;
              my $AT_site = $2;
              $db_AT_AMB_hash{$AT_site} .= 'A';
              $db_AT_hash{$AT_site} = $AT_tmp;
            }
          }
        }
  
      }
      if(/^AT   METAL: (.+)$/){
        if($db_AT_count == 1){
          my $tmp = $1;
          my @db_AT_array = ($tmp =~ /[A-Za-z]\d+/g);
          foreach my $AT_tmp (@db_AT_array){
            if($AT_tmp =~ /^([A-Za-z])(\d+)$/){
              my $AT_AA = $1;
              my $AT_site = $2;
              $db_AT_AMB_hash{$AT_site} .= 'M';
              $db_AT_hash{$AT_site} = $AT_tmp;
            }
          }
        }
      }
      if(/^AT   BINDING: (.+)$/){
        if($db_AT_count == 1){
          my $tmp = $1;
          my @db_AT_array = ($tmp =~ /[A-Za-z]\d+/g);
          foreach my $AT_tmp (@db_AT_array){
            if($AT_tmp =~ /^([A-Za-z])(\d+)$/){
              my $AT_AA = $1;
              my $AT_site = $2;
              $db_AT_AMB_hash{$AT_site} .= 'B';
              $db_AT_hash{$AT_site} = $AT_tmp;
            }
          }
        }
      }
      if(/^UP   Uniport_id: (.+)$/){
        if($db_AT_count == 1){
          my $tmp = $1;
          print OUT6 "$db_query_count_hash{$name_all}#$tmp\n";
        }
      }
    }
    if(/^EV   / && $db_AT_count == 1){
      if(%db_AT_AMB_hash){
        my @AT_AMB_array = sort{$a <=> $b} (keys %db_AT_AMB_hash);
        print OUT3 "$db_query_count_hash{$name_all}\#\"";
        my $count_AT_AMB = scalar @AT_AMB_array;
        my $count_AT_AMB_tmp;
        foreach my $AT_AMB_tmp(@AT_AMB_array){
          $count_AT_AMB_tmp++;
          if($db_AT_AMB_hash{$AT_AMB_tmp} =~ /^A+$/){
            print OUT3 "A$AT_AMB_tmp";
          }elsif($db_AT_AMB_hash{$AT_AMB_tmp} =~ /^B+$/){
            print OUT3 "B$AT_AMB_tmp";
          }elsif($db_AT_AMB_hash{$AT_AMB_tmp} =~ /^M+$/){
            print OUT3 "M$AT_AMB_tmp";
          }elsif($db_AT_AMB_hash{$AT_AMB_tmp} =~ /^A+B+$/){
            print OUT3 "J$AT_AMB_tmp";
          }elsif($db_AT_AMB_hash{$AT_AMB_tmp} =~ /^A+M+$/){
            print OUT3 "O$AT_AMB_tmp";
          }elsif($db_AT_AMB_hash{$AT_AMB_tmp} =~ /^A+M+B+$/){
            print OUT3 "X$AT_AMB_tmp";
          }elsif($db_AT_AMB_hash{$AT_AMB_tmp} =~ /^M+B+$/){
            print OUT3 "Z$AT_AMB_tmp";
          }else{
            print STDERR "FATAL ERROR: Please redownload \'analysis.pl\'.\n";
            exit(1);
          }
          if($count_AT_AMB_tmp < $count_AT_AMB){
            print OUT3 "\, ";
          }
          if($count_AT_AMB_tmp == $count_AT_AMB){
            print OUT3 "\"\,\n";
          }
        }
        %db_AT_AMB_hash = ();
      }
      if(%db_AT_hash){
        my @AT_AMB_array = sort{$a <=> $b} (keys %db_AT_hash);
        print OUT2 "$db_query_count_hash{$name_all}\#\"";
        my $count_AT_AMB = scalar @AT_AMB_array;
        my $count_AT_AMB_tmp;
        foreach my $AT_AMB_tmp(@AT_AMB_array){
          $count_AT_AMB_tmp++;
          print OUT2 "$db_AT_hash{$AT_AMB_tmp}";
          if($count_AT_AMB_tmp < $count_AT_AMB){
            print OUT2 ", ";
          }
          if($count_AT_AMB_tmp == $count_AT_AMB){
            print OUT2 "\",\n";
          }
        }
        %db_AT_hash = ();
      }
    }
    if(/^DE   /){
      $db_AT_count = 2;
    }
  }
  close IN;
  close OUT2;
  close OUT3;
  close OUT6;
}

sub write_output_file{
  if($active_list){
    print STDOUT "'--active_list' result file: '$output/actsites'.\n"
  }
  if($metal_list){
    print STDOUT "'--metal_list' result file: '$output/actsites'.\n"
  }
  if($binding_list){
    print STDOUT "'--binding_list' result file: '$output/actsites'.\n"
  }
  if($active_binding_list){
    print STDOUT "'--active_binding_list' result file: '$output/actsites'.\n"
  }
  if($metal_binding_list){
    print STDOUT "'--metal_binding_list' result file: '$output/actsites'.\n"
  }
  if($active_metal_binding_list){
    print STDOUT "'--active_metal_binding_list' result file: '$output/actsites'.\n"
  }
  if($EC_list){
    print STDOUT "'--EC_list' result file: '$output/EC_number'.\n"
  }
  if($family_list){
    print STDOUT "'--family_list' result file: '$output/family'.\n"
  }
  if($name_list){
    print STDOUT "'--name_list' result file: '$output/name'.\n"
  }
  if($query_list){
    print STDOUT "'--query_list' result file: '$output/query'.\n"
  }
}



sub write_query{
  my $query_name_tmp = $_[0];
  my $i;
  print OUT "\nQN   $query_name_tmp\n";
  for($i = 1;$i <= $count_QN_max;$i++){
    my $tmp = $query_name_tmp . '_' . $i;
    if($ALL){
      if($query_result_all_hash{$tmp}){
        print OUT "$query_result_all_hash{$tmp}";
      }
    }
    my $uniprot_id = $uniprot_id_hash{$tmp};
    if($uniprot_id){
      if($reference || $note || %annotation_cross_hash || %active_site_hash || %metal_site_hash || %binding_site_hash){
        print OUT "UP   $uniprot_id\n";
      }
      if($reference){
        print OUT "$reference_hash{$uniprot_id}" if $reference_hash{$uniprot_id};
      }
      if($note){
        print OUT "$note_hash{$uniprot_id}" if $note_hash{$uniprot_id};
      }
      if(%annotation_cross_hash){
        print OUT "$annotation_cross_hash{$uniprot_id}" if $annotation_cross_hash{$uniprot_id};
      }
      if(%active_site_hash){
        print OUT "$active_site_hash{$uniprot_id}" if $active_site_hash{$uniprot_id};
      }
      if(%metal_site_hash){
        print OUT "$metal_site_hash{$uniprot_id}" if $metal_site_hash{$uniprot_id};
      }
      if(%binding_site_hash){
        print OUT "$binding_site_hash{$uniprot_id}" if $binding_site_hash{$uniprot_id};
      }
    }
  }
  

}


sub read_uniprot{
  my $in = $_[0];
  my $uniprot_id;
  open IN, $in;
  while(<IN>){
    chomp;
    $_ =~ s/\r//g;
    if(/^AC   (.+)$/){
      my $id_tmp = $1;
      my @id_array = ($id_tmp =~ /\w+/g);
      foreach my $id (@id_array){
        if($uniprot_name_hash{$id}){
          $uniprot_id = $id;
        }
      }
    }
    if($reference){
      if(/^RN   /){
        $reference_hash{$uniprot_id} .= "$_\n";
      }
      if(/^RX   /){
        $reference_hash{$uniprot_id} .= "$_\n";
      }
      if(/^RT   /){
        $reference_hash{$uniprot_id} .= "$_\n";
      }
      if(/^RL   /){
        $reference_hash{$uniprot_id} .= "$_\n";
      }
    }
    if($note){
      if(/^CC   /){
        $note_hash{$uniprot_id} .= "$_\n";
      }
    }
    if($annotation_all){
      if(/^DR   /){
        $annotation_cross_hash{$uniprot_id} .= "$_\n";
      }
    }else{
      if(/^DR   ([^\;]+)\;/){
        my $DR_tmp = $1;
        $annotation_list_hash{$DR_tmp} += 1;
        if($annotation_hash{$DR_tmp}){
          $annotation_cross_hash{$uniprot_id} .= "$_\n";
        }
      }
    }    
    if($active_site){
      if(/FT   ACT_SITE/){
        $active_site_hash{$uniprot_id} .= "$_\n";
      }
    }
    if($metal_site){
      if(/FT   METAL/){
        $metal_site_hash{$uniprot_id} .= "$_\n";
      }
    }
    if($binding_site){
      if(/FT   BINDING/){
        $binding_site_hash{$uniprot_id} .= "$_\n";
      }
    }
  }
}








#next if $next_switch == 1;
sub match_name{
  my $tmp_name = $_[0];
  $next_switch = 1;
  if($tmp_all =~ /#N#[^#]*$tmp_name/i){
    $next_switch = 2;
  }
}

sub match_family{
  my $tmp_family = $_[0];
  $next_switch = 1;
  if($tmp_family =~ /^[ASMCGPUNT]00$/){
    $tmp_family =~ s/0//g;
  }
  if($tmp_family =~ /^([ASMCGPUNT])([1-9])$/){
    $tmp_family = $1 . '0' . $2;
  }
  if($tmp_family =~ /^([ASMCGPUNT])[0]+([1-9]\d+)/){
    $tmp_family = $1 . $2;
  }
  $tmp_family = '#F#' . $tmp_family . '_';
  if($tmp_all =~ /$tmp_family/i){
    $next_switch = 2;
  }
}

sub match_EC{
  my $tmp_EC = $_[0];
  $next_switch = 1;
  $tmp_EC = '#E#' . $tmp_EC . '_';
  if($tmp_all =~ /$tmp_EC/i){
    $next_switch = 2;
  }
}

sub match_active{
  my $tmp_active = $_[0];
  $next_switch = 1;
  $tmp_active = '#A#' . $tmp_active . '_';
  if($tmp_all =~ /$tmp_active/i){
    $next_switch = 2;
  }
}

sub match_metal{
  my $tmp_metal = $_[0];
  $next_switch = 1;
  $tmp_metal = '#M#' . $tmp_metal . '_';
  if($tmp_all =~ /$tmp_metal/i){
    $next_switch = 2;
  }
}

sub match_binding{
  my $tmp_binding = $_[0];
  $next_switch = 1;
  $tmp_binding = '#B#' . $tmp_binding . '_';
  if($tmp_all =~ /$tmp_binding/i){
    $next_switch = 2;
  }
}










sub search_query{
  my $in = $_[0];
  my $query_name;
  open IN, $in;
  while(<IN>){
    chomp;
    $_ =~ s/\r//g;
    if(/^QN   (.+)$/){
      $query_name = $1;
      $query_name_search_all_hash{$query_name} += 1;
    }
    if(/^AT   ACT_SITE: (.+)$/){
      my $tmp = $1;
      $tmp =~ s/\s//g;
      $tmp =~ s/\d+//g;
      $tmp = '#A#' . $tmp . '_';
      $query_search_all_hash{$query_name} .= $tmp;
    }
    if(/^AT   METAL: (.+)$/){
      my $tmp = $1;
      $tmp =~ s/\s//g;
      $tmp =~ s/\d+//g;
      $tmp = '#M#' . $tmp . '_';
      $query_search_all_hash{$query_name} .= $tmp;
    }
    if(/^AT   BINDING: (.+)$/){
      my $tmp = $1;
      $tmp =~ s/\s//g;
      $tmp =~ s/\d+//g;
      $tmp = '#B#' . $tmp . '_';
      $query_search_all_hash{$query_name} .= $tmp;
    }
    if(/^FM   Family: ([^\.]+)/){
      my $tmp = $1;
      $tmp = '#F#' . $tmp . '_';
      $query_search_all_hash{$query_name} .= $tmp;
    } 
    if(/^DE\s+EC=([^\s\;\{]+)/){
      my $tmp = $1;
      $tmp = '#E#' . $tmp . '_';
      $query_search_all_hash{$query_name} .= $tmp;
    }
    if(/^DE\s+RecName: Full=([^\;\{]+)/){
      my $tmp = $1;
      $tmp =~ s/\s+/,/g;
      $tmp =~ s/_/,/g;
      $tmp = '#N#' . $tmp . '_';
      $query_search_all_hash{$query_name} .= $tmp;
    }
    if(/^DE\s+AltName: Full=([^\;\{]+)/){
      my $tmp = $1;
      $tmp =~ s/\s+/,/g;
      $tmp =~ s/_/,/g;
      $tmp = '#N#' . $tmp . '_';
      $query_search_all_hash{$query_name} .= $tmp;
    }
  }
  close IN;
}

sub read_query{
  my $in = $_[0];
  $count_QN_max = 1;
  my $count_AT = 1;
  my %AM_hash;
  my %AB_hash;
  my %BM_hash;
  my %ABM_hash;
  my $query_name;
  my $count_QN;
  my $count_end;
  my $query_name_num;
  open IN, $in;
  while(<IN>){
    chomp;
    $_ =~ s/\r//g;
    if(/^QN   (.+)$/){
      $query_name = $1;
      $count_QN = 1;
    }
    if(/^AT   ACT_SITE: (.+)$/){
      my $tmp = $1;
      $query_name_num = $query_name . '_' . $count_QN;
      $query_result_all_hash{$query_name_num} .= "$_\n";
      if($count_QN > $count_QN_max){
        $count_QN_max = $count_QN;
      }
      $tmp =~ s/\s//g;
      $tmp =~ s/\d+//g;
      $active_list_hash{$tmp} += 1;
      $AM_hash{$count_AT} .= $tmp . '_';
      $AB_hash{$count_AT} .= $tmp . '_';
      $ABM_hash{$count_AT} .= $tmp . '_';
      $count_end = 1;
    }
    if(/^AT   METAL: (.+)$/){
      my $tmp = $1;
      $query_result_all_hash{$query_name_num} .= "$_\n";
      $tmp =~ s/\s//g;
      $tmp =~ s/\d+//g;
      $metal_list_hash{$tmp} += 1;
      $BM_hash{$count_AT} .= $tmp . '_';
      $AM_hash{$count_AT} .= $tmp . '_';
      $ABM_hash{$count_AT} .= $tmp . '_';
    }
    if(/^AT   BINDING: (.+)$/){
      my $tmp = $1;
      $query_result_all_hash{$query_name_num} .= "$_\n";
      $tmp =~ s/\s//g;
      $tmp =~ s/\d+//g;
      $binding_list_hash{$tmp} += 1;
      $AB_hash{$count_AT} .= $tmp . '_';
      $BM_hash{$count_AT} .= $tmp . '_';
      $ABM_hash{$count_AT} .= $tmp . '_';
    }
    if(/^FM   Family: (.+)$/){
      my $tmp = $1;
      $query_result_all_hash{$query_name_num} .= "$_\n";
      $family_list_hash{$tmp} += 1;
    }
    if(/^DE   /){
      $query_result_all_hash{$query_name_num} .= "$_\n";
    }
    if(/^DE\s+EC=([^\s\;\{]+)/){
      my $tmp = $1;
      $EC_list_hash{$tmp} += 1;
    }
    if(/^DE\s+RecName: Full=([^\;\{]+)/){
      my $tmp = $1;
      $tmp =~ s/\s+/_/g;
      $name_list_hash{$tmp} += 1;
    }
    if(/^DE\s+AltName: Full=([^\;\{]+)/){
      my $tmp = $1;
      $tmp =~ s/\s+/_/g;
      $name_list_hash{$tmp} += 1;
    }
    if(/^MD   /){
      $count_AT++;
      $count_QN++;
      $query_result_all_hash{$query_name_num} .= "$_\n";
    }
    if(/^EV   /){
      $query_result_all_hash{$query_name_num} .= "$_\n";
    }
    if(/^SC   /){
      $query_result_all_hash{$query_name_num} .= "$_\n";
    }
    if(/^UP   Uniport_id: (.+)$/){
      my $id = $1;
      $uniprot_id_hash{$query_name_num} = $id;
      $query_result_all_hash{$query_name_num} .= "$_\n";
    }  

  }
  close IN;
  foreach my $count_tmp (keys %AM_hash){
    my $tmp = $AM_hash{$count_tmp};
    my $count_num = ($AM_hash{$count_tmp} =~ s/_/_/g);
    if($count_num == 2){
      $active_metal_list_hash{$AM_hash{$count_tmp}} += 1;
    }
  }
  foreach my $count_tmp (keys %AB_hash){
    my $tmp = $AB_hash{$count_tmp};
    my $count_num = ($AB_hash{$count_tmp} =~ s/_/_/g);
    if($count_num == 2){
      $active_binding_list_hash{$AB_hash{$count_tmp}} += 1;
    }
  }
  foreach my $count_tmp (keys %BM_hash){
    my $tmp = $BM_hash{$count_tmp};
    my $count_num = ($BM_hash{$count_tmp} =~ s/_/_/g);
    if($count_num == 2){
      $metal_binding_list_hash{$BM_hash{$count_tmp}} += 1;
    }
  }
  foreach my $count_tmp (keys %ABM_hash){
    my $tmp = $ABM_hash{$count_tmp};
    my $count_num = ($ABM_hash{$count_tmp} =~ s/_/_/g);
    if($count_num == 3){
      $active_metal_binding_list_hash{$ABM_hash{$count_tmp}} += 1;
    }
  }

}


sub statistics{
  my $in = $_[0];
  my $tag = $_[1];
  my $out = "$output/statistics";
  my $same_key = 2;
  if($tag eq 'statistics'){
    open OUT, '>>', $out;
  }
  $count_sum = scalar(keys %query_name_hash);
  print STDOUT "$count_sum queries are annotated with putative peptidase activity by PepTidy.\n";
  if($tag eq 'statistics'){
    print OUT "Results Summery:\n";
    print OUT "----------------------------------------------------------------\n";
    print OUT "$count_sum queries are annotated with putative peptidase activity by PepTidy.\n";
  }
  if(%bais_hash){
    $count_sum = scalar(keys %bais_hash);
    print STDOUT "$count_sum queries have more than one type of peptidase active site.\n";
    if($tag eq 'statistics'){
      print OUT "$count_sum queries have more than one type of peptidase active site.\n";
    }
  }
  if(%query_same_hash){
    $count_sum = scalar(keys %query_same_hash);
    print STDOUT "$count_sum queries are same with sequences in MEROPS database or custom database.\n";
    if($tag eq 'statistics'){
      print OUT "$count_sum queries are same with sequences in MEROPS database or custom database.\n";
    }
    if($count_sum == scalar(keys %query_name_hash)){
      $same_key = 1;
    }
  }
  if(%query_family_hash){
    foreach my $family(sort (keys %query_family_hash)){
      print STDOUT "$query_family_hash{$family} queries belongs to $family_sum_hash{$family}.\n";
      if($tag eq 'statistics'){
        print OUT "$query_family_hash{$family} queries belongs to $family_sum_hash{$family}.\n";
      }
    }
  }else{
    print STDERR "ERROR13: File '$input/*_annotation_results_query' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
    exit(1);
  }
  if($same_key == 2){
    if($evalue_max && $evalue_min){
      print STDOUT "The E_Value scope is from $evalue_min to $evalue_max\n\n\n";
      if($tag eq 'statistics'){
        print OUT "The E_Value scope is from $evalue_min to $evalue_max\n";
      }
    }else{
      print STDERR "ERROR14: File '$input/*_annotation_results_query' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
      exit(1);
    }
  }elsif($same_key == 1){
    print STDOUT "ALL queries are same with sequences in MEROPS database or custom database.\n";
    if($tag eq 'statistics'){
      print OUT "ALL queries are same with sequences in MEROPS database or custom database.\n";
    }
  }else{
    print STDERR "ERROR15: File '$input/*_annotation_results_query' is modified.\n Please run PepTidy angain to obtain the annotation result files.\n";
    exit(1);
  }

  if($tag eq 'statistics' && %bais_hash){
    print OUT "\n\nQueries with more than one type of putative peptidase active site:\n";
    print OUT "----------------------------------------------------------------\n";
    open IN, $in;
    my $query_name;
    my $print_key = 1;
    while(<IN>){
      chomp;
      $_ =~ s/\r//g;
      if(/^QN   (.+)$/){
        $query_name = $1;
        $print_key = 1;    
        if($bais_hash{$query_name}){
          $print_key = 2;
        }
      }
      if($print_key == 2){
        print OUT "$_\n";
      }
    }
    close IN;
  }
  close OUT if $tag eq 'statistics';

}














