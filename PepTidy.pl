#! /usr/bin/perl

eval 'exec perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

############################################################
#                                                          #
# PepTidy                                                  #
# Peptidases mining and annotation                         #
#                                                          #
# Authors: Yi Jiang                                        #
#                                                          #
# Contact   : yijiang@peptidy.cn                           #
# Bug report: bugs@peptidy.cn                              #   
#                                                          #
# Release date: April 15th, 2019                           #
#                                                          #
# This script is under the Artistic Licence                #
# https://opensource.org/licenses/artistic-license-2.0     #
#                                                          #
# Usage:                                                   #
# PepTidy.pl [OPTIONS] --in                                #
#                                                          #
############################################################

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

PepTidy is to mine and annotate peptidases.

SYNOPSIS

perl PepTidy.pl [OPTIONS] --in 

  --in    Protein sequences input in fasta format
          --in <String>


OPTIONS

    --help               Print this help message
                         --help
    --version            print version number of name_of_script.pl
                         --version

    Options related to annotation and mining
    -------------------------------------------------------------------------------------------
    --workingdir         Working directory
                         --workingdir <String>
                         Default is current directory
    --skip               To skip ckecking the standard PepTidy database or custom database files 
                         --skip
    --merops             All annotated results prefer to be consistent with the standard PepTidy
                         database or custom database.
                         (If protein A in database has the same sequence with query B, PepTidy.pl
                         will directly annotate B as A.)
                         --merops
    --db_build           To create a custom database for PepTidy, if set.
                         --db_build
    --db_name            To give a name to custom database. Only contain number, letter and
                         underline are vaild..
                         --db_name <String>
    --annotate           To create annotation files while using option '--db_build'
                         --annotate
                         Compatible with --db_build
    --custom             To run PepTidy by custom database.
                         --custom
    --database           The custom database name and the directory path to custom database files.
                         --database <String>
                         Example [--database path/to/database_files_directory/database_name]
    --plus               To run PepTidy by both MEROPS and custom database.
                         --plus

    Options related to Blastp
    --------------------------------------------------------------------------------------------
    --cpu                Number of CPUs that can be used
                         --cpu <Integer, >=1>
                         Default = '12'
    --evalue             E_value for blastp reslut filter out
                         --evalue <Real>
                         Default = '1e-5'
    --num_descriptions   Number of database sequences to show one-line descriptions for
                         --num_descriptions <Integer, >=0>
                         Default = '3'
    --num_alignments     Number of database sequences to show alignments for
                         --num_alignments <Integer, >=0>
                         Default = '300'
    --word_size          Word size for wordfinder algorithm
                         --word_size <Integer, >=2>
                         Default = '3'
    --gapopen            Cost to open a gap
                         --gapopen <Integer>
                         Default = '11'
    --gapextend          Cost to extend a gap
                         --gapextend <Integer>
                         Default = '1'
    --threshold          Minimum word score such that the word is added to the BLAST lookup table
                         --threshold <Real, >=0>
                         Default = '11'
    --matrix             Scoring matrix name (normally BLOSUM62)
                         --matrix <String>
                         Default is 'BLOSUM62'


DESCRIPTION

  Example:
    
    perl PepTidy.pl [OPTIONS] --in

ENDUSAGE

my $version = '1.0.0';                     # name_of_script.pl version number
my $printVersion = 0;                      # print version number, if set
my $query_input;                           # query protein sequences in fasta fromat
my $query_input_tmp;                       # query protein sequences in fasta fromat
my %query_all_name_hash;                   # query names of protein sequences
my $query_file_name;                       # query file name
my $currentDir = cwd();                    # working superdirectory where programme is called from
my $CPU = 12;                              # number of CPUs that can be used
my $workDir;                               # in the working directory results and temporary files are stored
my $help;                                  # to print usage
my $skip;
my $plus;
my $progress_now;                          # to accept standard input
my $merops_switch;                         # to prefer to be consistent with MEROPS
my $db_build;                              # to create MEROPS-like custom database for PepTidy, if set
my $db_name;                               # MEROPS-like custom database name. It can only contain numbers, letters and underlines
my $custom;                                # to run Peptidy by custom database;
my $database;                              # the directory containing all custom database files;
my $database_name;                         # database name for blastp searching
my $database_tmp;                          # blastp searching database
my %database_all_name_hash;                # sequence names in database
my $database_MEROPS_protease = $directory . '/MEROPS_database/protease_lib';
   # MEROPS proteases database for blastp searching
my $database_MEROPS_pepunit  = $directory . '/MEROPS_database/pepunit_lib';
   # MEROPS pepunit database for blastp searching
my $evalue = 1e-5;                         # e value for blastp reslut filter out
my $num_descriptions = 3;                  # In blastp reslut, Number of database sequences to show one-line descriptions for
my $num_alignments = 300;                  # In blastp reslut, Number of database sequences to show alignments for
my $word_size = 3;                         # In blastp alingment, Word size for wordfinder algorithm
my $gapopen = 11;                          # In blastp alingment, Cost to open a gap
my $gapextend = 1;                         # In blastp alingment, Cost to extend a gap
my $threshold = 11;                        # In blastp alingment, Minimum word score such that the word is added to the BLAST lookup table
my $matrix = 'BLOSUM62';                   # In blastp alingment, Scoring matrix name (normally BLOSUM62)
my %matrix_hash;                           # all valid parameters of option 'matrix';
%matrix_hash = ('BLOSUM62' => 1,
                'BLOSUM80' => 1,
                'BLOSUM45' => 1,
                'PAM30'    => 1,
                'PAM70'    => 1);
my $blastp_protease_out;                   # blastp result by searching against MEROPS 'protease_lib'
my $blastp_pepunit_out;                    # blastp result by searching against MEROPS 'pepunit_lib'
my $blastp_pepunit_out_seq;                
my $blastp_pepunit_out_first;
my $blastp_database_out_seq;
my $blastp_database_out;                   # blastp result by searching against custom database created by option '-database'
my %actsite_match_query_id_hash;
my $database_same_query;                   # Queries same with database sequeces
my $blastp_database_out_left;              # blastp result excluding all queries same with database sequeces
my $protease_same_query;                   # Queries same with MEROPS proteases
my $blastp_protease_out_left;              # blastp result excluding all queries same with MEROPS proteases
my $database_actsite_file;                 # custom actsite_file 
my $actsite_out_database;                  # actsites matching reslut from $blastp_database_out
my $actsite_out_protease;                  # actsites matching reslut from $blastp_protease_out
my $actsite_out_pepunit;                   # actsites matching reslut from $blastp_pepunit_out
my $switch_next = 0;                       # switch to excute next sub function
my $switch_option = 0;                     # switch for domain annotation peplines
my $actsite_out_protease_filter;           # actsites matching reslut after filtering
my $actsite_out_database_filter;           # actsites matching reslut after filtering
my $actsite_out_pepunit_filter;            # actsites matching reslut after filtering
my $actsite_out_database_rearrange;        # filtered actsites matching reslut after rearranging
my $actsite_out_protease_rearrange;        # filtered actsites matching reslut after rearranging
my $actsite_out_pepunit_rearrange;         # filtered actsites matching reslut after rearranging
my $actsite_out_database_delete_repeat;    # rearranged actsites matching reslut after deleting repeat
my $actsite_out_protease_delete_repeat;    # rearranged actsites matching reslut after deleting repeat
my $actsite_out_pepunit_delete_repeat;     # rearranged actsites matching reslut after deleting repeat
my $actsite_result_database;               # switch for actsites matching reslut
my $actsite_result_protease;               # switch for actsites matching reslut
my $actsite_result_pepunit;                # switch for actsites matching reslut
my $actsite_out_database_result1;          # actsites matching reslut for database
my $actsite_out_protease_result1;          # actsites matching reslut for protease
my $actsite_out_pepunit_result1;           # actsites matching reslut for pepunit
my $actsite_out_database_left1;            # actsites matching left
my $actsite_out_protease_left1;            # actsites matching left
my $actsite_out_pepunit_left1;             # actsites matching left
my $actsite_out_database_result2;
my $actsite_out_protease_result2;
my $actsite_out_pepunit_result2;
my $actsite_out_database_left2;
my $actsite_out_protease_left2;
my $actsite_out_pepunit_left2;
my $actsite_out_database_result3;
my $actsite_out_protease_result3;
my $actsite_out_pepunit_result3;
my $actsite_out_database_proof;
my $actsite_out_protease_proof;
my $actsite_out_pepunit_proof;
my $actsites_annotation_database;          # all actsites matching resluts for database
my $actsites_annotation_protease;          # all actsites matching resluts for protease
my $actsites_annotation_pepunit;           # all actsites matching resluts for pepunit
my $actsite_result;                        # actsite switch
my %MEROPS_code_hash;                      # MEROPS proteases code
my $actsite_protease_key;                  # switch for actsites annotation at protease lever
my $actsite_pepunit_key;                   # switch for actsites annotation at pepunit lever
my $balstp_result_for_domain_analysis;     # balstp result for domain annotation
my $count_max_actsite_repeat = 1;          # max annotation times for one query id
my $actsites_annotation_all;               # all actsite annotation result
my $actsites_rearrange_all;                # all actsite rearrange result
my $actsites_filter_all;
my $blastp_result_for_miss_domain_analysis;# balstp result with missed merops id for domain annotation
my $domain_annotation_match_result;        # domain annotation result for unmissed merops
my $domain_annotation_miss_match_input;    # 
my $domain_annotation_miss_match_result;
my $ii;
my $domain_annotation_miss_match_result_tmp;
my $MEROPS_protease_actsites_all = $directory . '/MEROPS_protease_actsites_all';
my $MEROPS_protease_length = $directory . '/MEROPS_protease_length';
my $active_metal_binding_AMB_sum = $directory . '/MEROPS_protease_actsites_AMB_sum';
my $mer_uniport_id_sum = $directory . '/MEROPS_uniport_annotation_sum';
my $uniport_annotation_description = $directory . '/uniport_proteases_description.txt';
my $MEROPS_miss_code = $directory . '/MEROPS_miss_code';
my $database_custom_database_seq_name;
my %custom_database_seq_name_hash;
my %actsite_AMB_result;
my %query_ori_id_hash;
my %uniport_id_DE_hash;
my %query_ori_uniport_id_hash;
my %query_id_annotation_seq_hash;
my $query_annotation_results;
my $uniport_annotation_results;
my $annotation_ori_count;
my $annotation_query_count;
my $annotation_uniport_count;
my $annotation_file;
my $custom_database_check_count;
my %active_metal_binding_hash = ('A' => 'A',
                                 'B' => 'B',
                                 'M' => 'M',
                                 'O' => 'AM',
                                 'J' => 'AB',
                                 'X' => 'AMB',
                                 'Z' => 'BM');

# A -> A; B -> B; M -> M; AM -> O; AB -> J; AMB -> X; BM -> Z;





if(@ARGV==0){
  print "$usage\n";
  exit(0);
}

GetOptions( 'help!'                         => \$help,
            'in=s'                          => \$query_input,
            'cpu=i'                         => \$CPU,#
            'workingdir=s'                  => \$workDir,
            'version!'                      => \$printVersion,
            'database=s'                    => \$database,
            'merops!'                       => \$merops_switch,
            'db_build!'                     => \$db_build,
            'db_name=s'                     => \$db_name,
            'custom!'                       => \$custom,
            'annotate!'                     => \$annotation_file,
            'skip!'                         => \$skip,
            'plus'                          => \$plus,

            'evalue=s'                      => \$evalue,#
            'num_descriptions=i'            => \$num_descriptions,#
            'num_alignments=i'              => \$num_alignments,#
            'word_size=i'                   => \$word_size,#
            'gapopen=i'                     => \$gapopen,#
            'gapextend=i'                   => \$gapextend,#
            'threshold=s'                   => \$threshold,#
            'matrix=s'                      => \$matrix) or die("$!");

if($help){
  print $usage;
  exit(0);
}

if($printVersion){
    print "PepTidy Version: $version\n";
    exit(0);
}


######################## make some regular checks ########################

print STDOUT "Check all files and directorys ...\n";
# if option '--annotate' is set, check whether it is set correct.
if($annotation_file){
  unless($db_build){  
    print STDERR "ERROR: option \'--annotate\' should be set with \'--db_build\'.\n$usage";
    exit(1);
  }
}

# if option '--plus' is set, check whether it is set correct.
if($plus){
  unless($custom && $database){
    print STDERR "ERROR: option \'--custom\' should be set with \'--custom\' and \'--database\'.\n$usage";
    exit(1);
  }
}

# if option '--custom' is set, check whether it is set correct.
if($custom){
  if($database){
    if($database =~ /^(.+)\/(\w+)$/){
      $database = $1;
      $database_name = $2;
      $database = rel2abs($database);
      if(! -e $database){
        print STDERR "ERROR: database directory \'$database\' does not exist. Please check.\n";
        exit(1);
      }
      if(! -d $database){
        print STDERR "ERROR: database directory \'$database\' does not exist. Please check.\n";
        exit(1);
      }
    }else{
      print STDERR "Option \'--database\' needs a argument contains the custom database name and the directory path to custom database files. Just like \'\/path_to_database_files_directory\/database_name\'.\n";
      exit(1);
    }
    print "Custom database directory: \'$database\'\nCustom database name: \'$database_name\'\n";
    unless($plus){
      $database_MEROPS_protease = $database . '/' . $database_name . '_database_seq';
      $database_custom_database_seq_name = $database . '/' . $database_name . '_seq_name';
      $MEROPS_protease_length = $database . '/' . $database_name . '_seq_length';
      $MEROPS_protease_actsites_all = $database . '/' . $database_name . '_seq_actsites_all';
      $active_metal_binding_AMB_sum = $database . '/' . $database_name . '_seq_actsites_AMB_sum';
      $mer_uniport_id_sum = $database . '/' . $database_name . '_seq_uniport_annotation_sum';
    }
  }
  unless($database){
    print STDERR "Please set option \'--database\' if running PepTidy.pl by custom database.\n$usage";
    exit(1);
  }
}
unless($custom){
  if($database){
    print STDERR "Please set option \'--custom\' if running PepTidy.pl by custom database.\n$usage";
    exit(1);
  }
}

#check all files belong to Peptidy

 #check whether files exist
sub Check_peptidy_files{
  my $file = $_[0];
  unless($file){
    print STDERR "FATAL ERROR: \'PepTidy.pl\' is modified. Please redownload this file and replace the old one.\n";
    exit(1);
  }
  if(! -e $file || -d $file){
    print STDERR "ERROR: No file \'$file\'. Please redownload this file.\n";
    exit(1);
  }
}

Check_peptidy_files($database_MEROPS_protease);
my $database_MEROPS_protease_blastdb = "$database_MEROPS_protease.phr";
Check_peptidy_files($database_MEROPS_protease_blastdb);
$database_MEROPS_protease_blastdb = "$database_MEROPS_protease.pin";
Check_peptidy_files($database_MEROPS_protease_blastdb);
$database_MEROPS_protease_blastdb = "$database_MEROPS_protease.psq";
Check_peptidy_files($database_MEROPS_protease_blastdb);

if(! defined($custom) || $plus){
  Check_peptidy_files($database_MEROPS_pepunit);
  my $database_MEROPS_pepunit_blastdb = "$database_MEROPS_pepunit.phr";
  Check_peptidy_files($database_MEROPS_pepunit_blastdb);
  $database_MEROPS_pepunit_blastdb = "$database_MEROPS_pepunit.pin";
  Check_peptidy_files($database_MEROPS_pepunit_blastdb);
  $database_MEROPS_pepunit_blastdb = "$database_MEROPS_pepunit.psq";
  Check_peptidy_files($database_MEROPS_pepunit_blastdb);
  Check_peptidy_files($MEROPS_miss_code);
}

Check_peptidy_files($MEROPS_protease_actsites_all);
Check_peptidy_files($MEROPS_protease_length);
Check_peptidy_files($active_metal_binding_AMB_sum);
Check_peptidy_files($mer_uniport_id_sum);
Check_peptidy_files($uniport_annotation_description);

if($custom){
  Check_peptidy_files($database_custom_database_seq_name);
}
print STDOUT "All Peptidy file paths are set correctly.\n";

unless($skip){
 #check the integrality
if($database_MEROPS_protease){
  open IN, $database_MEROPS_protease;
  my $check_count_all;
  my $check_count_seq;
  while(<IN>){
    $check_count_all++;
    chomp;
    if(/^>MER\d+\s/){
      $check_count_seq++;
    }
  }
  close IN;
  unless($check_count_all && $check_count_seq){
    print STDERR "ERROR: File \'$database_MEROPS_protease\' is modified. Please redownload or recreate this file.\n";
    exit(1);
  }
  unless($custom){
    unless($check_count_all == 6834201 && $check_count_seq == 710456){
      print STDERR "ERROR: File \'$database_MEROPS_protease\' is modified. Please redownload this file.\n";
      exit(1);
    }
  }
  if($custom){
    $custom_database_check_count = $check_count_seq;
  }
}
print STDOUT "File \'$database_MEROPS_protease\' is OK.\n";

 #check the integrality
unless($custom){
  if($database_MEROPS_pepunit){
    open IN, $database_MEROPS_pepunit;
    my $check_count_all;
    my $check_count_seq;
    while(<IN>){
      $check_count_all++;
      chomp;
      if(/^>MER\d+\s/){
        $check_count_seq++;
      }
    }
    close IN;
    unless($check_count_all && $check_count_seq){
      print STDERR "ERROR: File \'$database_MEROPS_pepunit\' is modified. Please redownload this file.\n";
      exit(1);
    }
    unless($check_count_all == 3981481 && $check_count_seq == 710405){
      print STDERR "ERROR: File \'$database_MEROPS_pepunit\' is modified. Please redownload this file.\n";
      exit(1);
    }
  }
  print STDOUT "File \'$database_MEROPS_pepunit\' is OK.\n";
}

 #check the integrality
if($uniport_annotation_description){
  open IN, $uniport_annotation_description;
  my $check_count_all;
  my $check_count_seq;
  while(<IN>){
    $check_count_all++;
    chomp;
    if(/^ID\s+/){
      $check_count_seq++;
    }
  }
  close IN;
  unless($check_count_all && $check_count_seq){
    print STDERR "ERROR: File \'$uniport_annotation_description\' is modified. Please redownload this file.\n";
    exit(1);
  }
  unless($check_count_all == 1412302 && $check_count_seq == 8771){
    print STDERR "ERROR: File \'$uniport_annotation_description\' is modified. Please redownload this file.\n";
    exit(1);
  }
}
print STDOUT "File \'$uniport_annotation_description\' is OK.\n";

 #check the integrality
if($MEROPS_protease_actsites_all){
  open IN, $MEROPS_protease_actsites_all;
  my $check_count_all;
  my $check_count_seq;
  while(<IN>){
    $check_count_all++;
    chomp;
    if(/^\d+\#\".+\"\,/){
      $check_count_seq++;
    }else{
      print STDERR "ERROR: File \'$MEROPS_protease_actsites_all\' is modified. Please redownload or recreate this file.\n";
      exit(1);
    }
  }
  close IN;
  unless($check_count_all && $check_count_seq){
    print STDERR "ERROR: File \'$MEROPS_protease_actsites_all' is modified. Please redownload or recreate this file.\n";
    exit(1);
  }
  unless($custom){
    unless($check_count_all == 710456 && $check_count_seq == 710456){
      print STDERR "ERROR: File \'$MEROPS_protease_actsites_all\' is modified. Please redownload this file.\n";
      exit(1);
    }
  }
  if($custom){
    unless($custom_database_check_count == $check_count_seq && $check_count_all == $check_count_seq){
      print STDERR "ERROR: File \'$MEROPS_protease_actsites_all\' is modified. Please redownload or recreate this file.\n";
    }
  }
}
print STDOUT "File \'$MEROPS_protease_actsites_all\' is OK.\n";

 #check the integrality
if($MEROPS_protease_length){
  open IN, $MEROPS_protease_length;
  my $check_count_all;
  my $check_count_seq;
  while(<IN>){
    $check_count_all++;
    chomp;
    if(/^\d+\s+\{peptidase\sunit:\s.+\}/){
      $check_count_seq++;
    }else{
      print STDERR "ERROR: File \'$MEROPS_protease_length\' is modified. Please redownload or recreate this file.\n";
      exit(1);
    }
  }
  close IN;
  unless($check_count_all && $check_count_seq){
    print STDERR "ERROR: File \'$MEROPS_protease_length\' is modified. Please redownload or recreate this file.\n";
    exit(1);
  }
  unless($custom){
    unless($check_count_all == 710456 && $check_count_seq == 710456){
      print STDERR "ERROR: File \'$MEROPS_protease_length\' is modified. Please redownload this file.\n";
      exit(1);
    }
  }
  if($custom){
    unless($custom_database_check_count == $check_count_seq && $check_count_all == $check_count_seq){
      print STDERR "ERROR: File \'$MEROPS_protease_actsites_all\' is modified. Please redownload or recreate this file.\n";
    }
  }
}
print STDOUT "File \'$MEROPS_protease_length\' is OK.\n";

}

 #check the integrality
my %actsite_AMB_mer_id_hash;
if($active_metal_binding_AMB_sum){
  open IN, $active_metal_binding_AMB_sum;
  my $check_count_all;
  my $check_count_seq;
  while(<IN>){
    $check_count_all++;
    chomp;
    if(/^(\d+)\#\"(.+)\"\,/){
      $check_count_seq++;
      my $mer_id = $1;
      my $actsite = $2;
      $actsite_AMB_mer_id_hash{$mer_id} = $actsite;
    }else{
      print STDERR "ERROR: File \'$active_metal_binding_AMB_sum\' is modified. Please redownload or recreate this file.\n";
      exit(1);
    }
  }
  close IN;
  unless($check_count_all && $check_count_seq){
    print STDERR "ERROR: File \'$active_metal_binding_AMB_sum\' is modified. Please redownload or recreate this file.\n";
    exit(1);
  }
  unless($custom){
    unless($check_count_all == 710456 && $check_count_seq == 710456){
      print STDERR "ERROR: File \'$active_metal_binding_AMB_sum\' is modified. Please redownload this file.\n";
      exit(1);
    }
  }
  unless($skip){
    if($custom){
      unless($custom_database_check_count == $check_count_seq && $check_count_all == $check_count_seq){
        print STDERR "ERROR: File \'$MEROPS_protease_actsites_all\' is modified. Please redownload or recreate this file.\n";
      }
    }
  }
}
print STDOUT "File \'$active_metal_binding_AMB_sum\' is OK.\n";

unless($skip){
#check the integrality
if($mer_uniport_id_sum){
  open IN, $mer_uniport_id_sum;
  my $check_count_all;
  my $check_count_seq;
  while(<IN>){
    $check_count_all++;
    chomp;
    if(/^\d+\#\w+$/){
      $check_count_seq++;
    }else{
      print STDERR "ERROR: File \'$mer_uniport_id_sum\' is modified. Please redownload or recreate this file.\n";
      exit(1);
    }
  }
  close IN;
  unless($check_count_all && $check_count_seq){
    print STDERR "ERROR: File \'$mer_uniport_id_sum\' is modified. Please redownload or recreate this file.\n";
    exit(1);
  }
  unless($custom){
    unless($check_count_all == 710456 || $check_count_seq == 710456){
      print STDERR "ERROR: File \'$mer_uniport_id_sum\' is modified. Please redownload this file.\n";
      exit(1);
    }
  }
  if($custom){
    unless($custom_database_check_count == $check_count_seq && $check_count_all == $check_count_seq){
      print STDERR "ERROR: File \'$MEROPS_protease_actsites_all\' is modified. Please redownload or recreate this file.\n";
    }
  }
}
print STDOUT "File \'$mer_uniport_id_sum\' is OK.\n";

#check the integrality
unless($custom){
  if($MEROPS_miss_code){
    open IN, $MEROPS_miss_code;
    my $check_count_all;
    my $check_count_seq;
    while(<IN>){
      $check_count_all++;
      chomp;
      if(/^\d+\s\w+\.\w+$/){
        $check_count_seq++;
      }else{
        print STDERR "ERROR: File \'$MEROPS_miss_code\' is modified. Please redownload this file.\n";
        exit(1);
      }
    }
    close IN;
    unless($check_count_all || $check_count_seq){
      print STDERR "ERROR: File \'$MEROPS_miss_code\' is modified. Please redownload this file.\n";
      exit(1);
    }
    unless($check_count_all == 11119 || $check_count_seq == 11119){
      print STDERR "ERROR: File \'$MEROPS_miss_code\' is modified. Please redownload this file.\n";
      exit(1);
    }
  }
  print STDOUT "File \'$MEROPS_miss_code\' is OK.\n";
}
}

#check the integrality
if($custom){
  if($database_custom_database_seq_name){
    open IN, $database_custom_database_seq_name;
    my $check_count_all;
    my $check_count_seq;
    while(<IN>){
      $check_count_all++;
      chomp;
      if(/^MER(\d+)\#(.+)$/){
        my $mer_id = $1;
        my $seq_name = $2;
        $custom_database_seq_name_hash{$mer_id} = $seq_name;
        $check_count_seq++;
      }else{
        print STDERR "ERROR: File \'$database_custom_database_seq_name\' is modified. Please recreate this file.\n";
        exit(1);
      }
    }
    close IN;
    unless($check_count_all && $check_count_seq){
      print STDERR "ERROR: File \'$database_custom_database_seq_name\' is modified. Please recreate this file.\n";
      exit(1);
    }
    unless($skip){
      unless($check_count_seq == $custom_database_check_count && $check_count_seq == $check_count_all){
        print STDERR "ERROR: File \'$database_custom_database_seq_name\' is modified. Please recreate this file.\n";
        exit(1);
      }
    }
  }
  print STDOUT "File \'$database_custom_database_seq_name\' is OK.\n";
}


# check whether query input is set
if(!defined $query_input){
  print STDERR "ERROR: No query protein sequences files specified. Please set one query file.\n$usage";
  exit(1);
}else{
  if($query_input =~ /\/$/){
    print STDERR "ERROR: Please set a file name without \'/\' at last.\n";
    exit(1);
  }
}

# check whether query input file exist
if($query_input){
  if(! -e $query_input){
    print STDERR "ERROR: query input file $query_input does not exist. Please check.\n";
    exit(1);
  }
  if($query_input =~ /\/([^\/]+)$/){
    $query_file_name = $1;
  }
  unless($query_input =~ /\//g){
    $query_file_name = $query_input;
  }
  $query_input = rel2abs($query_input);
}
print STDOUT "File \'$query_input\' is OK.\n";


# if no working directory is set, use current directory
if(!defined $workDir){
  $workDir = $currentDir;
}else{
  my $last_char = substr($workDir, -1);
  if($last_char eq "\/"){
    chop($workDir);
  }
  if(! -d $workDir){
    print STDERR "ERROR: No \'--out $workDir\'. Please specify an existed working directory.\n";
    exit(1);
  }
}

# check the write permission of $workDir before building of the work directory
if(! -w $workDir){
  print STDERR "ERROR: Do not have write permission for $workDir.\nPlease use command 'chmod' to reset permission or specify another working directory\n";
  exit(1);
}

# if '--db_build' is set, check whether it is set correct.
if($db_build){
  #if($database){
  #  print STDERR "ERROR: option \'--db_build\' is incompatible with \'--custom\' or \'--database\'.\n$usage";
  #  exit(1);
  #}
  unless($db_name){
    print STDERR "ERROR: Option \'--db_name\' should be set and requires an argument.\n";
    exit(1);
  }
  if($db_name){
    unless($db_name =~ /^\w+$/){
      print STDERR "ERROR: Argument of \'--db_name\' can only contain numbers, letters and underlines. Please check \'--db_name $db_name\'.\n";
      exit(1);
    }
    my $tmp_name = $query_file_name . '_tmp';
    if($db_name =~ /^$tmp_name$/){
      print STDERR "ERROR: Argument of \'--db_name\' can't be like \'$db_name\'. Please change it.\n";
      exit(1);
    }
  }
  $workDir = $workDir . '/' . $db_name;
  if(-e $workDir){
    if(-d $workDir){
      print STDOUT "WARNING: The custom database name is \'$db_name\'. Peptidy need create a directory \'$workDir\' to save resluts. Directory \'$workDir\' is already existed. Please make sure the directory contains no file with \'$db_name\' in name.\nIf you want exit PepTidy and make a check first, please enter 'no'.\nIf you want skip this check, please enter 'yes' by keyboard.\n";
      print STDOUT "Please enter [yes/no]:";
      chomp($progress_now = <STDIN>);
      $progress_now =~ s/\n//g;
      $progress_now =~ s/\r//g;
      unless($progress_now){
        $progress_now = 'enter';
      }
      if($progress_now =~ /enter/){
        print STDOUT "Please enter [yes/no]:";
        chomp($progress_now = <STDIN>);
      }
      if($progress_now =~ /^no$/i || $progress_now =~ /^\'no\'$/i){
        print STDOUT "Please make sure the directory \'$workDir\' contains no file with \'$db_name\' in name.\n";
        exit(0);
      }elsif($progress_now =~ /^yes$/i || $progress_now =~ /^\'yes\'$/i){
        print STDOUT "Skip this check.\n";
      }else{
        print STDERR "ERROR: Sorry, can't read your choice.\n";
        exit(1);
      }
    }
    if(! -d $workDir){
      my $cmd_mkdir_db = "mkdir $workDir";
      system("$cmd_mkdir_db")==0 or die("failed to execute command \'$cmd_mkdir_db\': $!\n");
    }
  }
  if(! -e $workDir){
    my $cmd_mkdir_db = "mkdir $workDir";
    system("$cmd_mkdir_db")==0 or die("failed to execute command \'$cmd_mkdir_db\': $!\n");
  }
}
unless($db_build){
  if($db_name){
    print STDERR "ERROR: If you need create a MEROPS-like custom database named \'$db_name\', please also set option \'--db_build\'.\n$usage";
    exit(1);
  }
}


print STDOUT "Done.\n";



         ########################### check whether options are set correctly ############################

print STDOUT "Check all options ...\n";

# check whether option 'cores' are set correctly
unless($CPU){
  print STDERR "ERROR: \'core\' can not be set 0.\n";
  exit(1);
}
if($CPU){
  check_integer($CPU, 'cores');
}
check_core();

# check whether option 'num_descriptions' are set correctly
if($num_descriptions){
  check_integer($num_descriptions, 'num_descriptions');
}

# check whether option 'num_alignments' are set correctly
unless($num_alignments){
  print STDERR "ERROR: \'num_alignments\' can not be set 0.\n";
  exit(1);
}
if($num_alignments){
  check_integer($num_alignments, 'num_alignments');
}

# check whether option 'gapopen' are set correctly
if($gapopen){
  check_integer($gapopen, 'gapopen');
}

# check whether option 'gapextend' are set correctly
if($gapextend){
  check_integer($gapextend, 'gapextend');
}

# check whether option 'threshold' are set correctly
if($threshold){
  check_numeric($threshold, 'threshold');
}
unless($threshold){
  print STDERR "ERROR: \'threshold\' can not be set 0.\n";
  exit(1);
}

# check whether option 'word_size' are set correctly
if($word_size){
  check_integer($word_size, 'word_size');
  if($word_size < 2 || $word_size > 7){
    print STDERR "ERROR: \'word_size\' should be set between 2 to 7.\n";
    exit(1);
  }
  print STDOUT "check option \'word_size\' OK.\n";
}

# check whether option 'evalue' are set correctly
if($evalue){
  check_numeric($evalue, 'evalue');
}
unless($evalue){
  print STDERR "ERROR: \'evalue\' can not be set 0.\n";
  exit(1);
}

# check whether option 'matrix' are set correctly
if($matrix){
  unless($matrix_hash{$matrix}){
    print STDERR "ERROR: \'$matrix\' is not a valid option for \'matrix\'.Please use \'BLOSUM62\\80\\45\' 'PAM30\\70'.\n";
    exit(1);
  }
}
unless($matrix){
  print STDERR "ERROR: \'matrix\' can not be set 0.\n";
  exit(1);
}

# check whether option 
print STDOUT "Dnoe.\n";

         ########################### some checks beforehand ############################

# check integer options are set correctly
sub check_integer{
  my $integer_in = $_[0];
  my $option_tag = $_[1];
  if($integer_in =~ /[^0-9]/){
    print STDERR "ERROR: $option_tag only accept integer.\n";
    exit(1);
  }
  if($option_tag ne 'word_size'){
    print STDOUT "check option \'$option_tag\' OK.\n";
  }
}
# check_integer($option, tag);

# check numeric options are set correctly
sub check_numeric{
  my $number = $_[0];
  my $option_tag = $_[1];
  unless($number > 0 || $number < 0){
    print STDERR "ERROR: $option_tag only accept number.\n";
    exit(1);
  }
  print STDOUT "check option \'$option_tag\' OK.\n";
}

# next step switch
#check_switch_next(); 
sub check_switch_next{
  unless($switch_next){
    print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
    exit(1);
  }
}


      ################## subfunction for Blastp ##################
sub delete_file{
  my $in = $_[0];
  my $in_tmp = $in;
  my $last_char = substr($in_tmp, -1);
  if($last_char eq '*'){
    chop($in_tmp);
  }
  if($in_tmp){
    if(-e $in_tmp){
      unless(-d $in_tmp){
        my $cmd_delete = "rm $in";
        system("$cmd_delete")==0 or print STDOUT "WARNING: PepTidy can't delete \'$in\'. Please check manually after the end of running of PepTidy.pl\n";
      }else{
        print STDOUT "WARNING: \'$in_tmp\' is a directory. PepTidy can't delete it. Please check manually after the end of running of PepTidy.pl\n";
      }
    }
  }
}

sub Check_fasta_format{
  $switch_next = 0;
  print STDOUT "Check fasta format ...\n";
  my $in = $_[0];
  my $count_line;
  my $count_start;
  my %query_name_hash;
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
      my @ERROR_chr = ($tmp_seq =~ /[^ACDEFGHIKLMNPQRSTVWY]/ig);
      if(@ERROR_chr){
        my $count_line_tmp;
        $count_line_tmp = $count_line + $count_add -1;
        print STDOUT "WARNING: Unexpected string \'@ERROR_chr\' at $count_line_tmp line of $in. It may affect the Blastp result.\n";
      }
    }
  }
  close IN;
  print STDOUT "Done.\n";
  $switch_next = 1;
}

sub Reformat_fasta_file{
  $switch_next = 0;
  print STDOUT "Reformat fasta file ...\n";
  my $in = $_[0];
  my $fasta_out;
  my $count_name;
  $fasta_out = $workDir  . '/' . $query_file_name . '_tmp';
  $query_input_tmp = $fasta_out;
  open IN, $in;
  open OUT, '>>', $fasta_out;
  while(<IN>){
    if(/^\n$/){
      next;
    }
    chomp;
    my $tmp_all = $_;
    if(/^>(.+)$/){
      my $name_tmp = $1;
      $count_name++;
      print OUT ">QER$count_name\n";
      $query_all_name_hash{$count_name} = $name_tmp;
      $switch_next = 1;
    }
    unless(/^>/){
      print OUT "$tmp_all\n";
      $switch_next = 1;
    }
  }
  close IN;
  close OUT;
  print STDOUT "Done.\n";
}
# Reformat_fasta_file();     if query file

sub MakeBlastdb{
  my $in = $workDir  . '/' . $db_name . '_database_seq';
  my $cmd_MakeBlastdb = "makeblastdb -in $in -dbtype prot -out $in";
  print STDOUT "Make blastp database from \'$in\' ...\n";
  system("$cmd_MakeBlastdb")==0 or die("failed to execute command \'$cmd_MakeBlastdb\': $!\n");
  print STDOUT "Done.\n";
}

sub Blastp{
  $switch_next = 0;
  my $query_in = $_[0];
  my $database_in = $_[1];
  my $key_tag = $_[2];
  my $blastp_out;
  if($key_tag eq 'database'){
    $blastp_out = $query_in . '_database_blastp';
    $blastp_database_out = $blastp_out;
  }
  if($key_tag eq 'protease'){
    $blastp_out = $query_in . '_protease_blastp';
    $blastp_protease_out = $blastp_out;
  }
  elsif($key_tag eq 'pepunit'){
    $blastp_out = $query_in . '_pepunit_blastp';
    $blastp_pepunit_out = $blastp_out;
  }
  else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }
  my $cmd_Blastp = "blastp -query $query_in -db $database_in -num_threads $CPU -evalue $evalue -num_descriptions $num_descriptions -num_alignments $num_alignments -out $blastp_out -word_size $word_size -gapopen $gapopen -gapextend $gapextend -matrix $matrix -threshold $threshold";
  print STDOUT "Blastp \'$query_in\' against database \'$database_in\' ...\nThis may take a long time. Please wait ...\n";
  system("$cmd_Blastp")==0 or die("failed to execute command \'$cmd_Blastp\': $!\n");
  print STDOUT "Done.\n";
  $switch_next = 1;
}
# Blastp( $query_input_tmp, $database_MEROPS_protease || $database_MEROPS_pepunit, 'protease' || 'pepunit');    if query file

my $count_hits_match;
#Fetch_hits_sequence($blastp_pepunit_out, $query_input_tmp, $custom_out_name);
#if $switch_next == 1 hits
#if $switch_next == 2 no hits
sub Fetch_hits_sequence{
  print STDOUT "Fetch balstp hits sequences ...\n";
  $switch_next = 0;
  my $in1 = $_[0];
  my $in2 = $_[1];
  my %query_id_hash;
  my $query_id;
  my $count_hits;
  my $print_key = 1;
  my $out = $in2 . '_pepunit_blastp_hits_seq';
  $_[2] = $out;
  open IN, $in1;
  while(<IN>){
    chomp;
    if(/^Query=\s(.+)/){
      $query_id = $1;
    }
    if(/No\shits\sfound/){
      next;
    }
    if(/Sequences\sproducing\ssignificant\salignments/){
      $query_id_hash{$query_id} += 1;
      $count_hits++;
    }
  }
  close IN;

  unless($count_hits){
    print STDOUT "No balstp hits.\n";
    $switch_next = 2;
  }

  if($count_hits){
    $count_hits_match = $count_hits;
    open IN, $in2;
    open OUT, '>>', $out;
    while(<IN>){
      chomp;
      my $tmp_all = $_;
      if(/^>(.+)/){
        $query_id = $1;
        if($query_id_hash{$query_id}){
          $print_key = 2;
        }
        unless($query_id_hash{$query_id}){
          $print_key = 1;
        }
      }
      if($print_key == 2){
        print OUT "$tmp_all\n";
        $switch_next = 1;
      }
    }
    close IN;
    close OUT;
    print STDOUT "Done.\n";
  }

}




      ################## subfunction for Actsites Annotation ################

#if $switch_next == 1 hits
#if $switch_next == 2 no hits
sub Check_blastp_input{

  $switch_next = 0;
  my $in = $_[0];
  print STDOUT "Check blastp result file: \'$in\' ...\n";

  # check whether blastp result file exit
  if($in){
    if(! -e $in){
      print STDERR "ERROR: blastp result file \'$in\'does not exit. Please check.\n";
      exit(1);
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  open IN, $in;
  my $count_line;
  my $count_seq;
  my $count_seq_sub;
  my $count_query = 1;
  my $count_aln;
  my $count_seq_sub_end;
  while(<IN>){
    chomp;
    $count_line++;
    if($count_line == 1){
      unless(/^BLASTP\s/){
        print STDERR "ERROR at line $count_line: blastp result dose not have a proper format. \n";
        exit(1);
      }
      if(/^BLASTP\s(\d)/){
        my $version = $1;
        if($version == 1){
          print STDERR "ERROR: blastp version is too old. Plesae reinstall a new version blast+.\n";
          exit(1);
        }
      }
    }
    if(/^Query=\s/){
      if($count_query){
        if($count_query > 1){
          if($count_seq == 4){
          }
          elsif($count_seq == 3 && $count_seq_sub_end == 4){  
          }else{
            print STDERR "ERROR at line $count_line: blastp result dose not have a proper format.\n";
            exit(1);
          }
        }
      }
      if(/^Query=\sQER\d+/){
        $count_seq = 1;
        $count_query++;
      }else{
        print STDERR "ERROR at line $count_line: blastp result dose not have a proper format. Query title is wrong.\n";
        exit(1);
      }
    }
    if($count_query > 0 && $count_seq){
      if($count_seq == 1){
        if(/^Length=\d+/){
          $count_seq = 2;
        }
      }
    }
    if(/No\shits\sfound/){
      if($count_seq == 2){
        $count_seq = 4;
      }else{
        print STDERR "ERROR at line $count_line: blastp result dose not have a proper format.\n";
        exit(1);
      }
    }
    if(/Sequences\sproducing\ssignificant\salignments/){
      $switch_next = 1;
      if($count_seq == 2){
        $count_seq = 3;
        $count_aln++;
      }else{
        print STDERR "ERROR at line $count_line: blastp result dose not have a proper format.\n";
        exit(1);
      }
    }
    if(/^>\s/){
      unless($count_seq == 3){
        print STDERR "ERROR at line $count_line: blastp result dose not have a proper format.\n";
        exit(1);
      }
      if($count_seq_sub){
        unless($count_seq_sub_end == 4){
          print STDERR "ERROR at line $count_line: blastp result dose not have a proper format.\n";
          exit(1);
        }
      }
      if(/^>\sMER\d+/){
        $count_seq_sub = 1;
        $count_seq_sub_end = 1;
      }else{
        print STDERR "ERROR at line $count_line: blastp result dose not have a proper format. Sbjct title is wrong.\n";
        exit(1);
      }
    }
    if($count_query > 1 && $count_seq_sub){
      if($count_seq_sub == 1 && $count_seq == 3){
        if(/^Length=\d+/){
          $count_seq_sub = 2;
        }
      }
    }
    if(/^Query\s+\d+.+?\s+\d+$/){
      if($count_seq_sub == 2 && $count_seq == 3){
        $count_seq_sub = 3;
      }else{
        print STDERR "ERROR at line $count_line: blastp result dose not have a proper format.\n";
        exit(1);
      }
    }
    if(/^Sbjct\s+\d+.+?\s+\d+$/){
      if($count_seq_sub == 3 && $count_seq == 3){
        $count_seq_sub = 2;
        $count_seq_sub_end = 4;
      }
    }
  }
  close IN;
  unless($switch_next){
    $switch_next = 2;
  }
  print "Done.\n";
}


my $switch_ori; #if $switch_ori == 1, query same with merops sequences. 
# Filter_MEROPS_protease(,'database'/'protease');
sub Filter_MEROPS_protease{
  
  print STDOUT "Filter sequences in database from blastp result ...\n";
  $switch_next = 0;
  my $in = $_[0];
  my $key_tag = $_[1];
  my $out1 = $in . '_MEROPS_Query';

  if($in){
    if(! -e $in){
      print STDERR "ERROR: blastp result file \'$in\' does not exit. Please check.\n";
      exit(1);
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  if($key_tag eq 'database'){
    $database_same_query = $out1;
  }elsif($key_tag eq 'protease'){
    $protease_same_query = $out1;
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }
  open IN, $in;
  open OUT, '>>', $out1;
  my %MEROPS_protease_in_hash;
  my %MEROPS_protease_in_sbjct_hash;
  my $query_id;
  my $key_query =3;
  my $length_query;
  my $length_sbjct;
  my $mer_id;
  my $count_hits;
  my $count_same_seq;
  while(<IN>){
    chomp;
    if(/^Query=\sQER(\d+)/){
      $query_id = $1;
      $key_query = 1;
    }
    if(/^Length=(\d+)/){
      if($key_query == 1){
        $length_query = $1;
      }
      if($key_query == 2){
        $length_sbjct = $1;
      }
    }
    if(/Sequences\sproducing\ssignificant\salignments/){
      $count_hits++;
    }
    if(/^>\sMER(\d+)/){
      $mer_id = $1;
      $key_query = 2;
    }
    if(/^\sIdentities\s=\s(\d+)\/(\d+)\s/){
      my $aln_identity = $1;
      my $aln_length = $2;
      if($aln_identity && $aln_length && $length_query && $length_sbjct){
        if($aln_identity == $aln_length && $aln_identity == $length_query && $aln_identity == $length_sbjct){
          print OUT "QER$query_id#MER$mer_id\n";
          $switch_ori = 1;
          $MEROPS_protease_in_hash{$query_id} += 1;
          if($MEROPS_protease_in_hash{$query_id} == 1){
            $count_same_seq++;
            $MEROPS_protease_in_sbjct_hash{$query_id} = $mer_id;
          }
          if($MEROPS_protease_in_hash{$query_id} > 1){
            my $tmp = ' ' . $mer_id;
            $MEROPS_protease_in_sbjct_hash{$query_id} .= $tmp;
          }
        }      
      }
    }
  }
  close IN;
  close OUT;

  unless($count_same_seq){
    $switch_next = 1;
  }
  if($count_same_seq){
    if($count_same_seq < $count_hits){
      $switch_next = 1;
    }
  }

  print STDOUT "Done.\n";

}
# Filter_MEROPS_protease(,'database'/'protease');


# Match_actsite(,'database'/'protease'/'pepunit');
sub Match_actsite{

  $switch_next = 0;
  my $in = $_[0];
  my $key_tag = $_[1];
  my $in1;
  my $in2;
  my $out;
  if($key_tag eq 'database'){
    $in1 = $database_actsite_file; #lack
    
    $out = $query_input_tmp . '_actsite_database';
    $actsite_out_database = $out;
  }
  elsif($key_tag eq 'protease'){
    $in1 = $MEROPS_protease_actsites_all;
    $in2 = $MEROPS_protease_length;
    $out = $query_input_tmp . '_actsite_protease';
    $actsite_out_protease = $out;
  }
  elsif($key_tag eq 'pepunit'){
    $in1 = $MEROPS_protease_actsites_all;
    $in2 = $database_MEROPS_pepunit;
    $out = $query_input_tmp . '_actsite_pepunit';
    $actsite_out_pepunit = $out;
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  print STDOUT "Match actsites in \'$in\'\n";
 
  # check whether blastp result file exit
  if($in){
    if(! -e $in){
      print STDERR "ERROR: blastp result file \'$in\' does not exit. Please check.\n";
      exit(1);
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  if($in1 && $in2){
    if(! -e $in1){
      print STDERR "ERROR: \'$in1\' does not exit. Please make sure 'PepTidy.pl' is not moved to another directory.\n";
      exit(1);
    }
    if(! -e $in2){
      print STDERR "ERROR: \'$in2\' does not exit. Please make sure 'PepTidy.pl' is not moved to another directory.\n";
      exit(1);
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }
  print STDOUT "Check all necessary files OK.\n";
 
  open IN1, $in1;
  my %actsite_hash;
  my $id; 
  #step1 Read actsites into hash '%actsite_hash'('mer_id' as key, 'actsite' as value)
  while(<IN1>){
    chomp;
    my @actsite_array = split /#/, $_;
    $actsite_hash{$actsite_array[0]} = $actsite_array[1];     
  }
  close IN1;

  #Check whether actsites are read into hash '%actsite_hash'
  my @number_of_actsite = keys %actsite_hash;
  if(@number_of_actsite){
    print STDOUT "Read ". @number_of_actsite ." active and metal binding sites in\n";
  }else{
    die "Cannot read file '$in1'. Please make sure the correct format.\n";
  }


  $/ = "\n";
  my ($tmp,$count_nu,$count_unit,$gap,$unit_end_middle,$unit_start_middle,$middle);
  my %unit_hash;
  open IN2, $in2;

  #step2 read proteases length or pepunit sites into hash '%unit_hash'
  while(<IN2>){
    chomp;
    if($key_tag eq 'protease'){
      if(/^(\d+)\s+\{peptidase\sunit:\s(.+)\}/){
        $id = $1;
        $unit_hash{$id} = $2;
      }
    }elsif($key_tag eq 'pepunit'){
      if(/^>MER(\d+)\s.*unit:\s(.+)\}/){
        $id = $1;$unit_hash{$id} = $2;
      }
    }elsif($key_tag eq 'database'){
      #lack
    }else{
      print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
      exit(1);
    }
  }
  $id = '';
  close IN2; 
 
  #Check whether protease length or pepunit sites are read into hash '%unit_hash'
  my @number_of_pepunit = keys (%unit_hash);
  if(@number_of_pepunit){
    print STDOUT "Read ". @number_of_pepunit ." pepunits or protreases in\n";
    print STDOUT "Actsite matching...\n";
  }else{
    die "Cannot read file \'$in2\'. Please make sure the correct format.\n";
  }

  #step3 match actsites in blastp result
  $/ = "\n";
  my ($query_id,$mer_id,$active_sites,$i,$i1,$unit_start,$unit_end);
  my ($query_start,$query,$query_end,$sbjct_start,$sbjct,$sbjct_end);
  my ($active_sbjct,$active_query,$active_sites_query,$length_active,$query_tmp,$count_tmp);
  my ($count1,$count2,$count3,$count4);
  my ($active_count,$active_count_copy,$count_);
  my $acsite_count = 0;
  my $count_query = 2;
  my @units_array;
  my @active_sites_array;
  my($active_sites_array,$active_reslut_sbjct,$active_reslut_query);
  my @active_reslut_sbjct_array;
  my @active_reslut_query_array;
  my $query_id_1;

  open IN, $in;
  open OUT, '>>', $out;

  while(<IN>){
    chomp;
    #read each query id
    if(/^Query=\s(.+)/){
      $query_id_1 = $1;
      $count_query = 0;
      $active_sites = '';
      @active_sites_array = ();
    }

    #match actsites in every pairwise alingment
    if(/^>\sMER(\d+)\s/){
      $mer_id = 0;
      #print out actsites matching result
      if(@active_reslut_sbjct_array){print OUT "MER$id\tactive_site\t@active_reslut_sbjct_array\t";}
      unless(@active_reslut_sbjct_array){print OUT "MER$id\tactive_site_MISS\t$query_id\n" if $id;}
      if(@active_reslut_query_array){
        print OUT "$query_id\tactive_site\t@active_reslut_query_array\t";
        $switch_next = 1;
      }

      $id = $1;
      if(@active_reslut_sbjct_array or @active_reslut_query_array){print OUT "\n";}

      $query_id = $query_id_1; 
      $active_sites = '';
      @active_sites_array = ();
      @active_reslut_sbjct_array = ();
      @active_reslut_query_array = ();
      $tmp = $unit_hash{$id};

      #read protease length or pepunit sites from hash '在%unit_hash'
      @units_array = ($unit_hash{$id} =~ /\d+\-\d+/g);
      $count_nu = scalar(@units_array);
      $count_unit = 0;
      $gap = 0;

      #if only one unit, read start and end site
      if($count_nu == 1){
        if($units_array[0] =~ /(\d+)\-(\d+)/){
          $unit_start = $1;
          $unit_end = $2;
        }
      }

      #read actsite from hash '%actsite_hash'
      if($actsite_hash{$id}){    
        #actsites
        if($actsite_hash{$id} =~ /^\".+\",$/){
          #single actsite
          if($actsite_hash{$id} =~ /^\"[A-Z](\d+)\",$/){
            $active_sites = $1;
            $acsite_count = 5;
          }
          #poly actsite
          if($actsite_hash{$id} =~ /^\"[A-Z](\d+,\s[A-Z].+?)\",$/){
            @active_sites_array = split /,\s[A-Z]/, $1;
            @active_sites_array = sort {$a <=> $b}@active_sites_array;
            $acsite_count = 6;
          }
        }
        #no actsites;
        if($actsite_hash{$id} =~ /^,$/){
          $count4++;
          $acsite_count = 9;
        }
      }

      #if more than one unit, read start and end site
      if($count_nu > 1){
        $unit_start = 0;
        $unit_end_middle = 0;
        $unit_start_middle = 0;
        $unit_end = 0;

        if($active_sites){
          $count_unit = 0;
          $gap = 0;
          $unit_start = 0;
          $unit_end_middle = 0;
          $unit_start_middle = 0;
          $unit_end = 0;

          #change unit sites if there is gap between each units
          foreach(@units_array){
            #read the start site from the first unit
            if($count_unit == 0){
              if(/(\d+)\-(\d+)/){
                $unit_start = $1;
                $unit_end_middle = $2;
                my $unit_start1 = $unit_start - 1;
                my $unit_end_middle1 = $unit_end_middle -1;
                #if actsites in first unit
                if($active_sites > $unit_start1 and $active_sites < $unit_end_middle1){
                  last;
                }
              }
            }

            #second unit to last two unit
            if($count_unit > 0 and $count_unit < $count_nu){
              if(/(\d+)\-(\d+)/){
                $unit_start_middle = $1;
                my $unit_start1 = $unit_start_middle - 1;
                $middle = $unit_start_middle - $unit_end_middle -1;
                $unit_end_middle = $2;
                my $unit_end_middle1 = $unit_end_middle -1;
                print "WARING: overlap\t$id\t$query_id\tA\n" if $middle < 0;
                $gap += $middle if $middle > 0;
                #if actsites in any other units but not first and last unit
                if($active_sites > $unit_start1 and $active_sites < $unit_end_middle1){
                  $active_sites -= $gap;
                  last;
                }
              }
            }

            #read end site from last unit site 
            if($count_unit == $count_nu){
              if(/(\d+)\-(\d+)/){
                $unit_start_middle = $1;
                $unit_end = $2;
                my $unit_start1 = $unit_start_middle - 1;
                my $unit_end_middle1 = $unit_end -1;
                $gap += ($unit_start_middle - $unit_end_middle  -1);
                #if actsites in last unit
                if($active_sites > $unit_start1 and $active_sites < $unit_end_middle1){
                  $active_sites -= $gap;
                  last;
                }
              }
            }
            #count unit
            $count_unit++;
          }
        } 
      
        if(@active_sites_array){
          my @tmp_active_sites_array;
          #read each actsite from array '@active_sites_array'
          foreach(@active_sites_array){
            $count_unit = 1;
            $gap = 0;
            $active_sites = $_;
            $unit_start = 0;
            $unit_end_middle = 0;
            $unit_start_middle = 0;
            $unit_end = 0;
 
            #change unit sites if there is gap between each units
            foreach(@units_array){
              #read the start site from the first unit
              if($count_unit == 1){
                if(/(\d+)\-(\d+)/){
                  $unit_start = $1;
                  $unit_end_middle = $2;
                  my $unit_start1 = $unit_start - 1;
                  my $unit_end_middle1 = $unit_end_middle -1;
                  #if actsites in first unit
                  if($active_sites > $unit_start1 and $active_sites < $unit_end_middle1){
                    push @tmp_active_sites_array, $active_sites;
                    last;
                  }
                }
              }
              #second unit to last two unit
              if($count_unit > 1 and $count_unit < $count_nu){
                if(/(\d+)\-(\d+)/){
                  $unit_start_middle = $1;
                  my $unit_start1 = $unit_start_middle - 1;
                  $middle = $unit_start_middle - $unit_end_middle -1;
                  $unit_end_middle = $2;
                  my $unit_end_middle1 = $unit_end_middle -1;
                  print "WARING: overlap\t$id\t$query_id\tA\n" if $middle < 0;
                  $gap += $middle if $middle > 0;
                  #if actsites in any other units but not first and last unit
                  if($active_sites > $unit_start1 and $active_sites < $unit_end_middle1){
                    $active_sites -= $gap;
                    push @tmp_active_sites_array, $active_sites;
                    last;
                  }
                }
              }
              #read end site from last unit site
              if($count_unit == $count_nu){
                if(/(\d+)\-(\d+)/){
                  $unit_start_middle = $1;
                  $unit_end = $2;
                  my $unit_start1 = $unit_start_middle - 1;
                  my $unit_end_middle1 = $unit_end -1;
                  $gap += ($unit_start_middle - $unit_end_middle  -1);
                  #if actsites in last unit
                  if($active_sites > $unit_start1 and $active_sites < $unit_end_middle1){
                    $active_sites -= $gap;
                    push @tmp_active_sites_array, $active_sites;
                    last;
                  }
                }
              }
              #count unit
              $count_unit++;
            }
          }
          #sort all changed actsites and store into array '@active_sites_array'
          @active_sites_array = sort {$a <=> $b}@tmp_active_sites_array;
          $active_sites = '';
        }


      } #if($count_nu > 1){
    } #if(/^>\sMER(\d+)\s/){


    #read pairwise alingment Query line start and end sites, acid amino sequences and gap sites
    if(/^Query\s+?(\d+)\s+(.+?)\s+(\d+)/){
      $query_start = $1;
      $query = $2;
      $query_end = $3;
    }

    #read pairwise alingment Sbjct line start and end sites, acid amino sequences and gap sites
    if(/^Sbjct\s+?(\d+)\s+(.+?)\s+(\d+)/){
      $sbjct_start = $1 + $unit_start - 1;
      $sbjct = $2;
      $sbjct_end = $3 + $unit_start - 1;

      if($active_sites){
        $active_sbjct = '';
        $active_reslut_sbjct = '';
        $active_reslut_query = '';
        $active_query = '';
        $active_sites_query = '';

        #count actsite site in pairwise alingment
        if($active_sites == $sbjct_start){
          if($sbjct =~ /^(\w)/){
            $active_sbjct = $1;
            if($query =~ /^(\w)/){
              $active_query = $1;
              $active_sites_query = $query_start;
            }else{
              $active_query = 'MISS';
              $active_sites_query = 0;
            }
          }
          if($sbjct =~ /^(\-+)(\w)/){
            $length_active = length($1);
            $active_sbjct = $2;
            if($query =~ /^(.{$length_active})(\w)/){
              $active_query = $2;
              $query_tmp = $1;
              $count_tmp = ($query_tmp =~ s/\-/X/g);
              $active_sites_query = $query_start + $length_active unless $count_tmp;
              $active_sites_query = $query_start + $length_active - $count_tmp if $count_tmp;
            }else{
              $active_query = 'MISS';
              $active_sites_query = 0;
            }
          }
        }

        #count actsite site in pairwise alingment
        if($active_sites == $sbjct_end){
          if($sbjct =~ /(\w)$/){
            $active_sbjct = $1;
            if($query =~ /(\w)$/){
              $active_query = $1;
              $active_sites_query = $query_end;
            }else{
              $active_query = 'MISS';
              $active_sites_query = 0;
            }
          }
          if($sbjct =~ /(\w)(\-+)$/){
            $length_active = length($2);
            $active_sbjct = $1;
            if($query =~ /(\w)(.{$length_active})$/){
              $active_query = $1;
              $query_tmp = $2;
              $count_tmp = ($query_tmp =~ s/\-/X/g);
              $active_sites_query = $query_end - $length_active unless $count_tmp;
              $active_sites_query = $query_end - $length_active + $count_tmp if $count_tmp;
            }else{
              $active_query = 'MISS';
              $active_sites_query = 0;
            }
          }
        }

        #count actsite site in pairwise alingment
        if($active_sites > $sbjct_start and $active_sites < $sbjct_end){
          $active_count = $active_sites - $sbjct_start + 1;
          $active_count_copy = $active_count;
          $count1 = 0;
          my @array = split /\-+/, $sbjct;
          my @array1 = split /\w+/, $sbjct;
          my @array_copy = @array;
          foreach(@array){
            if($array_copy[0]){
              $count1++;
              my $count = length($_);
              if($active_count > $count){
                $active_count -= $count;
                next;
              }else{
                my $count_tmp = $active_count - 1;
                if(/.{$count_tmp}(\w)/){
                  $active_sbjct = $1;
                  last;
                }
              }
            }
            unless($array_copy[0]){
              next unless $_;
              $count1++;
              my $count = length($_);
              if($active_count > $count){
                $active_count -= $count;
                next;
              }else{
                my $count_tmp = $active_count - 1;
                if(/.{$count_tmp}(\w)/){
                  $active_sbjct = $1;
                  last;
                } 
              }
            }
          }
          $count_ = 0;
          foreach(@array1){
            last if $count1 == 0;
            $count1 -= 1;
            my $count = length($_);
            $count_ += $count;
          }
          $length_active = $active_count_copy + $count_ - 1;
          if($query =~ /^(.{$length_active})(\w)/){
            $active_query = $2;
            $query_tmp = $1;
            $count_tmp = ($query_tmp =~ s/\-/X/g);
            $active_sites_query = $query_start + $length_active unless $count_tmp;
            $active_sites_query = $query_start + $length_active - $count_tmp if $count_tmp;
          }else{
            $active_query = 'MISS';
            $active_sites_query = 0;
          }
        }
        #store actsites into array '@active_reslut_sbjct_array' and '@active_reslut_query_array'
        $active_reslut_sbjct = "$active_sbjct"."$active_sites" if $active_sbjct;
        @active_reslut_sbjct_array = ($active_reslut_sbjct) if $active_reslut_sbjct;
        $active_reslut_query = "$active_query"."$active_sites_query" if $active_query;
        @active_reslut_query_array = ($active_reslut_query) if $active_reslut_query;
      }
     
      if(@active_sites_array){
        #read each actsite in array '@active_sites_array'
        foreach(@active_sites_array){
          $active_sbjct = '';
          $active_reslut_sbjct = '';
          $active_reslut_query = '';
          $active_query = '';
          $active_sites_query = '';
          $active_sites_array = $_;
          next if $active_sites_array < $sbjct_start;
          last if $active_sites_array > $sbjct_end;
 
          #count actsite site in pairwise alingment
          if($active_sites_array == $sbjct_start){
            if($sbjct =~ /^(\w)/){
              $active_sbjct = $1;
              if($query =~ /^(\w)/){
                $active_query = $1;  
                $active_sites_query = $query_start;
              }else{
                $active_query = 'MISS';
                $active_sites_query = 0;
              }
            }
            if($sbjct =~ /^(\-+)(\w)/){
              $length_active = length($1);
              $active_sbjct = $2;
              if($query =~ /^(.{$length_active})(\w)/){
                $active_query = $2;
                $query_tmp = $1;
                $count_tmp = ($query_tmp =~ s/\-/X/g);
                $active_sites_query = $query_start + $length_active unless $count_tmp;
                $active_sites_query = $query_start + $length_active - $count_tmp if $count_tmp;
              }else{
                $active_query = 'MISS';
                $active_sites_query = 0;
              }
            }
          }

          #count actsite site in pairwise alingment
          if($active_sites_array == $sbjct_end){
            if($sbjct =~ /(\w)$/){
              $active_sbjct = $1;
              if($query =~ /(\w)$/){
                $active_query = $1;
                $active_sites_query = $query_end;
              }else{
                $active_query = 'MISS';
                $active_sites_query = 0;
              }
            }
            if($sbjct =~ /(\w)(\-+)$/){
              $length_active = length($2);
              $active_sbjct = $1;
              if($query =~ /(\w)(.{$length_active})$/){
                $active_query = $1;
                $query_tmp = $2;
                $count_tmp = ($query_tmp =~ s/\-/X/g);
                $active_sites_query = $query_end - $length_active unless $count_tmp;
                $active_sites_query = $query_end - $length_active + $count_tmp if $count_tmp;
              }else{
                $active_query = 'MISS';
                $active_sites_query = 0;
              }
            }
          }
  
          #count actsite site in pairwise alingment
          if($active_sites_array > $sbjct_start and $active_sites_array < $sbjct_end){
            $active_count = $active_sites_array - $sbjct_start + 1;
            $active_count_copy = $active_count;
            $count1 = 0;
            my @array = split /\-+/, $sbjct;
            my @array1 = split /\w+/, $sbjct;
            my @array_copy = @array;
            foreach(@array){
              if($array_copy[0]){
                $count1++;
                my $count = length($_);
                if($active_count > $count){
                  $active_count -= $count;
                  next;
                }else{
                  my $count_tmp = $active_count - 1;
                  if(/.{$count_tmp}(\w)/){
                    $active_sbjct = $1;
                    last;
                  }
                }
              }
              unless($array_copy[0]){
                next unless $_;
                $count1++;
                my $count = length($_);
                if($active_count > $count){
                  $active_count -= $count;
                  next;
                }else{
                  my $count_tmp = $active_count - 1;
                  if(/.{$count_tmp}(\w)/){
                    $active_sbjct = $1;
                    last;
                  }
                }
              }
            }
            $count_ = 0;
            foreach(@array1){
              last if $count1 == 0;
              $count1 -= 1;
              my $count = length($_);
              $count_ += $count;
            }
            $length_active = $active_count_copy + $count_ - 1;
            if($query =~ /^(.{$length_active})(\w)/){
              $active_query = $2;
              $query_tmp = $1;
              $count_tmp = ($query_tmp =~ s/\-/X/g);
              $active_sites_query = $query_start + $length_active unless $count_tmp;
              $active_sites_query = $query_start + $length_active - $count_tmp if $count_tmp;
            }else{
              $active_query = 'MISS';
              $active_sites_query = 0;
            }
          }

          #store actsites into array '@active_reslut_sbjct_array' and '@active_reslut_query_array'
          $active_reslut_sbjct = "$active_sbjct"."$active_sites_array" if $active_sbjct;
          push @active_reslut_sbjct_array, $active_reslut_sbjct if $active_reslut_sbjct;
          $active_reslut_query = "$active_query"."$active_sites_query";
          push @active_reslut_query_array, $active_reslut_query if $active_reslut_query;
        }
      }
    } #if(/^Sbjct\s+?(\d+)\s+(.+?)\s+(\d+)/){
  } #while(<IN>){

  close IN;

  #print out last actsites matching result
  if(@active_reslut_sbjct_array){print OUT "MER$id\tactive_site\t@active_reslut_sbjct_array\t";}
  unless(@active_reslut_sbjct_array){print OUT "MER$id\tactive_site_MISS\t$query_id\n" if $id;}
  if(@active_reslut_query_array){
    print OUT "$query_id\tactive_site\t@active_reslut_query_array\t";
    $switch_next = 1;
  }
  if(@active_reslut_sbjct_array or @active_reslut_query_array){print OUT "\n";}
  print STDOUT "Done.\n";

}

my $count_actsites_match;
my $database_blastp_all_query_id;
my $protease_blastp_all_query_id;
my $pepunit_blastp_all_query_id;
#Filter_actsite_matching_result(,'database' / 'protease' / 'pepunit')
sub Filter_actsite_matching_result{

  print STDOUT "Filter actsite matching result ...\n";
  $switch_next = 0;
  my $in1 = $_[0];
  my $key_tag = $_[1];
  my $in2;
  my $out1;
  my $out2;
  my $key_protease;
  if($key_tag eq 'database'){
    $in2 = $database_actsite_file;
    $out1 = $in1 . '_filter';
    $out2 = $in1 . '_all_query_id';
    $actsite_out_database_filter = $out1;
    $database_blastp_all_query_id = $out2;
    $key_protease = 2;
  }
  elsif($key_tag eq 'protease'){
    $in2 = $MEROPS_protease_actsites_all;
    $out1 = $in1 . '_filter';
    $out2 = $in1 . '_all_query_id';
    $actsite_out_protease_filter = $out1;
    $protease_blastp_all_query_id = $out2;
    $key_protease = 1;
  }
  elsif($key_tag eq 'pepunit'){
    $in2 = $MEROPS_protease_actsites_all;
    $out1 = $in1 . '_filter';
    $out2 = $in1 . '_all_query_id';
    $actsite_out_pepunit_filter = $out1;
    $pepunit_blastp_all_query_id = $out2;
    $key_protease = 2;
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  my $count;
  my %actsite_hash;
  my %mer_id_hash;
  open IN, $in2;
  while(<IN>){
    chomp;
    if(/^(\d+)\#\"(.+)\"\,/){
      my $uniport_id = $1;
      my $tmp = $2;
      $actsite_hash{$uniport_id} = $tmp;
    }
  }
  close IN;

  open IN, $in1;
  open OUT1, '>>', $out1;
  open OUT2, '>>', $out2;
  while(<IN>){
    chomp;
    if(/MISS/){
      next;
    }
    my $tmp_all = $_;
    if(/^MER(.+)\s+active\_site\s+(.+)\s+QER(\d+)\s+active\_site\s+(.+)\s$/){
      my $uniport_id = $1;
      next unless $actsite_hash{$uniport_id};
      my $tmp_sbjct = $2;
      my $mer_id = $3;
      my $tmp_query = $4;
      my %hash_tmp;
      my $count_next = 1;
      my $tmp_sbjct_chra = $tmp_sbjct;
      $tmp_sbjct_chra =~ s/\d+//g;
      $tmp_sbjct_chra =~ s/\s//g;
      my @actsite_sbjct_array = ($tmp_sbjct =~ /\d+/g);
      foreach my $actsite_tmp (@actsite_sbjct_array){
        $hash_tmp{$actsite_tmp} += 1;

        if($key_protease == 1){
          if($hash_tmp{$actsite_tmp} > 1){
            $count_next = 2;
            last;
          }
        }
      }
      if($key_protease == 1 && $count_next == 2){
        next;
      }
    
      @actsite_sbjct_array = keys %hash_tmp;
      @actsite_sbjct_array = sort {$a <=> $b} @actsite_sbjct_array;
      my $acsite_sbjct_count = scalar @actsite_sbjct_array;
      my $tmp_query_chra = $tmp_query;
      $tmp_query_chra =~ s/\d+//g;
      $tmp_query_chra =~ s/\s//g;
      my @actsite_query_array = ($tmp_query =~ /\d+/g);
      if($tmp_sbjct_chra ne $tmp_query_chra){
        next;
      }
      my %actsite_query_hash;
      foreach my $actsite_query (@actsite_query_array){
        $actsite_query_hash{$actsite_query} += 1;
        if($actsite_query_hash{$actsite_query} > 1){
          $count_next = 2;
          last;
        }
      }
      if($count_next == 2){
        next;
      }
      if($tmp_sbjct_chra eq $tmp_query_chra){
        my @acsite_uniport_array = ($actsite_hash{$uniport_id} =~ /\d+/g);
        my $acsite_uniport_count = scalar @acsite_uniport_array;
        my $count_tmp = 1;
        my %hash_tmp_tmp;
        if($acsite_sbjct_count > $acsite_uniport_count){
          next;
        }
        foreach my $actsite_tmp (@acsite_uniport_array){
          if($hash_tmp{$actsite_tmp}){
            $hash_tmp_tmp{$actsite_tmp} = 1;
          }else{
            $count_tmp = 2;
            last;
          }
        }
        if($count_tmp == 1){
          my @array_tmp = keys %hash_tmp_tmp;
          my $protease_actsite_count = scalar @array_tmp;
          if($protease_actsite_count == $acsite_uniport_count){
            print OUT1 "$tmp_all\n";
            $switch_next = 1;
            $mer_id_hash{$mer_id} += 1;
            print OUT2 "$mer_id\n" if $mer_id_hash{$mer_id} == 1;
          }
          if($protease_actsite_count < $acsite_uniport_count){
            next;
          }
        }
        if($count_tmp == 2){
          next;
        }
      }
    }
  }
  close IN;
  close OUT1;
  close OUT2;
  if($switch_next == 1){
    my @mer_id_count_array = keys %mer_id_hash;
    $count_actsites_match = scalar @mer_id_count_array;
  }
  unless($switch_next){
    $switch_next = 1010;
  }
  print STDOUT "Done.\n";
}

# Rearrange_actsite_matching_result(,'database'/'protease'/'pepunit');
sub Rearrange_actsite_matching_result{

  print STDOUT "Rearrange actsite matching result ...\n";
  my $in1 = $_[0];
  my $key_tag = $_[1];
  my $out1;

  if($in1){
    if(! -e $in1){
      print STDERR "ERROR: Filtered actsites mating result file does not exit. Please check.\n";
      exit(1);
    }   
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  if($key_tag eq 'database'){
    $out1 = $in1 . '_rearrange';
    $actsite_out_database_rearrange = $out1;
    }
  elsif($key_tag eq 'protease'){
    $out1 = $in1 . '_rearrange';
    $actsite_out_protease_rearrange = $out1;
    }
  elsif($key_tag eq 'pepunit'){
    $out1 = $in1 . '_rearrange';
    $actsite_out_pepunit_rearrange = $out1;
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  $switch_next = 0;
  my ($count);
  my %mer_id_hash;

  open IN, $in1;
  open OUT1, '>>', $out1;
  while(<IN>){
    chomp;
    my $tmp_all = $_;
    if(/^MER(.+)\s+active\_site\s+(.+)\s+QER(\d+)\s+active\_site\s+(.+)\s$/){
      my $uniport_id = $1;
      my $tmp_sbjct = $2;
      my $mer_id = $3;
      my $tmp_query = $4;

      my @actsite_sbjct_array = ($tmp_sbjct =~ /\w+/g);
      my $tmp_query_chra = $tmp_query;
      my @actsite_query_site_array = ($tmp_query =~ /\d+/g);
      my @actsite_query_full_array = ($tmp_query =~ /\w+/g);

      my %actsite_query_full_hash;
      my %actsite_query_order_hash;
      my %actsite_sbjct_full_hash;
      my $order_num;

      foreach my $actsite_query (@actsite_query_site_array){
        $order_num++;
        $actsite_query_order_hash{$actsite_query} = $order_num;
        my $full_query = shift @actsite_query_full_array;
        $actsite_query_full_hash{$actsite_query} = $full_query;
        my $full_sbjct = shift @actsite_sbjct_array;
        $actsite_sbjct_full_hash{$actsite_query} = $full_sbjct;
      }

      $order_num = 0;
      my $key_count = 1;
      my @actsite_query_site_array_sorted = sort{$a <=> $b} @actsite_query_site_array;
      foreach my $actsite_query (@actsite_query_site_array_sorted){
        $order_num++;
        if($actsite_query_order_hash{$actsite_query} != $order_num){
          $key_count = 2;
          last;
        }
      }
    
      if($key_count == 1){
        print OUT1 "$tmp_all\n";
        $switch_next = 1;
        next;
      }
      if($key_count == 2){
        my @sub_actsite_query_site_array;
        my @sub_actsite_sbjct_site_array;
        foreach my $actsite_query (@actsite_query_site_array_sorted){
          push @sub_actsite_query_site_array, $actsite_query_full_hash{$actsite_query};
          push @sub_actsite_sbjct_site_array, $actsite_sbjct_full_hash{$actsite_query};
        }
        print OUT1 "MER$uniport_id\tactive\_site\t@sub_actsite_sbjct_site_array\tQER$mer_id\tactive\_site\t@sub_actsite_query_site_array \n";
        $switch_next = 1;
      }
    }
  }
  close IN;
  close OUT1;
  print STDOUT "Done.\n";

}

#sub Delete_repeat_in_actsite_matching_result(, 'database'/'protease'/'pepunit')
sub Delete_repeat_in_actsite_matching_result{

  $switch_next = 0;
  print STDOUT "Delete repeat in actsite matching result ...\n";
  my $in = $_[0];
  my $key_tag = $_[1];
  my $out1;
  my ($count1,$count2);
  my %mer_id_hash;
  my %actsite_match_hash;
  my %mer_id_actsite_match_hash;

  if($in){
    if(! -e $in){
      print STDERR "ERROR: Rerranged actsites result file \'$in\' does not exit. Please check.\n";
      exit(1);
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  if($key_tag eq 'database'){
    $out1 = $in . '_delete_repeat';
    $actsite_out_database_delete_repeat = $out1;
  }
  elsif($key_tag eq 'protease'){
    $out1 = $in . '_delete_repeat';
    $actsite_out_protease_delete_repeat = $out1;
  }
  elsif($key_tag eq 'pepunit'){
    $out1 = $in . '_delete_repeat';
    $actsite_out_pepunit_delete_repeat = $out1;
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  open IN, $in;
  open OUT, '>>', $out1;
  while(<IN>){
    chomp;
    if(/active\_site.+\s+QER(\d+)\s+active\_site\s+(.+)\s/){
      my $mer_id = $1;
      my $actsite_match = $2;
      #count the repeat times of one mer_id from $out1
      $mer_id_hash{$mer_id} += 1;
      if($mer_id_hash{$mer_id} == 1){
        %actsite_match_hash = ();
        $count1++;
      }
      $actsite_match_hash{$actsite_match} += 1;
      if($actsite_match_hash{$actsite_match} == 1){
        print OUT "$_\n";
        $switch_next = 1;
        $count2++;
      }
    }else{
      print "$_";
    }
  }
  close IN;
  close OUT;
  print STDOUT "Done.\n";

}

# Separate_actsite_matching_result(,'database'/'protease'/'pepunit')
sub Separate_actsite_matching_result{

  print STDOUT "Separate actsite matching result ...\n";
  $switch_next = 0;
  my $in = $_[0];
  my $key_tag = $_[1];
  my $out1;
  my $out2;

  if($in){ 
    if(! -e $in){
      print STDERR "ERROR: actsites matching result file \'$in\' does not exit. Please check.\n";
      exit(1);
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  if($key_tag eq 'database'){
    $out1 = $in . '_result1';
    $out2 = $in . '_left1';
    $actsite_out_database_result1 = $out1;
    $actsite_out_database_left1 = $out2;
  }
  elsif($key_tag eq 'protease'){
    $out1 = $in . '_result1';
    $out2 = $in . '_left1';
    $actsite_out_protease_result1 = $out1;
    $actsite_out_protease_left1 = $out2;
  }
  elsif($key_tag eq 'pepunit'){
    $out1 = $in . '_result1';
    $out2 = $in . '_left1';
    $actsite_out_pepunit_result1 = $out1;
    $actsite_out_pepunit_left1 = $out2;
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }
    
  my ($count1,$count2);
  my %mer_id_hash;
  my %actsite_match_hash;
  my %mer_id_actsite_match_hash;

  open IN, $in;
  while(<IN>){
    chomp;
    if(/active\_site.+\s+QER(\d+)\s+active\_site\s+(.+)\s/){
      my $mer_id = $1;
      my $actsite_match = $2;
      #count the repeat times of one mer_id from $in
      $mer_id_hash{$mer_id} += 1; 
      if($mer_id_hash{$mer_id} == 1){
        %actsite_match_hash = ();
        $count1++;
      }
      $actsite_match_hash{$actsite_match} += 1;
      if($actsite_match_hash{$actsite_match} == 1){
        $count2++;
        if($mer_id_hash{$mer_id} == 1){
          $mer_id_actsite_match_hash{$mer_id} = $actsite_match;
        }
        if($mer_id_hash{$mer_id} > 1){
          my $actsite_match_tmp = ';' . $actsite_match;
          #collect all actsites within one mer_id 
          $mer_id_actsite_match_hash{$mer_id} .= $actsite_match_tmp;
        }
      }
    }
  }
  close IN;

  print STDOUT "$in mer_id sum: $count1 sum: $count2\n";
  %actsite_match_hash = ();
  $count1 = 0;
  $count2 = 0;
  open IN, $in;
  open OUT1, '>>', $out1;
  open OUT2, '>>', $out2;
  my %count_hash;
  my @count_array;
  my %actsite_full_len_hash;
  my %delete_count_hash_tmp;
  while(<IN>){
    chomp;
    my $tmp_all = $_;

    if(/active\_site.+\s+QER(\d+)\s+active\_site\s+(.+)\s/){
      my $mer_id = $1;
      my $actsite_match = $2;
      #count the repeat times of one mer_id from $in
      $count_hash{$mer_id} += 1;
      #if only one type of actsites with one mer_id, print the item to $out1;
      if($mer_id_hash{$mer_id} == 1){
        print OUT1 "$tmp_all\n";
        $actsite_result = 1;
      }

      #if only one type of actsites with one mer_id and the <IN> reading the first actsites \
      #collect all actsites by the sequence and collect the sequence of actsites;
      if($mer_id_hash{$mer_id} > 1 && $count_hash{$mer_id} == 1){
        $count1 = 0;
        %actsite_full_len_hash = ();
        %delete_count_hash_tmp = ();
        my @actsite_array_tmp = split /\;/, $mer_id_actsite_match_hash{$mer_id};
        foreach my $actsite_full_len (@actsite_array_tmp){
          #if more than one type of actsites with one mer_id, count the sequence of actsites;
          $count1++;
          #collect all actsites within one mer_id separately by the sequence of actsites;
          $actsite_full_len_hash{$count1} = $actsite_full_len;
        }
        #collect the sequence of actsites
        @count_array = keys %actsite_full_len_hash;
        @count_array = sort{$a <=> $b} @count_array;
      }

      #deal with every actsites belongs to one type of actsites with one mer_id one by one;
      if($mer_id_hash{$mer_id} > 1){
        my $num = scalar @count_array;
        if($num == 1){
          if($count_array[0] == $count_hash{$mer_id}){
            print OUT1 "$tmp_all\n";
            $actsite_result = 1;
          }
          next;
        }
        next if $delete_count_hash_tmp{$count_hash{$mer_id}};
        my $count_tmp;
        #fetch actsite of one type actsites;
        my @actsite_match_array = ($actsite_match =~ /\d+/g);
        my %actsite_tmp_hash;
        #count the number of actsite in one type actsites;
        my $actsite_match_array_num = scalar @actsite_match_array;
        #collect actsite of one type actsites by hash for inclusion or being inclusion relationship \
        #checking and check repeat of every actsite
        foreach(@actsite_match_array){
          $actsite_tmp_hash{$_} += 1;
          if($actsite_tmp_hash{$_} > 1){
            print "ERROR at $tmp_all\n";
          }
        }

        my @delete_count_array;
        foreach my $mer_id_count (@count_array){
          my $actsite_num;
          my $type = 1;
          next if $mer_id_count == $count_hash{$mer_id};
          #fetch actsite of one type actsites except the current sequence type of actsites
          my @actsite_full_len_array = ($actsite_full_len_hash{$mer_id_count} =~ /\d+/g);
          #check inclusion or being inclusion relationship between current sequence type of actsites and others    
          foreach my $actsite_full_len_tmp (@actsite_full_len_array){
            #check whether the actsite is repeat in current sequence type of actsites or not
            if($actsite_tmp_hash{$actsite_full_len_tmp}){
              $actsite_num++;
            }
            unless($actsite_tmp_hash{$actsite_full_len_tmp}){
              $type = 2;
            }
          }

          if($type == 1){
            #the actsites is included by current sequence type of actsites
            if($actsite_num < $actsite_match_array_num){
              push @delete_count_array, $mer_id_count;
              next;
            }
            #the actsites is equal to current sequence type of actsites
            if($actsite_num == $actsite_match_array_num){
              push @delete_count_array, $mer_id_count;
              next;
            }
          }

          if($type == 2){
            #the actsites is completely different from current sequence type of actsites
            unless($actsite_num){
              next;
            }
            #the actsites include current sequence type of actsites
            if($actsite_num == $actsite_match_array_num){
              push @delete_count_array, $count_hash{$mer_id};
              last;
            }
            #the actsites is different from current sequence type of actsites
            if($actsite_num < $actsite_match_array_num){
              next;
            }
          }
        }
        my %delete_count_hash = ();
        if(@delete_count_array){
          foreach(@delete_count_array){
            $delete_count_hash{$_} += 1;
            if($delete_count_hash{$_} > 1){
              print "ERROR $tmp_all\n";
            }
          }
          my @count_array_tmp = @count_array;
          @count_array = ();
          foreach(@count_array_tmp){
            unless($delete_count_hash{$_}){
              push @count_array, $_;
            }
          }
        }
        %delete_count_hash_tmp = %delete_count_hash;
        unless(%delete_count_hash){
          if(scalar(@count_array) == 1){
            if($count_array[0] == $count_hash{$mer_id}){
              print OUT1 "$tmp_all\n";
              $actsite_result = 1;
            }
          }
          if(scalar(@count_array) > 1){
            print OUT2 "$tmp_all\n";
            $switch_next = 1;
          }
        }
        if(%delete_count_hash && !$delete_count_hash{$count_hash{$mer_id}}){
          if(scalar(@count_array) == 1){
            if($count_array[0] == $count_hash{$mer_id}){
              print OUT1 "$tmp_all\n";
              $actsite_result = 1;
            }
          }
          if(scalar(@count_array) > 1){
            print OUT2 "$tmp_all\n";
            $switch_next = 1;
          }
        }
      } 
    }
  }
  close IN;
  close OUT1;
  close OUT2;

  if($key_tag eq 'database'){
    $actsite_result_database = $actsite_result; 
  }
  elsif($key_tag eq 'protease'){
    $actsite_result_protease = $actsite_result;
  }
  elsif($key_tag eq 'pepunit'){
    $actsite_result_pepunit = $actsite_result;
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1); 
  }
  print STDOUT "Done.\n";

}

# Separate_AA_actsite_matching_result(,'database'/'protease'/'pepunit')
sub Separate_AA_actsite_matching_result{

  print STDOUT "Separate AA_actsite matching result ...\n";
  $switch_next = 0;
  my $in = $_[0];
  my $key_tag = $_[1];
  my $actsite_result_next;
  my $out1;
  my $out2;
 
  if($in){
    if(! -e $in){
      print STDERR "ERROR: actsites matching result file \'$in\' does not exit. Please check.\n";
      exit(1);
    }
  }
  unless($in){
    print STDERR "FATAL ERROR1: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  if($key_tag eq 'database'){
    $out1 = $in . '_result2';
    $out2 = $in . '_left2';
    $actsite_out_database_result2 = $out1;
    $actsite_out_database_left2 = $out2;
    if($actsite_result_database){
      $actsite_result = $actsite_result_database;
      if($actsite_result_database == 1){
        $actsite_result_next = 2;
      }else{
        print STDERR "FATAL ERROR2: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
        exit(1);
      }
    }
    unless($actsite_result_database){
      $actsite_result_next = 21;
    }
  }
  elsif($key_tag eq 'protease'){
    $out1 = $in . '_result2';
    $out2 = $in . '_left2';
    $actsite_out_protease_result2 = $out1;
    $actsite_out_protease_left2 = $out2;
    if($actsite_result_protease){
      $actsite_result = $actsite_result_protease;
      if($actsite_result_protease == 1){
        $actsite_result_next = 2;
      }else{
        print STDERR "FATAL ERROR3: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
        exit(1);
      }
    }
    unless($actsite_result_protease){
      $actsite_result_next = 21;
    }
  }
  elsif($key_tag eq 'pepunit'){
    $out1 = $in . '_result2';
    $out2 = $in . '_left2';
    $actsite_out_pepunit_result2 = $out1;
    $actsite_out_pepunit_left2 = $out2;
    if($actsite_result_pepunit){
      $actsite_result = $actsite_result_pepunit;
      if($actsite_result_pepunit == 1){
        $actsite_result_next = 2;
      }else{
        print STDERR "FATAL ERROR4: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
        exit(1);
      }
    }
    unless($actsite_result_pepunit){
      $actsite_result_next = 21;
    }
  }else{
    print STDERR "FATAL ERROR5: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  my %mer_id_hash;
  my %mer_id_actsites_hash;
  my %actsite_match_hash;
  my ($count1, $count2);
  open IN, $in;
  while(<IN>){
    chomp;
    if(/active\_site.+\s+QER(\d+)\s+active\_site\s+(.+)\s/){
      my $mer_id = $1;
      my $actsite_match = $2;
      #count the repeat times of one mer_id from $out1
      $mer_id_hash{$mer_id} += 1;
      if($mer_id_hash{$mer_id} == 1){
        %actsite_match_hash = ();
      }
      $actsite_match =~ s/\d+//g;
      $actsite_match =~ s/\s+//g;
      $actsite_match_hash{$actsite_match} += 1;
      if($actsite_match_hash{$actsite_match} == 1){
        $mer_id_actsites_hash{$mer_id} += 1;
      }
    }
  }
  close IN;

  %mer_id_hash = ();
  open IN, $in;
  open OUT1, '>>', $out1;
  open OUT2, '>>', $out2;
  my %mer_id_poly_type_hash;
  while(<IN>){
    chomp;
    my $tmp_all = $_;
    if(/active\_site.+\s+QER(\d+)\s+active\_site/){
      my $mer_id = $1;
      $mer_id_hash{$mer_id} += 1;
      if($mer_id_actsites_hash{$mer_id} == 1 && $mer_id_hash{$mer_id} == 1){
        print OUT1 "$tmp_all\n";
        $actsite_result = $actsite_result_next;
        $count1++;
      }
      if($mer_id_actsites_hash{$mer_id} == 1 && $mer_id_hash{$mer_id} > 1){
        next;
      }
      if($mer_id_actsites_hash{$mer_id} > 1){
        print OUT2 "$tmp_all\n";
        $switch_next = 1;
        $mer_id_poly_type_hash{$mer_id} = 1;
      }
    }
  }
  close IN;
  close OUT1;
  close OUT2;


  if($key_tag eq 'database'){
    $actsite_result_database = $actsite_result;
    }
    elsif($key_tag eq 'protease'){
      $actsite_result_protease = $actsite_result;
    }
    elsif($key_tag eq 'pepunit'){
      $actsite_result_pepunit = $actsite_result;
    }else{
      print STDERR "FATAL ERROR6: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
      exit(1);
  }
  print STDOUT "Done.\n";

}

# Delete_inclusion_actsite_matching_result(,'database'/'protease'/'pepunit')
sub Delete_inclusion_actsite_matching_result{

  print STDOUT "Delete inclusion actsite matching result ...\n";
  my $in = $_[0];
  my $key_tag = $_[1];
  my $actsite_result_next;
  my $out1;
  my $out2;
  my $out3;

  if($in){
    if(! -e $in){
      print STDERR "ERROR: actsites matching result file \'$in\' does not exit. Please check.\n";
      exit(1);
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  if($key_tag eq 'database'){
    $out1 = $in . '_tmp';
    $out3 = $in . '_result3';
    $out2 = $in . '_proof';
    $actsite_out_database_result3 = $out3;
    $actsite_out_database_proof = $out2;
    if($actsite_result_database){
      $actsite_result = $actsite_result_database;
      if($actsite_result_database == 1){
        $actsite_result_next = 32;
      }
      elsif($actsite_result_database == 2){
       $actsite_result_next = 3;
      }
      elsif($actsite_result_database == 21){
       $actsite_result_next = 33;
      }
      else{
        print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
        exit(1);
      }
    }
    unless($actsite_result_database){
      $actsite_result_next = 31;
    }
  }
  elsif($key_tag eq 'protease'){
    $out1 = $in . '_tmp';
    $out3 = $in . '_result3';
    $out2 = $in . '_proof';
    $actsite_out_protease_result3 = $out3;
    $actsite_out_protease_proof = $out2;
    if($actsite_result_protease){
      $actsite_result = $actsite_result_protease;
      if($actsite_result_protease == 1){
        $actsite_result_next = 32;
      }
      elsif($actsite_result_protease == 2){
        $actsite_result_next = 3;
      }
      elsif($actsite_result_protease == 21){
        $actsite_result_next = 33;
      }
      else{
        print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
        exit(1);
      }
    }
    unless($actsite_result_protease){
      $actsite_result_next = 31;
    }
  }
  elsif($key_tag eq 'pepunit'){
    $out1 = $in . '_tmp';
    $out3 = $in . '_result3';
    $out2 = $in . '_proof';
    $actsite_out_pepunit_result3 = $out3;
    $actsite_out_pepunit_proof = $out2;
    if($actsite_result_pepunit){
      $actsite_result = $actsite_result_pepunit;
      if($actsite_result_pepunit == 1){
        $actsite_result_next = 32;
      }
      elsif($actsite_result_pepunit == 2){
        $actsite_result_next = 3;
      }
      elsif($actsite_result_pepunit == 21){
        $actsite_result_next = 33;
      }
      else{
        print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
        exit(1);
      }
    }
    unless($actsite_result_pepunit){
      $actsite_result_next = 31;
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  my %mer_id_hash;
  my %actsites_type_hash;
  my $count_sum;
  my $count_all;
  open IN, $in;
  open OUT, '>>', $out1;
  while(<IN>){
    chomp;
    my $tmp_all = $_;
    if(/^(.+active\_site.+\s+QER(\d+)\s+active\_site\s+)(.+)\s/){
      $count_sum++;
      my $tmp_part = $1;
      my $mer_id = $2;
      my $actsite_match = $3;
      $mer_id_hash{$mer_id} += 1;
      my @actsites_unsorted_array = ($actsite_match =~ /\d+/g);
      my @actsites_sorted_array = sort{$a <=> $b} @actsites_unsorted_array;
      my $count = 0;
      my $count1 = 1;
      foreach my $actsites_tmp (@actsites_unsorted_array){
        if($actsites_tmp != $actsites_sorted_array[$count]){
          $count1 = 2;
          last;
        }
        $count++;
      }
      if($count1 == 1){
        print OUT "$tmp_all\n";
        $count_all++;
      }
      $count = 0;
      my %actsites_all_hash;
      if($count1 == 2){
        my @actsites_all_array = ($actsite_match =~ /\w+/g);
        foreach my $actsite_all (@actsites_all_array){
          $actsites_all_hash{$actsites_unsorted_array[$count]} = $actsite_all;
          $count++;
        }
          @actsites_unsorted_array = ();
        foreach my $actsite_all (@actsites_sorted_array){
          push @actsites_unsorted_array, $actsites_all_hash{$actsite_all};
        }
        print OUT "$tmp_part@actsites_unsorted_array\t\n";
        $count_all++;
      }
    }
  }
  close IN;
  close OUT;

  %mer_id_hash = ();
  open IN, $out1;
  open OUT, '>>', $out2;
  while(<IN>){
    chomp;
    my $tmp_all = $_;
    if(/active\_site.+\s+QER(\d+)\s+active\_site\s+(.+)\s/){
      my $mer_id = $1;
      my $actsite_match = $2;
      $mer_id_hash{$mer_id} += 1;
      if($mer_id_hash{$mer_id} == 1){
        %actsites_type_hash = ();
      }
      $actsite_match =~ s/\d+//g;
      $actsite_match =~ s/\s+//g;
      $actsites_type_hash{$actsite_match} += 1;
      if($actsites_type_hash{$actsite_match} == 1){
        print OUT "$tmp_all\n";
      } 
    }
  }
  close IN;
  close OUT;
  `rm $out1`;
  %mer_id_hash = ();
  %actsites_type_hash = ();

  open IN, $out2;
  while(<IN>){
    chomp;
    if(/active\_site.+\s+QER(\d+)\s+active\_site\s+(.+)\s/){
      my $mer_id = $1;
      my $actsite_match = $2;
      $mer_id_hash{$mer_id} += 1;
      $actsite_match =~ s/\d+//g;
      $actsite_match =~ s/\s+//g;
      if($mer_id_hash{$mer_id} == 1){
        $actsites_type_hash{$mer_id} = $actsite_match;
      }
      if($mer_id_hash{$mer_id} > 1){
        my $actsite_match_tmp = ';' . $actsite_match;
        $actsites_type_hash{$mer_id} .= $actsite_match_tmp;
      }
    }
  }
  close IN;
  my @actsites_type_array;
  %mer_id_hash = ();
  open IN, $out2;
  open OUT, '>>', $out3;
  while(<IN>){
    chomp;
    my $tmp_all = $_;
    my $count_actsites_num;
    if(/active\_site.+\s+QER(\d+)\s+active\_site\s+(.+)\s/){
      my $mer_id = $1;
      my $actsite_match = $2;
      $actsite_match =~ s/\d+/.*/g;
      $actsite_match =~ s/\s+//g;
      $mer_id_hash{$mer_id} += 1;
      if($mer_id_hash{$mer_id} == 1){
        @actsites_type_array = split /\;/, $actsites_type_hash{$mer_id};
      }
      my $count_tmp = 1;
      foreach my $actsites_type (@actsites_type_array){
        $count_actsites_num++;
        if($mer_id_hash{$mer_id} == $count_actsites_num){
          next;
        }
        if($actsites_type =~ /$actsite_match/){
          $count_tmp = 2;
          last;
        }
      }
      if($count_tmp == 1){
        print OUT "$tmp_all\n";
        $actsite_result = $actsite_result_next;
      }
    }
  }
  close IN;
  close OUT;
  
  if($key_tag eq 'database'){
    $actsite_result_database = $actsite_result;
  }
  elsif($key_tag eq 'protease'){
    $actsite_result_protease = $actsite_result;
  }
  elsif($key_tag eq 'pepunit'){
    $actsite_result_pepunit = $actsite_result;
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }
  print STDOUT "Done.\n";
  
}


# combine_actsites('database'/'protease'/'pepunit')
sub combine_actsites{

  print STDOUT "Combine actsites results...\n";
  $switch_next = 0;
  my $key_tag = $_[0];
  my $cmd_cat;
  if($actsite_result){
    if($key_tag eq 'database'){
      $actsites_annotation_database = $query_input_tmp . '_actsite_annotation_database';
    }
    elsif($key_tag eq 'protease'){
      $actsites_annotation_protease = $query_input_tmp . '_actsite_annotation_protease';
    }
    elsif($key_tag eq 'pepunit'){
      $actsites_annotation_pepunit = $query_input_tmp . '_actsite_annotation_pepunit';
    }
    else{
      print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
      exit(1);
    }
    if($actsite_result == 1){
      if($key_tag eq 'database'){
        $cmd_cat = "cat $actsite_out_database_result1 > $actsites_annotation_database";
      }
      elsif($key_tag eq 'protease'){
        $cmd_cat = "cat $actsite_out_protease_result1 > $actsites_annotation_protease";
      }
      elsif($key_tag eq 'pepunit'){
        $cmd_cat = "cat $actsite_out_pepunit_result1 > $actsites_annotation_pepunit"; 
      }
      else{
        print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
        exit(1);
      }
      system("$cmd_cat")==0 or die("failed to execute command \'$cmd_cat\': $!\n");
      $switch_next = 1;
    }
    elsif($actsite_result == 2){
      if($key_tag eq 'database'){
        $cmd_cat = "cat $actsite_out_database_result1 $actsite_out_database_result2 > $actsites_annotation_database";
      }
      elsif($key_tag eq 'protease'){
        $cmd_cat = "cat $actsite_out_protease_result1 $actsite_out_protease_result2 > $actsites_annotation_protease";
      }
      elsif($key_tag eq 'pepunit'){
        $cmd_cat = "cat $actsite_out_pepunit_result1 $actsite_out_pepunit_result2 > $actsites_annotation_pepunit";
      }
      else{
        print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
        exit(1);
      }
      system("$cmd_cat")==0 or die("failed to execute command \'$cmd_cat\': $!\n");
      $switch_next = 1;
    }
    elsif($actsite_result == 21){
      if($key_tag eq 'database'){
        $cmd_cat = "cat $actsite_out_database_result2 > $actsites_annotation_database";
      }
      elsif($key_tag eq 'protease'){
        $cmd_cat = "cat $actsite_out_protease_result2 > $actsites_annotation_protease";
      }
      elsif($key_tag eq 'pepunit'){
        $cmd_cat = "cat $actsite_out_pepunit_result2 > $actsites_annotation_pepunit";
      }
      else{
        print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
        exit(1);
      }
      system("$cmd_cat")==0 or die("failed to execute command \'$cmd_cat\': $!\n");
      $switch_next = 1;
    }
    elsif($actsite_result == 3){
      if($key_tag eq 'database'){
        $cmd_cat = "cat $actsite_out_database_result1 $actsite_out_database_result2 $actsite_out_database_result3 > $actsites_annotation_database";
      }
      elsif($key_tag eq 'protease'){
        $cmd_cat = "cat $actsite_out_protease_result1 $actsite_out_protease_result2 $actsite_out_protease_result3 > $actsites_annotation_protease";
      }
      elsif($key_tag eq 'pepunit'){
        $cmd_cat = "cat $actsite_out_pepunit_result1 $actsite_out_pepunit_result2 $actsite_out_pepunit_result3 > $actsites_annotation_pepunit";
      }
      else{
        print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
        exit(1);
      }
      system("$cmd_cat")==0 or die("failed to execute command \'$cmd_cat\': $!\n");
      $switch_next = 1;
    }
    elsif($actsite_result == 31){
       if($key_tag eq 'database'){
        $cmd_cat = "cat $actsite_out_database_result3 > $actsites_annotation_database";
      }
      elsif($key_tag eq 'protease'){
        $cmd_cat = "cat $actsite_out_protease_result3 > $actsites_annotation_protease";
      }
      elsif($key_tag eq 'pepunit'){
        $cmd_cat = "cat $actsite_out_pepunit_result3 > $actsites_annotation_pepunit";
      }
      else{
        print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
        exit(1);
      }
      system("$cmd_cat")==0 or die("failed to execute command \'$cmd_cat\': $!\n");
      $switch_next = 1;
    }
    elsif($actsite_result == 32){
       if($key_tag eq 'database'){
        $cmd_cat = "cat $actsite_out_database_result1 $actsite_out_database_result3 > $actsites_annotation_database";
      }
      elsif($key_tag eq 'protease'){
        $cmd_cat = "cat $actsite_out_protease_result1 $actsite_out_protease_result3 > $actsites_annotation_protease";
      }
      elsif($key_tag eq 'pepunit'){
        $cmd_cat = "cat $actsite_out_pepunit_result1 $actsite_out_pepunit_result3 > $actsites_annotation_pepunit";
      }
      else{
        print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' a
nd replace the modified 'PepTidy.pl'.\n";
        exit(1);
      }
      system("$cmd_cat")==0 or die("failed to execute command \'$cmd_cat\': $!\n");
      $switch_next = 1;
    }
    elsif($actsite_result == 33){
       if($key_tag eq 'database'){
        $cmd_cat = "cat $actsite_out_database_result2 $actsite_out_database_result3 > $actsites_annotation_database";
      }
      elsif($key_tag eq 'protease'){
        $cmd_cat = "cat $actsite_out_protease_result3 $actsite_out_protease_result3 > $actsites_annotation_protease";
      }
      elsif($key_tag eq 'pepunit'){
        $cmd_cat = "cat $actsite_out_pepunit_result3 $actsite_out_pepunit_result3 > $actsites_annotation_pepunit";
      }
      else{
        print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
        exit(1);
      }
      system("$cmd_cat")==0 or die("failed to execute command \'$cmd_cat\': $!\n");
      $switch_next = 1;
    }else{
      print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    }
 
  }else{
    print STDOUT "No actsites annotation result.\n";
  }
  print STDOUT "Done.\n";
}

my $balstp_actsites_left_seq;
#if($count_hits_match > $count_actsites_match)
#Fetch_actsites_sequences($protease_blastp_all_query_id ,$blastp_pepunit_out_seq)
sub Fetch_actsites_sequences{

    print STDOUT "Fetch unmatched sequences ...\n";
    $switch_next = 0;
    my $in1 = $_[0];
    my $in2 = $_[1];
    my $out = $query_input_tmp . '_left_seq';
    $balstp_actsites_left_seq = $out;
    my %query_id_hash;
    open IN, $in1;
    while(<IN>){
      chomp;
      if(/^(\d+)/){
        my $query_id = $1;
        $query_id_hash{$query_id} = 1;
      }
    }
    close IN;
   
    open IN, $in2;
    open OUT, '>>', $out;
    my $count_print = 1;
    while(<IN>){
      chomp;
      my $tmp_all = $_;
      if(/^>QER(\d+)/){
        my $query_id = $1;
        if($query_id_hash{$query_id}){
          $count_print = 1;
        }
        unless($query_id_hash{$query_id}){
          $count_print = 2;
        }
      }
      if($count_print == 2){
        print OUT "$tmp_all\n";
        $switch_next = 1;
      }
    }
    print "Done.\n";
}


#my %actsite_AMB_result;
sub Transform_actsites{

  my %uniport_id_hash;
  if($db_build){
    my $out1 = $workDir  . '/' . $db_name  . '_seq_actsites_all';
    my $out2 = $workDir  . '/' . $db_name  . '_seq_actsites_AMB_sum';
    my $out3 = $workDir  . '/' . $db_name  . '_seq_uniport_annotation_sum';
    open OUT1, '>>', $out1;
    open OUT2, '>>', $out2;
    open OUT3, '>>', $out3;

    open IN, $mer_uniport_id_sum;
    while(<IN>){
      chomp;
      if(/^(\d+)\#(.+)$/){
        my $mer_id = $1;
        my $uniport_id = $2;
        $uniport_id_hash{$mer_id} = $uniport_id;
      }
    }
    close IN;
  }
  #$actsite_AMB_result = $out;
  print STDOUT "Transform actsites ...\n";
  $actsites_filter_all = $query_input_tmp . '_actsite_filter_all';
  if($actsite_protease_key){
    if(! -e $actsite_out_protease_filter){
      print STDERR "ERROR: no \'$actsite_out_protease_filter\' file.\n";
      exit(1);
    }
    if($actsite_pepunit_key){
      if(! -e $actsite_out_pepunit_filter){
        print STDERR "ERROR: no \'$actsite_out_pepunit_filter\' file.\n";
        exit(1);
      }
      my $cmd_combine_all = "cat $actsite_out_protease_filter $actsite_out_pepunit_filter > $actsites_filter_all";
      system("$cmd_combine_all")==0 or die("failed to execute command \'$cmd_combine_all\': $!\n");
    }
    unless($actsite_pepunit_key){
      $actsites_filter_all = $actsite_out_protease_filter;
    }
  }
  unless($actsite_protease_key){
    if(! -e $actsite_out_pepunit_filter){
      print STDERR "ERROR: no \'$actsite_out_pepunit_filter\' file.\n";
      exit(1);
    }
    $actsites_filter_all = $actsite_out_pepunit_filter;
  }

  my %query_mer_id_hash;
  my %mer_id_hash;
  my %query_id_build_hash;
  my %query_mer_id_build_hash;
  open IN, $actsites_annotation_all;
  while(<IN>){
    chomp;
    if(/^MER(\d+).+\s(QER(\d+))\tactive_site\t(.+)\s$/){
      my $mer_id = $1;
      my $query_id = $2;
      my $query_id_num = $3;
      my $acsites = $4;
      $mer_id_hash{$mer_id} += 1;
      $query_id_build_hash{$query_id} += 1; 
      my $mer_query_id = 'MER' . $mer_id . $query_id;
      $query_mer_id_hash{$mer_query_id} += 1;
      if($query_mer_id_hash{$mer_query_id} > 1){
        print STDERR "ERROR: $query_id MER$mer_id repeat in actsite annotation result.\nThe primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
        exit(1);
      }
      if($db_build){
        if($merops_switch){
          if($query_ori_id_hash{$query_id}){
            next;
          }
        }
        if($query_id_build_hash{$query_id} == 1){
          $query_id_annotation_seq_hash{$query_id_num} += 1;
          $query_mer_id_build_hash{$mer_query_id} += 1;
          my $count_acsite_tmp;
          my @actsites_array = ($acsites =~ /\w+/g);
          my $count_array = scalar @actsites_array;
          print OUT1 "$query_id_num\#\"";
          foreach my $acsite_tmp (@actsites_array){
            $count_acsite_tmp++;
            if($count_acsite_tmp == $count_array){
              print OUT1 "$acsite_tmp\"\,\n";
              last;
            }
            if($count_acsite_tmp < $count_array){
              print OUT1 "$acsite_tmp\, ";
            }
          }
          print OUT3 "$query_id_num\#$uniport_id_hash{$mer_id}\n";
        }
      }
    }
  }
  close IN;


  #my %actsite_AMB_query_mer_id_hash;
  open IN, $actsites_filter_all;
  #open OUT, '>>', $out;
  while(<IN>){
    chomp;
    if(/^MER(\d+)\tactive_site\t(.+)\t(QER(\d+))\tactive_site\t(.+)\t/){
      my $mer_id = $1;
      my $query_id = $3;
      my $actsite_mer = $2;
      my $query_id_num = $4;
      my $actsite_qer = $5;
      my $mer_query_id = 'MER' . $mer_id . $query_id;
      if($query_mer_id_hash{$mer_query_id}){
        my @actsite_mer_array = split /\s/, $actsite_mer;
        my @actsite_qer_array = split /\s/, $actsite_qer;
        my @actsite_AMB_mer_cha_array = ($actsite_AMB_mer_id_hash{$mer_id} =~ /[AMBJOXZ]/g);
        my @actsite_AMB_mer_num_array = ($actsite_AMB_mer_id_hash{$mer_id} =~ /\d+/g);
        unless(@actsite_AMB_mer_cha_array || @actsite_AMB_mer_num_array){
          print STDERR "ERROR: Wrong format at MER$mer_id in file \'$active_metal_binding_AMB_sum\'.\nPlease redownload the file \'$active_metal_binding_AMB_sum\' and replace it.\n";
          exit(1);
        }
        my $count_cha = scalar @actsite_AMB_mer_cha_array;
        my $count_num = scalar @actsite_AMB_mer_num_array;
        if($count_cha != $count_num){
          print STDERR "ERROR: Wrong format at MER$mer_id in file \'$active_metal_binding_AMB_sum\'.\nPlease redownload the file \'$active_metal_binding_AMB_sum\' and replace it.\n";
          exit(1);
        }
        my %actsite_AMB_mer_hash;
        my $iii;
        for($iii = 0; $iii < $count_cha; $iii++){
           $actsite_AMB_mer_hash{$actsite_AMB_mer_num_array[$iii]} = $actsite_AMB_mer_cha_array[$iii];
        }
        my %actsite_AMB_mer_seq_hash;
        my $count_seq = 1;
        foreach my $actsite_AMB_mer (@actsite_mer_array){
          $actsite_AMB_mer =~ s/[^0-9]//g;
          $actsite_AMB_mer_seq_hash{$count_seq} = $actsite_AMB_mer_hash{$actsite_AMB_mer};
          $count_seq++;
        }
        $count_num++;
        if($count_num != $count_seq){
          print STDOUT "Repeat domain: $mer_query_id MER_actsite: @actsite_mer_array MER_AMB: $actsite_AMB_mer_id_hash{$mer_id}\n";
        }
        $count_seq = 1;
        foreach my $actsite_qer (@actsite_qer_array){
          $actsite_qer =~ s/[^0-9]//g;
          my $actsite_AMB_qer = $actsite_AMB_mer_seq_hash{$count_seq} . $actsite_qer . ' ';
          $actsite_AMB_result{$mer_query_id} .= $actsite_AMB_qer;
          $count_seq++;
        }
        if($db_build){
          if($merops_switch){
            if($query_ori_id_hash{$query_id}){
              next;
            }
          }
          if($query_mer_id_build_hash{$mer_query_id}){
            my $AMB_actsites = $actsite_AMB_result{$mer_query_id};
            my @AMB_actsites_array = ($AMB_actsites =~ /\w+/g);
            my $count_array = scalar @AMB_actsites_array;
            my $count_acsite_tmp;
            print OUT2 "$query_id_num\#\"";
            foreach my $AMB_acsite_tmp (@AMB_actsites_array){
              $count_acsite_tmp++;
              if($count_acsite_tmp == $count_array){
                print OUT2 "$AMB_acsite_tmp\"\,\n";
                last;
              }
              if($count_acsite_tmp < $count_array){
                print OUT2 "$AMB_acsite_tmp\, ";
              }
            }  
          }
        }
      }
    }
  }
  close IN;

  if($db_build){
    close OUT1;
    close OUT2;
    close OUT3;
  }

}


#create_merops_custom_database($protease_same_query);
sub create_merops_custom_database{
  my $in = $_[0];
  my $out1 = $workDir  . '/' . $db_name  . '_seq_actsites_all';
  my $out2 = $workDir  . '/' . $db_name  . '_seq_actsites_AMB_sum';
  my $out3 = $workDir  . '/' . $db_name  . '_seq_uniport_annotation_sum';
  open OUT1, '>>', $out1;
  open OUT2, '>>', $out2;
  open OUT3, '>>', $out3;

  my %uniport_id_hash;
  open IN, $mer_uniport_id_sum;
  while(<IN>){
    chomp;
    if(/^(\d+)\#(.+)$/){
      my $mer_id = $1;
      my $uniport_id = $2;
      $uniport_id_hash{$mer_id} = $uniport_id;
    }
  }
  close IN;

  my %actsite_hash;
  open IN, $MEROPS_protease_actsites_all;
  while(<IN>){
    chomp;
    if(/^(\d+)\#(\".+\"\,)/){
      my $mer_id = $1;
      my $tmp = $2;
      $actsite_hash{$mer_id} = $tmp;
    }
  }
  close IN;

  open IN, $in;
  while(<IN>){
    chomp;
    if(/^(QER(\d+))\#MER(\d+)$/){
      my $query_id = $1;
      my $query_id_num = $2;
      my $mer_id = $3;
      $query_ori_id_hash{$query_id} += 1;
      $query_id_annotation_seq_hash{$query_id_num} += 1;
      my $uniport_id = $uniport_id_hash{$mer_id};
      my $AMB_actsites = $actsite_AMB_mer_id_hash{$mer_id};
      my $actsites = $actsite_hash{$mer_id};
      print OUT1 "$query_id_num\#$actsites\n";
      print OUT2 "$query_id_num\#\"$AMB_actsites\"\,\n";
      print OUT3 "$query_id_num\#$uniport_id\n";
    }
  }
  close IN;
  close OUT1;
  close OUT2;
  close OUT3;

}







       ################## subfunction for Domain Annotation ################

# Extract_mer_seq($actsites_annotation_database, $actsite_out_database_rearrange, $database_MEROPS_pepunit) 
# Extract_mer_seq($actsites_annotation_protease, $actsite_out_protease_rearrange ,$database_MEROPS_pepunit)
# Extract_mer_seq($actsites_annotation_pepunit, $actsite_out_pepunit_rearrange ,$database_MEROPS_pepunit)
my $domain_annotation_pepuint_fasta;
sub Extract_mer_seq{

  #$switch_next = 0;
  my $in1 = $_[0];
  my $in2 = $_[1];
  my $in3 = $_[2];
  my $out = $in1 . '_pepuint_fasta';
  $domain_annotation_pepuint_fasta = $out;

  if($in1){
    if(! -e $in1){
      print STDERR "ERROR: file \'$in1\' does not exit. Please check.\n";
      exit(1);
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  if($in2){
    if(! -e $in2){
      print STDERR "ERROR: file \'$in2\' does not exit. Please check.\n";
      exit(1);
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  my %query_mer_actsite_hash;
  my %mer_actsite_count_hash;
  my %mer_id_hash;
  open IN, $in1;
  while(<IN>){
    chomp;
    my $tmp_all = $_;
    if(/^MER.+\s+active\_site\s+(.+)\s+(\d+)\s+active\_site/){
      my $query_id = $2;
      my $mer_actsite = $1;
      $mer_actsite =~ s/\d+//g;
      $mer_actsite =~ s/\s+//g;
      $mer_actsite_count_hash{$query_id} += 1;
      if($mer_actsite_count_hash{$query_id} > 1){
        print "ERROR $query_id\n";
      }
      $query_mer_actsite_hash{$query_id} = $mer_actsite;
    }else{
      print "ERROR $tmp_all\n";
    }
  }
  close IN;

  open IN, $in2;
  open OUT, '>>', $out;
  while(<IN>){
    chomp;
    my $tmp_all = $_;
    if(/^MER(.+)\s+active\_site\s+(.+)\s+(\d+)\s+active\_site/){
      my $query_id = $3;
      my $mer_actsite = $2;
      my $mer_id = $3;
      $mer_actsite =~ s/\d+//g;
      $mer_actsite =~ s/\s+//g;
      if($mer_actsite eq $query_mer_actsite_hash{$query_id}){
        $mer_id_hash{$mer_id} += 1;
      }
    }
  } 
  close IN;

  my %pepunit_hash;
  $/ = '>MER';
  open IN, $in3;
  while(<IN>){
    chomp;
    my $tmp_all = $_;
    if(/^(\d+)\s/){
      my $mer_id = $1;
      $pepunit_hash{$mer_id} = $tmp_all;
    }else{
      print "ERROR $tmp_all\n";
    }
  }

  my @mer_id_array = keys %mer_id_hash;
  foreach my $mer_id (@mer_id_array){
    if($pepunit_hash{$mer_id}){
      print OUT ">MER$pepunit_hash{$mer_id}";
    }else{
      print "ERROR last $mer_id\n";
    }
  }
  close OUT;
}







my $domain_annotation_hits_end;
sub Filter_hits_end{

  $switch_next = 0;
  my $in = $_[0];
  my $out = $in . '_hits_end';
  my $count_match;
  $domain_annotation_hits_end = $out;

  if($in){
    if(! -e $in){
      print STDERR "ERROR: blastp result file \'$in\' does not exit. Please check.\n";
      exit(1);
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  open IN, $in;
  open OUT, '>>', $out;
  $/ = 'Query=';
  while(<IN>){
    chomp;
    if(/\n>\sMER\d+\s/){
      print OUT "Query=$_";
      $switch_next = 1;
    }
  }
  $/ = "\n";
  close IN;
  print OUT "\nQuery= end";
  close OUT;

}

my $domain_annotation_change;
# Change_hits_format($domain_annotation_hits_end);
sub Change_hits_format{

  print STDOUT "Change hits format ...\n";
  $switch_next = 0;
  my $in = $_[0];
  my $out = $in . '_change';
  $domain_annotation_change = $out;

  if($in){
    if(! -e $in){
      print STDERR "ERROR: blastp result file \'$in\' does not exit. Please check.\n";
      exit(1);
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  open IN, $in;
  open OUT, '>>', $out;
  my $count;
  while(<IN>){
    if(/^Query=/){
      $count = 0;
    }
    if(/>\sMER/){
      $count++;
      if($count == 1){
        print OUT "$_";
        next;
      }
      if($count > 1){
        print OUT "Tidy\n\n$_";
        $switch_next++;
        next;
      }
    }
    if(/^Effective\ssearch/){
      print OUT "Tidy\n\n$_";
      $switch_next++;
      next;
    }
    print OUT "$_";
  }
  close IN;
  close OUT;
  print STDOUT "Done.\n";

}

# $in1 = $database_MEROPS_pepunit;
# $in2 = $domain_annotation_hits_end;
my $domain_annotation_rank;
# Rank_all_MEROPS_id($domain_annotation_hits_end);
sub Rank_all_MEROPS_id{

  $switch_next = 0;
  my $in1 = $database_MEROPS_pepunit;
  my $in2 = $_[0];
  my $out = $query_input_tmp . '_pepunit_rank';
  $domain_annotation_rank = $out;

  if($in1){
    if(! -e $in1){
      print STDERR "ERROR: blastp result file \'$in1\' does not exit. Please check.\n";
      exit(1);
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  if($in2){
    if(! -e $in2){
      print STDERR "ERROR: blastp result file \'$in2\' does not exit. Please check.\n";
      exit(1);
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }
  my $rank_pl = $directory . '/rank_unit_hits.pl';
  my $cmd_Rank = "perl $rank_pl $in1 $in2 $out";

  system("$cmd_Rank")==0 or die("failed to execute command \'$cmd_Rank\': $!\n");

  $switch_next++;
}


my $domain_annotation_splited;
#Splited_domain_assembly($domain_annotation_change);
sub Splited_domain_assembly{

  $switch_next = 0;
  my $in = $_[0];
  my $out = $in . '_splited';
  $domain_annotation_splited = $out;

  if($in){
    if(! -e $in){
      print STDERR "ERROR: blastp result file \'$in\' does not exit. Please check.\n";
      exit(1);
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  open IN, $in;
  open OUT, '>>', $out;
  my $mer_id;
  my $score_count;
  my $alines_count;
  my $end;
  my ($query_start, $query_end);
  while(<IN>){
    chomp;

    if(/^Query=\s(.+)/){
      my $query_id  = $1;
      print OUT "Query=$query_id\n";
    }

    if(/^>\sMER(\d+)\s/){
      $mer_id = $1;
      $score_count = 0;
      print OUT "code=$mer_id\n";
    }

    if(/^\sScore\s=/){
      $score_count++;
      $alines_count = 1;
      if($score_count > 1){
        $query_end = $end;
        print OUT "$mer_id $query_start\-$query_end\n";
        $switch_next++;
      }
    }

    if(/^Query\s+(\d+).+?(\d+)/){
      my $start = $1;
      $end = $2;
      if($alines_count == 1){
        $query_start = $start;
      }
      $alines_count++;
    }

    if(/^Tidy/){
      if($score_count == 1){
        $query_end = $end;
        print OUT "$mer_id $query_start\-$query_end\n";
        $switch_next++;
      }
      if($score_count > 1){
        $query_end = $end;
        print OUT "$mer_id $query_start\-$query_end\n";
        $switch_next++;
      }
    }
  }
  close IN;
  close OUT;

}

my $domain_annotation_input;
my $domain_annotation_combine;
#Domain_region_analysis($domain_annotation_splited, $domain_annotation_combine);
my $domain_annotation_result; # = $domain_annotation_input . 'reslut';
#Domain_region_analysis($domain_annotation_input, $domain_annotation_result);
sub Domain_region_analysis{

  $switch_next = 0;
  my $in = $_[0];
  my $out = $_[1];

  if($in){
    if(! -e $in){
      print STDERR "ERROR: blastp result file \'$in\' does not exit. Please check.\n";
      exit(1);
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  open IN, $in;
  open OUT, '>>', $out;

  my ($query_start,$query_end,$mer_id);
  my @array_start;
  my @array_end;
  my %hash_start;
  my %hash_end;

  my $MAX_end_0;
  my $Query_id;
  my $code_id;

  while(<IN>){
    chomp;
    if(/^Query=(.+)/){
 
      $Query_id = $1;

      if(@array_start){

        %hash_start = ();
        %hash_end = ();
        my @array_start_tmp = @array_start;
        my @array_end_tmp = @array_end;
        foreach(@array_start){
          my $start_value = $_;
          my $end_value = shift @array_end_tmp;
          if($hash_start{$start_value}){
            if($hash_start{$start_value} == $end_value or $hash_start{$start_value} > $end_value){

            }else{
              $hash_start{$start_value} = $end_value;
            }
          }
          unless($hash_start{$start_value}){
            $hash_start{$start_value} = $end_value;
          }
          if($hash_end{$end_value}){
            if($hash_end{$end_value} == $start_value or $hash_end{$end_value} < $start_value){

            }else{
              $hash_end{$end_value} = $start_value;
            }          
          }
          unless($hash_end{$end_value}){
            $hash_end{$end_value} = $start_value;
          }
        }

        my $count_code = scalar @array_start;
        if($count_code == 1){
          print OUT "$query_start\-$query_end\n";
          $switch_next++;
        }
        if($count_code > 1){
          my @array_start_sorted = sort{$a <=> $b} @array_start;
          my @array_end_sorted = sort{$a <=> $b} @array_end;
          my @array_start_sorted_tmp = @array_start_sorted;
          my @array_end_sorted_tmp = @array_end_sorted;
          my ($MIN_start,$MAX_end,$start_tmp,$end_tmp);
          my @array_mer_start;
          my @array_mer_end;
          my $start_jump;

          foreach(@array_start_sorted){
            $MIN_start = $_;
            if($start_jump){
              if($start_jump == $_){
                $start_jump = '';
                my $tmp = shift @array_start_sorted_tmp;
                next;
              }
            }
            my $tmp = shift @array_start_sorted_tmp;
            my $middle_start = $array_start_sorted_tmp[0];
            my $MIN_start_0 = $array_start_sorted[0];
            $MAX_end_0 = $array_end_sorted[-1];
 
            $MAX_end = pop @array_end_sorted_tmp;
 
            if($hash_start{$MIN_start_0} == $MAX_end_0){ 
              print OUT "#$MIN_start\-$MAX_end\n";
              $switch_next++;
              last;
            }
 
            if($middle_start){
          
              if($middle_start == $MIN_start){
                next;
              }
 
              if($end_tmp){
  
                if($MIN_start == $end_tmp){
                  $end_tmp = $hash_start{$MIN_start};
                  next;
                }
 
                if(($MIN_start-1) == $end_tmp){
                  $end_tmp = $hash_start{$MIN_start};
                  next;
                }

                if($MIN_start < $end_tmp){
                  if($hash_start{$MIN_start} == $end_tmp){
                    next;
                  }
 
                  if($hash_start{$MIN_start} < $end_tmp){
                    next;
                  }
                  if($hash_start{$MIN_start} > $end_tmp){
                    $end_tmp = $hash_start{$MIN_start};
                    next;
                  }
                }
                if($MIN_start > $end_tmp){
                  push @array_mer_end, $end_tmp;
                  $end_tmp = '';
                }
              }
             
              if($middle_start == $hash_start{$MIN_start}){
                $start_jump = $middle_start;
                push @array_mer_start, $MIN_start;
                $end_tmp = $hash_start{$middle_start};
                next;
              }

              if(($middle_start-1) == $hash_start{$MIN_start}){
                $start_jump = $middle_start;
                push @array_mer_start, $MIN_start;
                $end_tmp = $hash_start{$middle_start};
                next;
              }
 
              if($middle_start < $hash_start{$MIN_start}){
 
                if($hash_start{$MIN_start} == $hash_start{$middle_start}){
                  push @array_mer_start, $MIN_start;
                  $end_tmp = $hash_start{$MIN_start};
                  $start_jump = $middle_start;
                }

                if($hash_start{$MIN_start} < $hash_start{$middle_start}){
                  $end_tmp = $hash_start{$middle_start};
                  push @array_mer_start, $MIN_start;
                  $start_jump = $middle_start;
                }
                if($hash_start{$MIN_start} > $hash_start{$middle_start}){
                  $start_jump = $middle_start;
                  push @array_mer_start, $MIN_start;
                  $end_tmp = $hash_start{$MIN_start};
                  next;
                }
              }
              if($middle_start > $hash_start{$MIN_start}){
                push @array_mer_start, $MIN_start;
                push @array_mer_end, $hash_start{$MIN_start};
              }
            }
            unless($middle_start){
  
              unless($end_tmp){
                unless(@array_mer_end){
                  if(@array_mer_start){
                    push @array_mer_end, $MAX_end_0;
                  }
                  unless(@array_mer_start){
                  push @array_mer_start, $MIN_start;
                  push @array_mer_end, $MAX_end_0;
                  next;
                  }
                }
                if($MIN_start > $array_mer_end[-1]){
                  push @array_mer_start, $MIN_start;
                  push @array_mer_end, $MAX_end_0;
                  next;
                }
              }
 
              if($end_tmp){
                if($hash_start{$MIN_start} == $end_tmp){
                  push @array_mer_end, $MAX_end_0;
                  $end_tmp = '';
                  next;
                }
                if($hash_start{$MIN_start} < $end_tmp){
                  push @array_mer_end, $MAX_end_0;
                  $end_tmp = '';
                  next;
                }
                if($MIN_start < $end_tmp and $hash_start{$MIN_start} > $end_tmp){
                  push @array_mer_end, $MAX_end_0;
                  $end_tmp = '';
                  next;
                }
      
                if($MIN_start == $end_tmp){
                  push @array_mer_end, $MAX_end_0;
                  $end_tmp = '';
                  next;
                }

                if(($MIN_start-1) == $end_tmp){
                  push @array_mer_end, $MAX_end_0;
                  $end_tmp = '';
                  next;
                }
                if($MIN_start > $end_tmp){
                  push @array_mer_end, $end_tmp;
                  $end_tmp = '';
                  push @array_mer_start, $MIN_start;
                  push @array_mer_end, $MAX_end_0;
                  next;
                }
              
              }

              unless($MIN_start){ print "MIN_start $Query_id $code_id\n";}
              unless(@array_mer_end){ print "array_mer_end $Query_id $code_id\nstart @array_start_sorted\nend @array_end_sorted\n";}

            }

          }
          if($end_tmp){
            push @array_mer_end, $MAX_end_0;
          }
          if(@array_mer_start){
            my $number_domain = scalar(@array_mer_start);
            my $count_domain;
            print OUT "#";
            foreach(@array_mer_start){
              $count_domain++;
              my $mer_start = $_;
              my $mer_end = shift @array_mer_end;
              print OUT "$mer_start\-$mer_end";
              print OUT "," if $count_domain < $number_domain;
              print OUT "\n" if $count_domain == $number_domain;
              $switch_next++;
            }
          }
        }
      }

      @array_start = ();
      @array_end = ();
      %hash_start = ();
      %hash_end = ();
      print OUT "Query=$Query_id\n";
    }
    if(/^code=(.+)$/){

      $code_id = $1;

      if(@array_start){

        %hash_start = ();
        %hash_end = ();
        my @array_start_tmp = @array_start;
        my @array_end_tmp = @array_end;
        foreach(@array_start){
          my $start_value = $_;
          my $end_value = shift @array_end_tmp;
          if($hash_start{$start_value}){
            if($hash_start{$start_value} == $end_value or $hash_start{$start_value} > $end_value){
  
            }else{
              $hash_start{$start_value} = $end_value;
            }
          }
          unless($hash_start{$start_value}){
            $hash_start{$start_value} = $end_value;
          }
          if($hash_end{$end_value}){
            if($hash_end{$end_value} == $start_value or $hash_end{$end_value} < $start_value){
 
            }else{
              $hash_end{$end_value} = $start_value;
            }
          }
          unless($hash_end{$end_value}){
            $hash_end{$end_value} = $start_value;
          }
        }

        my $count_code = scalar @array_start;
        if($count_code == 1){
          print OUT "$query_start\-$query_end\n";
          $switch_next++;
        }
        if($count_code > 1){
          my @array_start_sorted = sort{$a <=> $b} @array_start;
          my @array_end_sorted = sort{$a <=> $b} @array_end;
          my @array_start_sorted_tmp = @array_start_sorted;
          my @array_end_sorted_tmp = @array_end_sorted;
          my ($MIN_start,$MAX_end,$start_tmp,$end_tmp);
          my @array_mer_start;
          my @array_mer_end;
          my $start_jump;

          foreach(@array_start_sorted){
            $MIN_start = $_;
            if($start_jump){
              if($start_jump == $_){
                $start_jump = '';
                my $tmp = shift @array_start_sorted_tmp;
                next;
              }
            }
            my $tmp = shift @array_start_sorted_tmp;
            my $middle_start = $array_start_sorted_tmp[0];
            my $MIN_start_0 = $array_start_sorted[0];
            $MAX_end_0 = $array_end_sorted[-1];


            $MAX_end = pop @array_end_sorted_tmp;

            if($hash_start{$MIN_start_0} == $MAX_end_0){
              print OUT "#$MIN_start\-$MAX_end\n";
              $switch_next++;
              last;
            }

            if($middle_start){
 
              if($middle_start == $MIN_start){
                next;
              }

              if($end_tmp){
 
               if($MIN_start == $end_tmp){      
                  $end_tmp = $hash_start{$MIN_start};
                  next;
                }

                if(($MIN_start-1) == $end_tmp){
                  $end_tmp = $hash_start{$MIN_start};
                  next;
                }
 
                if($MIN_start < $end_tmp){
                  if($hash_start{$MIN_start} == $end_tmp){
                    next;
                  }
 
                  if($hash_start{$MIN_start} < $end_tmp){
                    next;
                  }
                  if($hash_start{$MIN_start} > $end_tmp){
                    $end_tmp = $hash_start{$MIN_start};
                    next;
                  }
                }
                if($MIN_start > $end_tmp){
                  push @array_mer_end, $end_tmp;
                  $end_tmp = '';
                }
              }
 
              if($middle_start == $hash_start{$MIN_start}){
                $start_jump = $middle_start;
                push @array_mer_start, $MIN_start;
                $end_tmp = $hash_start{$middle_start};
                next;
              }

              if(($middle_start-1) == $hash_start{$MIN_start}){
                $start_jump = $middle_start;
                push @array_mer_start, $MIN_start;
                $end_tmp = $hash_start{$middle_start};
                next;
              }

              if($middle_start < $hash_start{$MIN_start}){

                if($hash_start{$MIN_start} == $hash_start{$middle_start}){
                  push @array_mer_start, $MIN_start;
                  $end_tmp = $hash_start{$MIN_start};
                  $start_jump = $middle_start;
                }

                if($hash_start{$MIN_start} < $hash_start{$middle_start}){
                  $end_tmp = $hash_start{$middle_start};
                  push @array_mer_start, $MIN_start;
                  $start_jump = $middle_start;
                }
                if($hash_start{$MIN_start} > $hash_start{$middle_start}){
                  $start_jump = $middle_start;
                  push @array_mer_start, $MIN_start;
                  $end_tmp = $hash_start{$MIN_start};
                  next;
                }
              }
              if($middle_start > $hash_start{$MIN_start}){
                push @array_mer_start, $MIN_start;
                push @array_mer_end, $hash_start{$MIN_start};
              }
            }
            unless($middle_start){
 
              unless($end_tmp){
                unless(@array_mer_end){
                  if(@array_mer_start){
                    push @array_mer_end, $MAX_end_0;
                  }
                  unless(@array_mer_start){
                  push @array_mer_start, $MIN_start;
                  push @array_mer_end, $MAX_end_0;
                  next;
                  }
                }
                if($MIN_start > $array_mer_end[-1]){
                  push @array_mer_start, $MIN_start;
                  push @array_mer_end, $MAX_end_0;
                  next;
                }
              }

              if($end_tmp){
                if($hash_start{$MIN_start} == $end_tmp){
                  push @array_mer_end, $MAX_end_0;
                  $end_tmp = '';
                  next;
                }
                if($hash_start{$MIN_start} < $end_tmp){
                  push @array_mer_end, $MAX_end_0;
                  $end_tmp = '';
                  next;
                }
                if($MIN_start < $end_tmp and $hash_start{$MIN_start} > $end_tmp){
                  push @array_mer_end, $MAX_end_0;
                  $end_tmp = '';
                  next;
                }
 
                if($MIN_start == $end_tmp){
                  push @array_mer_end, $MAX_end_0;                
                  $end_tmp = '';
                  next;
                }
  
                if(($MIN_start-1) == $end_tmp){
                  push @array_mer_end, $MAX_end_0;
                  $end_tmp = '';
                  next;
                }
                if($MIN_start > $end_tmp){
                  push @array_mer_end, $end_tmp;
                  $end_tmp = '';
                  push @array_mer_start, $MIN_start;
                  push @array_mer_end, $MAX_end_0;
                  next;
                }
 
              }

              unless($MIN_start){ print "MIN_start $Query_id $code_id\n";}
              unless(@array_mer_end){ print "array_mer_end $Query_id $code_id\nstart @array_start_sorted\nend @array_end_sorted\n";}
  
            }
  
          }
          if($end_tmp){
            push @array_mer_end, $MAX_end_0;
          }
          if(@array_mer_start){
            my $number_domain = scalar(@array_mer_start);
            my $count_domain;
            print OUT "#";
            foreach(@array_mer_start){
              $count_domain++;
              my $mer_start = $_;
              my $mer_end = shift @array_mer_end;
              print OUT "$mer_start\-$mer_end";
              print OUT "," if $count_domain < $number_domain;
              print OUT "\n" if $count_domain == $number_domain;
              $switch_next++;
            }
          }
        }
      }


      @array_start = ();
      @array_end = ();
      %hash_start = ();
      %hash_end = ();
      print OUT "$code_id\t";
    }
    if(/^\d+\s(\d+)\-(\d+)/){
      $query_start = $1;
      $query_end = $2;
     
      if($hash_start{$query_start}){
        if($query_end < $hash_start{$query_start}){
          next;
        }
        if($query_end == $hash_start{$query_start}){
          next;
        }
      }
      if($hash_end{$query_end}){
        if($query_start > $hash_end{$query_end}){
          next;
        }
        if($query_start == $hash_end{$query_end}){
          next;
        }
      }

      my $count_end = 1;

      if($hash_start{$query_start}){
        if($query_end > $hash_start{$query_start}){
          my @array_end_tmp = @array_end;
          @array_end = ();
          foreach(@array_end_tmp){
            my $end = $_;
            if($end == $hash_start{$query_start}){
              push @array_end, $query_end;
              $count_end = 2;
            }else{
              push @array_end, $end;
            }
          }

          my $count_end_end = 1;

          if($count_end == 1){
            my $count_array;
            foreach(@array_start){
              $count_array++;
              if($_ == $query_start){
                $count_array = $count_array - 1;
                if($array_end[$count_array] < $query_end){
                  $array_end[$count_array] = $query_end;
                }
                $count_end_end = 2;
                last;
              }
            }
            if($count_end_end == 1){
              push @array_start, $query_start;
              push @array_end, $query_end;
            }
          }

          $hash_start{$query_start} = $query_end;
        }
        unless($hash_end{$query_end}){
          $hash_end{$query_end} = $query_start;
          next;
        }
      }

      my $count_start = 1;

      if($hash_end{$query_end}){
        if($query_start < $hash_end{$query_end}){
          my @array_start_tmp = @array_start;
          @array_start = ();
          foreach(@array_start_tmp){
            my $start = $_;
            if($start == $hash_end{$query_end}){
              push @array_start, $query_start;
              $count_start = 2;
            }else{
              push @array_start, $start;
            }
          }

          my $count_start_start = 1;

          if($count_start == 1){
            my $count_array;
            foreach(@array_end){
              $count_array++;
              if($_ == $query_end){
                $count_array = $count_array - 1;
                if($array_start[$count_array] > $query_start){
                  $array_start[$count_array] = $query_start;
                }
                $count_start_start = 2;
                last;
              }
            }
            if($count_start_start == 1){
              push @array_start, $query_start;
              push @array_end, $query_end;
            }
          }

          $hash_end{$query_end} = $query_start;
        }

        unless($hash_start{$query_start}){
          $hash_start{$query_start} = $query_end;
          next;
        }
        next;
      }
      push @array_start, $query_start;
      push @array_end, $query_end;
      $hash_end{$query_end} = $query_start;
      $hash_start{$query_start} = $query_end;
    }
    if(/^\d+\s#(.+)/){
      my $domain_tmp = $1;
      if(/,/){
        my @domain_array_tmp = split /,/, $domain_tmp;
        foreach(@domain_array_tmp){
          my $domain_tmp_tmp = $_;
          if($domain_tmp_tmp =~ /(\d+)\-(\d+)/){
            $query_start = $1;
            $query_end = $2;

            if($hash_start{$query_start}){
              next if $query_end < $hash_start{$query_start};
              next if $query_end == $hash_start{$query_start};
            }
            if($hash_end{$query_end}){
              next if $query_start > $hash_end{$query_end};
              next if $query_start == $hash_end{$query_end};
            }
 
            my $count_end = 1;

            if($hash_start{$query_start}){
              if($query_end > $hash_start{$query_start}){
                my @array_end_tmp = @array_end;
                @array_end = ();
                foreach(@array_end_tmp){
                  my $end = $_;
                  if($end == $hash_start{$query_start}){
                    push @array_end, $query_end;
                    $count_end = 2;
                  }else{
                    push @array_end, $end;
                  }
                }
  
                my $count_end_end = 1;

                if($count_end == 1){
                  my $count_array;
                  foreach(@array_start){
                    $count_array++;
                    if($_ == $query_start){
                      $count_array = $count_array - 1;
                      if($array_end[$count_array] < $query_end){
                        $array_end[$count_array] = $query_end;
                      }
                      $count_end_end = 2;
                      last;
                    }
                  }
                  if($count_end_end == 1){
                    push @array_start, $query_start;
                    push @array_end, $query_end;
                  }
                }

                $hash_start{$query_start} = $query_end;
              }
              unless($hash_end{$query_end}){
                $hash_end{$query_end} = $query_start;
                next;
              }
            }

            my $count_start = 1;

            if($hash_end{$query_end}){
              if($query_start < $hash_end{$query_end}){
                my @array_start_tmp = @array_start;
                @array_start = ();
                foreach(@array_start_tmp){
                  my $start = $_;
                  if($start == $hash_end{$query_end}){
                    push @array_start, $query_start;
                    $count_start = 2;
                  }else{
                    push @array_start, $start;
                  }
                }
  
                my $count_start_start = 1;

                if($count_start == 1){
                  my $count_array;
                  foreach(@array_end){
                    $count_array++;
                    if($_ == $query_end){
                      $count_array = $count_array - 1;
                      if($array_start[$count_array] > $query_start){
                        $array_start[$count_array] = $query_start;
                      }
                      $count_start_start = 2;
                      last;
                    }
                  }
                  if($count_start_start == 1){
                    push @array_start, $query_start;
                    push @array_end, $query_end;
                  }
                }

                $hash_end{$query_end} = $query_start;
              }
              unless($hash_start{$query_start}){
                $hash_start{$query_start} = $query_end;
                next;
              }
              next;
            }

            push @array_start, $query_start;
            push @array_end, $query_end;
            $hash_end{$query_end} = $query_start;
            $hash_start{$query_start} = $query_end;

          }
        }
      }else{
        if($domain_tmp =~ /(\d+)\-(\d+)/){
          $query_start = $1;
          $query_end = $2;

          if($hash_start{$query_start}){
            next if $query_end < $hash_start{$query_start};
            next if $query_end == $hash_start{$query_start};
          }
          if($hash_end{$query_end}){
            next if $query_start > $hash_end{$query_end};
            next if $query_start == $hash_end{$query_end};
          }
    

          my $count_end = 1;

          if($hash_start{$query_start}){
            if($query_end > $hash_start{$query_start}){
              my @array_end_tmp = @array_end;
              @array_end = ();
              foreach(@array_end_tmp){
                my $end = $_;
                if($end == $hash_start{$query_start}){
                  push @array_end, $query_end;
                  $count_end = 2;
                }else{
                  push @array_end, $end;
                }
              }

              my $count_end_end = 1;

              if($count_end == 1){
                my $count_array;
                foreach(@array_start){
                  $count_array++;
                  if($_ == $query_start){
                    $count_array = $count_array - 1;
                    if($array_end[$count_array] < $query_end){
                      $array_end[$count_array] = $query_end;
                    }
                    $count_end_end = 2;
                    last;
                  }
                }
                if($count_end_end == 1){
                  push @array_start, $query_start;
                  push @array_end, $query_end;
                }
              }
  
              $hash_start{$query_start} = $query_end;
            }
            unless($hash_end{$query_end}){
              $hash_end{$query_end} = $query_start;
              next;
            }
          }

          my $count_start = 1;
 
          if($hash_end{$query_end}){
            if($query_start < $hash_end{$query_end}){
              my @array_start_tmp = @array_start;
              @array_start = ();
              foreach(@array_start_tmp){
                my $start = $_;
                if($start == $hash_end{$query_end}){
                  push @array_start, $query_start;
                  $count_start = 2;
                }else{
                  push @array_start, $start;
                }
              }

              my $count_start_start = 1;

              if($count_start == 1){
                my $count_array;
                foreach(@array_end){
                  $count_array++;
                  if($_ == $query_end){
                    $count_array = $count_array - 1;
                    if($array_start[$count_array] > $query_start){
                      $array_start[$count_array] = $query_start;
                    }
                    $count_start_start = 2;
                    last;
                  }
                }
                if($count_start_start == 1){
                  push @array_start, $query_start;
                  push @array_end, $query_end;
                }
              }

              $hash_end{$query_end} = $query_start;
            }

            unless($hash_start{$query_start}){
              $hash_start{$query_start} = $query_end;
              next;
            }
            next;
          }
 
          push @array_start, $query_start;
          push @array_end, $query_end;
          $hash_end{$query_end} = $query_start;
          $hash_start{$query_start} = $query_end;
        }
      }
    }
  }
  close IN;
  close OUT;  

}

#$domain_annotation_input;
#Create_input_for_domain_analysis($domain_annotation_combine, $domain_annotation_rank);
sub Create_input_for_domain_analysis{
  
  $switch_next = 0;
  my $in1 = $_[0];
  my $in2 = $_[1];
  my $out = $in1 . '_input';
  $domain_annotation_miss_match_input = $out;

  if($in1){
    if(! -e $in1){
      print STDERR "ERROR: file \'$in1\' does not exit. Please check.\n";
      exit(1);
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }

  if($in2){
    if(! -e $in2){
      print STDERR "ERROR: file \'$in2\' does not exit. Please check.\n";
      exit(1);
    }
  }else{
    print STDERR "FATAL ERROR: The primary code of the main program is modified. Please redownload the main program 'PepTidy.pl' and replace the modified 'PepTidy.pl'.\n";
    exit(1);
  }
  
  my $input_pl = $directory . '/create_input_for_domain_region_analysis.pl';
  my $cmd_Input = "perl $input_pl $in1 $in2 $out";
 
  system("$cmd_Input")==0 or die("failed to execute command \'$cmd_Input\': $!\n");
 
  $switch_next++;
  
}


           ###################### Actsites Annotation Peplines ##########################

#if at protease lever
#my $key_protease = 1;

#if at pepunit lever
#$key_protease = 0;

#Step_1: check the format of Blastp result file
#if Done and next $switch_next = 2;
#Check_blastp_input();

#Step_1.1: if blastp at protease or database lever, filter sequences in database from blastp results
#if($switch_next ==2);
#Filter_MEROPS_protease(,'database'/'protease');

#Step_2: match Actsites
#Match_actsite();

#Step_3: filter Actsite matching results to find out all correct results
#Filter_actsite_matching_result();

#Step_4: rearrrange query actsite by numerical value from small to large
#Rearrange_actsite_matching_result(); 

#step_5: delete repeated actsite type in one query_id
#Delete_repeat_in_actsite_matching_result();

#step_6: separate one query_id with single actsite type from one query_id with poly actsite type
#Separate_actsite_matching_result();

#step_7: separate one query_id with single amino acid type from one query_id with poly actsite amino acid type
#Separate_AA_actsite_matching_result();

#step_8: check the inclusion relationship of actsite types within one query_id and delete the inclued actsite type
#Delete_inclusion_actsite_matching_result();

#step_9: combine all actsites annotation results
#combine_actsites();









            ###################### Domain Annotation Peplines ##########################

my $domain_annotation_match_input;
#create_input_for_domain_match_analysis($domain_annotation_combine, $actsites_annotation_all);
sub create_input_for_domain_match_analysis{
  
  my $in1 = $_[0];
  my $in2 = $_[1];
  my $out = $domain_annotation_combine . '_match_input';
  $domain_annotation_match_input = $out; 
  my %query_id_match_hash;
  open IN, $in1;
  while(<IN>){
    chomp;
    if(/^Query=(.+)/){
      my $query_id = $1;
      $query_id_match_hash{$query_id} += 1;
    }
  }
  close IN;

  my %query_id_match_repeat_hash;
  open IN, $in2;
  while(<IN>){
    chomp;
    if(/^MER(\d+)\s.+(QER\d+)\s/){
      my $mer_id = $1;
      my $query_id = $2;
      $query_id_match_repeat_hash{$query_id} += 1;
      if($query_id_match_repeat_hash{$query_id} == $ii){
        if($query_id_match_hash{$query_id}){
          my $code_domain = $MEROPS_code_hash{$mer_id};
          if($code_domain =~ /(.+)\./){
            $code_domain = $1;
          }
          $query_id_match_hash{$query_id} = $code_domain;
        }
        #unless($query_id_match_hash{$query_id}){

        #}
      }
    }
  }
  close IN;
  
  open IN, $in1;
  open OUT, '>>', $out;
  while(<IN>){
    my $tmp_all = $_;
    print OUT "$tmp_all";
    chomp;
    if(/^Query=(.+)/){
      my $query_id = $1;
      if($query_id eq 'end'){
        next;
      }
      print OUT "code=$query_id_match_hash{$query_id}\n";
    }
  }
  close IN;
  close OUT;

}



#Fetch_domain_blastp_result();
sub Fetch_domain_blastp_result{

  print STDOUT "Fetch blastp result for domain analysis ...\n";
  $switch_option = 0;
  #$blastp_pepunit_out = $query_input_tmp . '_pepunit_blastp';
  $actsites_annotation_all = $query_input_tmp . '_actsite_annotation_all';
  $actsites_rearrange_all = $query_input_tmp . '_actsite_rearrange_all';
  if($actsite_protease_key){
    if(! -e $actsites_annotation_protease){
      print STDERR "ERROR: no \'$actsites_annotation_protease\' file.\n";
      exit(1);
    }
    if(! -e $actsite_out_protease_rearrange){
      print STDERR "ERROR: no \'$actsite_out_protease_rearrange\' file.\n";
      exit(1);
    }
    if($actsite_pepunit_key){
      if(! -e $actsites_annotation_pepunit){
        print STDERR "ERROR: no \'$actsites_annotation_pepunit\' file.\n";
        exit(1);
      }
      if(! -e $actsite_out_pepunit_rearrange){
        print STDERR "ERROR: no \'$actsite_out_pepunit_rearrange\' file.\n";
        exit(1);
      }
      my $cmd_combine_all = "cat $actsites_annotation_protease $actsites_annotation_pepunit > $actsites_annotation_all";
      system("$cmd_combine_all")==0 or die("failed to execute command \'$cmd_combine_all\': $!\n");
      $cmd_combine_all = "cat $actsite_out_protease_rearrange $actsite_out_pepunit_rearrange > $actsites_rearrange_all";
      system("$cmd_combine_all")==0 or die("failed to execute command \'$cmd_combine_all\': $!\n");
    }
    unless($actsite_pepunit_key){
      $actsites_annotation_all = $actsites_annotation_protease;
      $actsites_rearrange_all = $actsite_out_protease_rearrange;
    }
  }
  unless($actsite_protease_key){
    if($actsite_pepunit_key){
      if(! -e $actsites_annotation_pepunit){
        print STDERR "ERROR: no \'$actsites_annotation_pepunit\' file.\n";
        exit(1);
      }
      if(! -e $actsite_out_pepunit_rearrange){
        print STDERR "ERROR: no \'$actsite_out_pepunit_rearrange\' file.\n";
        exit(1);
      }
      $actsites_annotation_all = $actsites_annotation_pepunit;
      $actsites_rearrange_all = $actsite_out_pepunit_rearrange;
    }
  }

  my $in_miss_code = $MEROPS_miss_code;
  my $in_pepunit = $database_MEROPS_pepunit;
  if(! -e $in_pepunit){
    print STDERR "ERROR: no \'$in_pepunit\' file.\n";
    exit(1);
  }
  if(! -e $in_miss_code){
    print STDERR "ERROR: no \'$in_pepunit\' file.\n";
    exit(1);
  }
  if(! -e $actsites_annotation_all){
    print STDERR "ERROR: no \'$actsites_annotation_all\' file.\n";
    exit(1);
  }
  if(! -e $actsites_rearrange_all){
    print STDERR "ERROR: no \'$actsites_rearrange_all\' file.\n";
    exit(1);
  }
  open IN, $in_pepunit;
  while(<IN>){
    chomp;
    if(/^>MER(\d+)\s.+\[(.+)\]\#/){
      my $mer_id = $1;
      $MEROPS_code_hash{$mer_id} = $2;
    }
  }
  close IN;

  open IN, $in_miss_code;
  while(<IN>){
    chomp;
    if(/^(\d+)\s+(.+)$/){
      my $mer_id = $1;
      $MEROPS_code_hash{$mer_id} = $2;
    }
  }
  close IN;

  unless(%MEROPS_code_hash){
    print STDERR "ERROR: Can't read MEROPS codes. File \'$in_pepunit\' dose not have a proper format.\n";
    exit(1);
  }

  my %query_id_repeat_hash;
  open IN, $actsites_annotation_all;
  while(<IN>){
    chomp;
    if(/(QER\d+)/){
      my $query_id = $1;
      $query_id_repeat_hash{$query_id} += 1;
      if($query_id_repeat_hash{$query_id} > $count_max_actsite_repeat){
        $count_max_actsite_repeat = $query_id_repeat_hash{$query_id};
      }
    }
  }
  $count_max_actsite_repeat++;
  close IN;
  unless(%query_id_repeat_hash){
    print STDERR "ERROR: Can't read MEROPS codes. File \'$actsites_annotation_all\' dose not have a proper format.\n";
    exit(1);
  }

  my $i;
  my %query_id_miss_hash;
  #$balstp_result_for_domain_analysis = $workDir . '/' . 'pro.fa_domain_blastp'; ########change#
  $balstp_result_for_domain_analysis = $query_input_tmp . '_domain_blastp';
  for($i = 1; $i < $count_max_actsite_repeat; $i++){

    my $out = "$balstp_result_for_domain_analysis" . "$i";
    %query_id_repeat_hash = ();
    my %sbjct_actsite_hash;
    open IN, $actsites_annotation_all;
    while(<IN>){
      chomp;
      if(/^MER.+\s+active\_site\s+(.+)\s+(QER\d+)\s/){
        my $tmp_sbjct = $1;
        my $query_id = $2;
        $query_id_repeat_hash{$query_id} += 1;
        if($query_id_repeat_hash{$query_id} == $i){
          $tmp_sbjct =~ s/\d+//g;
          $tmp_sbjct =~ s/\s//g;
          $sbjct_actsite_hash{$query_id} = $tmp_sbjct;
        }
        else{
          next;
        }
      }
    }
    close IN;
    unless(%sbjct_actsite_hash){
      print STDERR "ERROR: Can't read MEROPS codes. File \'$actsites_annotation_all\' dose not have a proper format.\n";
      exit(1);
    }

    %query_id_repeat_hash = ();
    my %mer_id_actsite_hash;
    open IN, $actsites_rearrange_all;
    while(<IN>){
      chomp;
      if(/^MER(.+)\s+active\_site\s+(.+)\s+(QER\d+)\s/){
        my $mer_id = $1;
        my $tmp_sbjct = $2;
        my $query_id = $3;
        unless($sbjct_actsite_hash{$query_id}){
          next;
        }
        $tmp_sbjct =~ s/\d+//g;
        $tmp_sbjct =~ s/\s//g;
        if($sbjct_actsite_hash{$query_id} eq $tmp_sbjct){
          $query_id_repeat_hash{$query_id} += 1;
          if($query_id_repeat_hash{$query_id} == 1){
            $mer_id_actsite_hash{$query_id} = $mer_id;
          }
          if($query_id_repeat_hash{$query_id} > 1){
            my $mer_id_tmp = '#' . $mer_id;
            $mer_id_actsite_hash{$query_id} .= $mer_id_tmp;
          }
        }
      }
    }
    close IN;
    unless(%mer_id_actsite_hash){
      print STDERR "ERROR: Can't read MEROPS codes. File \'$actsites_rearrange_all\' dose not have a proper format.\n";
      exit(1);
    }

    open IN, $blastp_pepunit_out_first;
    open OUT, '>>', $out;
    my $count_all;
    my $count_k;
    my $print_balstp_key = 1;
    my $query_match_key = 1;
    my %mer_id_match_hash;
    my @mer_id_match_array;
    my $empty_key = 3;
    my $query_id_old;
    while(<IN>){
      my $tmp_all = $_;
      chomp;
      if(/^Query=\s(.+)/){
        my $query_id = $1;
        %mer_id_match_hash = ();
        if($mer_id_actsite_hash{$query_id}){
          $count_all++;
          if($empty_key == 2){
            #print "$query_id_old $mer_id_match_array[0] $MEROPS_code_hash{$mer_id_match_array[0]}\n";
            $query_id_miss_hash{$query_id_old} += 1;
          }
          if($empty_key == 1){
            $empty_key = 2;
            $count_k++;
          }
          $query_match_key = 2;
          $print_balstp_key = 2;
          @mer_id_match_array = ();
          @mer_id_match_array = split /#/, $mer_id_actsite_hash{$query_id};
          foreach my $mer_id_tmp (@mer_id_match_array){
            $mer_id_match_hash{$mer_id_tmp} += 1;
          }
          $empty_key = 2;
          $query_id_old = $query_id;
        }
        unless($mer_id_actsite_hash{$query_id}){
          $query_match_key = 1;
          $print_balstp_key = 1;
        }
      }
      if($query_match_key == 2){
        if(/>\sMER(\d+)\s/){
          my $mer_id = $1;
          if($mer_id_match_hash{$mer_id}){
            $print_balstp_key = 2;
            $empty_key = 1;
            $switch_option = 1;
          }
          unless($mer_id_match_hash{$mer_id}){
            $print_balstp_key = 1;
          }
        }
        if(/^Lambda\s+K\s+H/){
          $print_balstp_key = 2;
        }
      }

      if($print_balstp_key == 2){
        print OUT "$tmp_all";
      }
    }
    if($empty_key == 2){
      #print "$query_id_old $mer_id_match_array[0] $MEROPS_code_hash{$mer_id_match_array[0]}\n";
      $query_id_miss_hash{$query_id_old} += 1;
    }
  }
  close IN;
  close OUT;

  if(%query_id_miss_hash){
    open IN, $blastp_pepunit_out_first;
    $blastp_result_for_miss_domain_analysis = $query_input_tmp . '_miss_domain_blastp';
    open OUT, '>>', $blastp_result_for_miss_domain_analysis;
    my $print_balstp_key = 1;
    while(<IN>){
      my $tmp_all = $_;
      chomp;
      if(/^Query=\s(.+)/){
        my $query_id = $1;
        if($query_id_miss_hash{$query_id}){
          $print_balstp_key = 2;
          unless($switch_option){
            $switch_option = 2;
          }
          if($switch_option == 1){
            $switch_option = 3;
          }
        }
        unless($query_id_miss_hash{$query_id}){
          $print_balstp_key = 1;
        }
      }
      if($print_balstp_key == 2){
        print OUT "$tmp_all";
      }
    }
    print OUT "\nQuery= end";
    close OUT;
  }
  print STDOUT "Done.\n";

}

#$domain_annotation_miss_match_input;
#create_input_for_miss_domain_analysis();
sub create_input_for_miss_domain_analysis{

  my %code_repeat_hash;
  open IN, $domain_annotation_miss_match_result_tmp;
  my $out = $domain_annotation_miss_match_result_tmp . '_domain_input';
  $domain_annotation_miss_match_input = $out;
  open OUT, '>>', $out;
  while(<IN>){
    my $tmp_all = $_;
    chomp;
    if(/^Query=/){
      print OUT "$tmp_all";
      %code_repeat_hash =();
    }
    if(/^(.+)\..+\s(.+)$/){
      my $code_id = $1;
      my $range = $2;
      $code_repeat_hash{$code_id} += 1;
      if($code_repeat_hash{$code_id} == 1){
        print OUT "code=$code_id\n";
      }
      print OUT "1\t$range\n";
    }
  }
  close IN;
  close OUT;

}


my %miss_annotation_tmp_hash;
sub create_miss_domain_annotation{


  open IN, $domain_annotation_miss_match_result;
  my $query_id_tmp;
  while(<IN>){
    chomp;
    my $tmp_all = $_;
    if(/^Query=(.+)/){
      $query_id_tmp = $1;
      next;
    }
    if(/^(.+)\t(.+)/){
      my $family = $1;
      my $domain = $2;
      my $query_id_family = $query_id_tmp . $family;
      $miss_annotation_tmp_hash{$query_id_family} = $domain;
    }
  }
  close IN;
  close OUT;
}

my %domain_annotation_hash;
my %family_annotation_hash;
sub create_domain_annotation{

  my %domain_annotation_tmp_hash;
  if($switch_next){
    #my %domain_annotation_tmp_hash;
    open IN, $domain_annotation_match_result;
    my $query_id_tmp;
    while(<IN>){
      chomp;
      my $tmp_all = $_;
      if(/^Query=(.+)/){
        $query_id_tmp = $1;
        next;
      }
      $domain_annotation_tmp_hash{$query_id_tmp} = $tmp_all;
    }
    close IN;
  }

  if($switch_option == 1 || $switch_option == 3){
    my %query_id_repeat_hash;
    open IN, $actsites_annotation_all;
    while(<IN>){
      chomp;
      if(/^MER(\d+)\s+.+(QER\d+)\s/){
        my $query_id = $2;
        my $mer_id = $1;
        my $mer_query_id = 'MER' . $mer_id . $query_id;
        $query_id_repeat_hash{$query_id} += 1;
        if($query_id_repeat_hash{$query_id} == $ii){
          if($domain_annotation_tmp_hash{$query_id}){
            if($domain_annotation_tmp_hash{$query_id} =~ /^(.+)\t(.+)/){
              my $family = $1;
              my $domain = $2;
              $domain_annotation_hash{$mer_query_id} = $domain;
              $family_annotation_hash{$mer_query_id} = $family;
            }
          }
          unless($domain_annotation_tmp_hash{$query_id}){
            if($MEROPS_code_hash{$mer_id} =~ /(.+)\./){
              my $family = $1;
              my $query_id_family = $query_id . $family;
              $domain_annotation_hash{$mer_query_id} = $miss_annotation_tmp_hash{$query_id_family};
              $family_annotation_hash{$mer_query_id} = $family;
            }else{
              print STDERR "ERROR: $MEROPS_code_hash{$mer_id} $query_id MER$mer_id\n"; ###test
              exit(1);
            }
          }
        }
      }
    }
    close IN;
  }

  if($switch_option == 2){
    open IN, $actsites_annotation_all;
    while(<IN>){
      chomp;
      if(/^MER(\d+)\s+.+(QER\d+)\s/){
        my $query_id = $2;
        my $mer_id = $1;
        my $mer_query_id = 'MER' . $mer_id . $query_id;
        if($MEROPS_code_hash{$mer_id} =~ /(.+)\./){
          my $family = $1;
          my $query_id_family = $query_id . $family;
          $domain_annotation_hash{$mer_query_id} = $miss_annotation_tmp_hash{$query_id_family};
          $family_annotation_hash{$mer_query_id} = $family;
        }else{
          print STDERR "ERROR: $query_id MER$mer_id\n";
          exit(1);
        }
      }
    }
    close IN;
  }

}


#create_merops_annotation($protease_same_query);
sub create_merops_annotation{

  $query_annotation_results = $workDir  . '/' . $query_file_name . '_annotation_results_query';
  $uniport_annotation_results = $workDir  . '/' . $query_file_name . '_annotation_results_uniport';
  my $in = $_[0];

  my %MEROPS_code_miss_hash;
  my %MEROPS_protease_code_hash;
  my %merops_domain_hash;
  unless($custom || $annotation_file){
    open IN, $MEROPS_miss_code;
    while(<IN>){
      chomp;
      if(/^(\d+)\s+(.+)$/){
        my $mer_id = $1;
        $MEROPS_code_miss_hash{$mer_id} = $2;
      }
    }
    close IN;

    open IN, $database_MEROPS_protease;
    while(<IN>){
      chomp;
      if(/^>MER(\d+)\s/){
        my $mer_id = $1;
        if(/^>MER\d+\s.*unit:\s(.+)\}/){
          my $domain = $1;
          $merops_domain_hash{$mer_id} = $domain;
        }else{
          $merops_domain_hash{$mer_id} = 'NA';
        }

        if(/^>MER\d+\s.+\[(.+)\]\#/){
          $MEROPS_protease_code_hash{$mer_id} = $1;
        }
      }
    }
    close IN;
  }

  my %uniport_id_hash;
  open IN, $mer_uniport_id_sum;
  while(<IN>){
    chomp;
    if(/^(\d+)\#(.+)$/){
      my $mer_id = $1;
      my $uniport_id = $2;
      $uniport_id_hash{$mer_id} = $uniport_id;
    }
  }
  close IN;

  open IN, $in;
  while(<IN>){
    chomp;
    if(/^QER\d+\#MER(\d+)$/){
      my $mer_id = $1;
      my $uniport_id = $uniport_id_hash{$mer_id};
      $query_ori_uniport_id_hash{$uniport_id} += 1;
    }
  }
  close IN;

  open OUT, '>>', $uniport_annotation_results;
  open IN, $uniport_annotation_description;
  my $tmp_uniport;
  my $print_key = 1;
  my $print_key_key = 1;
  my $rec_name_key = 1;
  my $uniport_id_DE;
  while(<IN>){
    my $tmp_all = $_;
    $print_key_key = 1;
    chomp;
    if(/^(ID\s+.+)/){
      $tmp_uniport = $1;
      $print_key = 1;
      $print_key_key = 1;
    }
    if(/^(AC\s+(.+))/){
      $tmp_uniport .= "\n$1";
      my $AC_tmp = $1;
      my @uniport_id_array = ($AC_tmp =~ /[A-Za-z0-9]+/g);
      foreach my $uniport_id (@uniport_id_array){
        if($query_ori_uniport_id_hash{$uniport_id}){
          $annotation_uniport_count++;
          print OUT "$tmp_uniport\n";
          $switch_ori = 2;
          $print_key = 2;
          $uniport_id_DE = $uniport_id;
          last;
        }
      }
      next;
    }
    if($print_key == 2){
      if(/^DE\s+RecName\:/){
        $print_key_key = 2;
        $rec_name_key = 2;
      }
      if(/^OS\s+/){
        $print_key_key = 2;
        $rec_name_key = 1;
      }
      if(/^OX\s+/){
        $print_key_key = 2;
        $rec_name_key = 1;
      }
      if(/^RN\s+/){
        $print_key_key = 2;
        $rec_name_key = 1;
      }
      if(/^RX\s+/){
        $print_key_key = 2;
      }
      if(/^RT\s+/){
        $print_key_key = 2;
      }
      if(/^RL\s+/){
        $print_key_key = 2;
      }
      if(/^CC\s+/){
        $print_key_key = 2;
      }
      if(/^DR\s+/){
        $print_key_key = 2;
      }
      if(/^PE\s+/){
        $print_key_key = 2;
      }
      if(/^FT\s+/){
        $print_key_key = 2;
      }
      if(/^SQ\s+/){
        $print_key_key = 2;
      }
      if(/^\/\//){
        $print_key_key = 2;
      }
    }
    if($print_key_key == 2){
      print OUT "$tmp_all";
    }
    if($rec_name_key == 2){
      $uniport_id_DE_hash{$uniport_id_DE} .= $tmp_all;
    }
  }
  close IN;
  close OUT;

  my %actsite_hash;
  open IN, $MEROPS_protease_actsites_all;
  while(<IN>){
    chomp;
    if(/^(\d+)\#\"(.+)\"\,/){
      my $mer_id = $1;
      my $tmp = $2;
      $actsite_hash{$mer_id} = $tmp;
    }
  }
  close IN;

  my %query_id_hash;
  my %mer_id_doman_hash;
  open IN, $in;
  open OUT, '>>', $query_annotation_results;
  my $count_end;
  while(<IN>){
    chomp;
    if(/^(QER(\d+))\#MER(\d+)$/){
      $count_end++;
      my $query_id = $1;
      my $query_id_num = $2;
      my $mer_id = $3;
      $query_id_annotation_seq_hash{$query_id_num} += 1;
      my $uniport_id = $uniport_id_hash{$mer_id};
      $query_ori_id_hash{$query_id} += 1;
      $query_id_hash{$query_id} += 1;
      if($query_id_hash{$query_id} == 1){
        $annotation_ori_count++;
        %mer_id_doman_hash = ();
        if($count_end > 1){
          print OUT "end\n\n\n";
        }
      }
      if($query_id_hash{$query_id} > 1){
        next;
      }
      my $query_name = $query_all_name_hash{$query_id_num};
      print OUT "QN   $query_name\n";
      if($switch_ori == 2){
        $switch_ori = 3;
      }
      if($switch_ori == 1){
        $switch_ori = 4;
      }
      my $actsites = $actsite_hash{$mer_id};
      my $AMB_actsites = $actsite_AMB_mer_id_hash{$mer_id};
      my %AMB_actsite_hash;
      my @AMB_actsite_cha_array = ($AMB_actsites =~ /[AMBJOXZ]/g);
      my @AMB_actsite_num_array = ($AMB_actsites =~ /\d+/g);
      my $count_cha = scalar @AMB_actsite_cha_array;
      my $iii;
      for($iii = 0; $iii < $count_cha; $iii++){
         $AMB_actsite_hash{$AMB_actsite_num_array[$iii]} = $AMB_actsite_cha_array[$iii];
      }
      my @actsites_array = ($actsites =~ /\w+/g);
      my @active_array;
      my @metal_array;
      my @binding_array;
      foreach my $actsites_tmp (@actsites_array){
        my $actsites_num = $actsites_tmp;
        $actsites_num =~ s/[A-Za-z]+//g;
        $actsites_num =~ s/\s+//g;
        if($AMB_actsite_hash{$actsites_num} eq 'A'){
          push @active_array, $actsites_tmp;
        }
        elsif($AMB_actsite_hash{$actsites_num} eq 'M'){
          push @metal_array, $actsites_tmp;
        }
        elsif($AMB_actsite_hash{$actsites_num} eq 'B'){
          push @binding_array, $actsites_tmp;
        }
        elsif($AMB_actsite_hash{$actsites_num} eq 'J'){
          push @active_array, $actsites_tmp;
          push @binding_array, $actsites_tmp;
        }
        elsif($AMB_actsite_hash{$actsites_num} eq 'O'){
          push @active_array, $actsites_tmp;
          push @metal_array, $actsites_tmp;
        }
        elsif($AMB_actsite_hash{$actsites_num} eq 'Z'){
          push @binding_array, $actsites_tmp;
          push @metal_array, $actsites_tmp;
        }
        elsif($AMB_actsite_hash{$actsites_num} eq 'X'){
          push @active_array, $actsites_tmp;
          push @metal_array, $actsites_tmp;
          push @binding_array, $actsites_tmp;
        }else{
          print STDERR "ERRPR: actsite changing is error. Please redowload software Peptidy.\n";
          exit(1);
        }
      }
      if(@active_array){
        print OUT "AT   ACT_SITE: @active_array\n";
      }
      if(@metal_array){
        print OUT "AT   METAL: @metal_array\n";
      }
      if(@binding_array){
        print OUT "AT   BINDING: @binding_array\n";
      }
      unless($custom){
        unless($annotation_file){
          my $domain_query = $merops_domain_hash{$mer_id};
          $mer_id_doman_hash{$domain_query} += 1;
          my $family = $MEROPS_protease_code_hash{$mer_id};
          unless($MEROPS_code_miss_hash{$mer_id}){
            print OUT "DO   Domain: $domain_query\n";
            print OUT "FM   Family: $family\n";

          }
          if($MEROPS_code_miss_hash{$mer_id}){
            print OUT "DO   Domain: $domain_query \#$family\n";
            print OUT "FM   Family: $MEROPS_code_miss_hash{$mer_id}\n";
          }
        }
        print OUT "MD   MEROPS_id: $mer_id\n";
      }
      if($custom){
        print OUT "MD   MEROPS_id: $custom_database_seq_name_hash{$mer_id}\n";
      }
      print OUT "EV   E-Value: same\n";
      print OUT "SC   Score: MAX\n";
      print OUT "UP   Uniport_id: $uniport_id\n";
      print OUT "$uniport_id_DE_hash{$uniport_id}";
    }

  }
  print OUT "end\n\n\n";

}


sub create_query_annotation{

  $query_annotation_results = $workDir  . '/' . $query_file_name . '_annotation_results_query';
  $uniport_annotation_results = $workDir  . '/' . $query_file_name . '_annotation_results_uniport';

  my %query_id_seq_hash;
  my @query_id_seq_array;
  my %query_acstites_hash;
  my %mer_id_hash;
  my %uniport_id_hash;
  my %query_uniport_id_hash;
  my %mer_query_uniport_id_hash;
  open IN, $mer_uniport_id_sum;
  while(<IN>){
    chomp;
    if(/^(\d+)\#(.+)$/){
      my $mer_id = $1;
      my $uniport_id = $2;
      $uniport_id_hash{$mer_id} = $uniport_id;
    }
  }
  close IN;

  open IN, $actsites_annotation_all;
  my $count_uniport_seq = 1;
  while(<IN>){
    chomp;
    if(/^MER(\d+)\t.+(QER(\d+))\tactive_site\t(.+)\s$/){
      my $mer_id = $1;
      my $query_id = $2;
      my $mer_query_id = 'MER' . $mer_id . $query_id;
      my $query_id_num = $3;
      my $query_acsite = $4;
      $query_id_annotation_seq_hash{$query_id_num} += 1;
      if($merops_switch){
        if($query_ori_id_hash{$query_id}){
          next;
        }
      }
      my $id_tmp = $mer_query_id . '#';
      my $uniport_id = $uniport_id_hash{$mer_id};
      $query_id_seq_hash{$query_id_num} .= $id_tmp;
      $query_acstites_hash{$mer_query_id} = $query_acsite;
      $mer_id_hash{$mer_query_id} = $mer_id;
      unless($query_ori_uniport_id_hash{$uniport_id}){
        $query_uniport_id_hash{$uniport_id} += 1;
      }
      $mer_query_uniport_id_hash{$mer_query_id} = $uniport_id;
    }
  }
  close IN;


  my %evalue_hash;
  my %score_hash;
  if($actsite_protease_key){
    open IN, $blastp_protease_out;
    my $evalue_key = 1;
    my $mer_query_id_tmp;
    my $query_id_tmp;
    while(<IN>){
      chomp;
      if(/^Query=\s+(.+)/){
        $query_id_tmp = $1;
        $evalue_key = 1;
      }
      if(/^>\s+MER(\d+)\s/){
        my $mer_id = $1;
        $mer_query_id_tmp = 'MER' . $mer_id . $query_id_tmp;
        if($mer_id_hash{$mer_query_id_tmp}){
          $evalue_key = 2;
        }
        unless($mer_id_hash{$mer_query_id_tmp}){
          $evalue_key = 1;
        }
      }
      if($evalue_key == 2){
        if(/^\s+Score\s+=\s+(.+)\,\s+Expect\s+=\s+(.+)\,\s+Method/){
          my $score = $1;
          my $evalue = $2;
          $evalue_hash{$mer_query_id_tmp} = $evalue;
          $score_hash{$mer_query_id_tmp} = $score;
        }
      }
    }
    close IN;
  }

  if($actsite_pepunit_key){
    open IN, $blastp_pepunit_out;
    my $evalue_key = 1;
    my $mer_query_id_tmp;
    my $query_id_tmp;
    while(<IN>){
      chomp;
      if(/^Query=\s+(.+)/){
        $query_id_tmp = $1;
        $evalue_key = 1;
      }
      if(/^>\s+MER(\d+)\s/){
        my $mer_id = $1;
        $mer_query_id_tmp = 'MER' . $mer_id . $query_id_tmp;
        if($mer_id_hash{$mer_query_id_tmp}){
          $evalue_key = 2;
        }
        unless($mer_id_hash{$mer_query_id_tmp}){
          $evalue_key = 1;
        }
      }
      if($evalue_key == 2){
        if(/^\s+Score\s+=\s+(.+)\,\s+Expect\s+=\s+(.+)\,\s+Method/){
          my $score = $1;
          my $evalue = $2;
          $evalue_hash{$mer_query_id_tmp} = $evalue;
          $score_hash{$mer_query_id_tmp} = $score;
        }
      }
    }
    close IN;
  }

  open OUT, '>>', $uniport_annotation_results;
  open IN, $uniport_annotation_description;
  my $tmp_uniport;
  my $print_key = 1;
  my $print_key_key = 1;
  my $rec_name_key = 1;
  my $uniport_id_DE;
  while(<IN>){
    my $tmp_all = $_;
    $print_key_key = 1;
    chomp;
    if(/^(ID\s+.+)/){
      $tmp_uniport = $1;
      $print_key = 1;
      $print_key_key = 1;
    }
    if(/^(AC\s+(.+))/){
      $tmp_uniport .= "\n$1";
      my $AC_tmp = $1;
      my @uniport_id_array = ($AC_tmp =~ /[A-Za-z0-9]+/g);
      foreach my $uniport_id (@uniport_id_array){
        if($query_uniport_id_hash{$uniport_id}){
          $annotation_uniport_count++;
          print OUT "$tmp_uniport\n";
          $print_key = 2;
          $uniport_id_DE = $uniport_id;
          last;
        }
      }
      next;
    }
    if($print_key == 2){
      if(/^DE\s+RecName\:/){
        $print_key_key = 2;
        $rec_name_key = 2;
      }
      if(/^OS\s+/){
        $print_key_key = 2;
        $rec_name_key = 1;
      }
      if(/^OX\s+/){
        $print_key_key = 2;
        $rec_name_key = 1;
      }
      if(/^RN\s+/){
        $print_key_key = 2;
        $rec_name_key = 1;
      }
      if(/^RX\s+/){
        $print_key_key = 2;
      }
      if(/^RT\s+/){
        $print_key_key = 2;
      }
      if(/^RL\s+/){
        $print_key_key = 2;
      }
      if(/^CC\s+/){
        $print_key_key = 2;
      }
      if(/^DR\s+/){
        $print_key_key = 2;
      }
      if(/^PE\s+/){
        $print_key_key = 2;
      }
      if(/^FT\s+/){
        $print_key_key = 2;
      }
      if(/^SQ\s+/){
        $print_key_key = 2;
      }
      if(/^\/\//){
        $print_key_key = 2;
      }
    }
    if($print_key_key == 2){
      print OUT "$tmp_all";
    }
    if($rec_name_key == 2){
      $uniport_id_DE_hash{$uniport_id_DE} .= $tmp_all;
    }
  }
  close IN;
  close OUT;

  @query_id_seq_array = keys %query_id_seq_hash;
  @query_id_seq_array = sort{$a <=> $b} @query_id_seq_array;

  open OUT, '>>', $query_annotation_results;
  foreach my $query_id_num (@query_id_seq_array){
    my $query_id = 'QER' . $query_id_num;
    my $query_name = $query_all_name_hash{$query_id_num};
    print OUT "QN   $query_name\n";
    $annotation_query_count++;
    my @mer_query_id_array = ($query_id_seq_hash{$query_id_num} =~ /\w+/g);

    foreach my $mer_query_id (@mer_query_id_array){
      my $actsites = $query_acstites_hash{$mer_query_id};
      my $AMB_actsites = $actsite_AMB_result{$mer_query_id};
      my $domain;
      my $family;
      unless($custom || $annotation_file){
        $domain = $domain_annotation_hash{$mer_query_id};
        $family = $family_annotation_hash{$mer_query_id};
      }
      my $mer_id = $mer_id_hash{$mer_query_id};
      my $uniport_id = $mer_query_uniport_id_hash{$mer_query_id};
      my $evalue = $evalue_hash{$mer_query_id};
      my $score = $score_hash{$mer_query_id};
      my %AMB_actsite_hash;
      my @AMB_actsite_cha_array = ($AMB_actsites =~ /[AMBJOXZ]/g);
      my @AMB_actsite_num_array = ($AMB_actsites =~ /\d+/g);
      my $count_cha = scalar @AMB_actsite_cha_array;
      my $iii;
      for($iii = 0; $iii < $count_cha; $iii++){
         $AMB_actsite_hash{$AMB_actsite_num_array[$iii]} = $AMB_actsite_cha_array[$iii];
      }
      my @actsites_array = ($actsites =~ /\w+/g);
      my @active_array;
      my @metal_array;
      my @binding_array;
      foreach my $actsites_tmp (@actsites_array){
        my $actsites_num = $actsites_tmp;
        $actsites_num =~ s/[A-Za-z]+//g;
        $actsites_num =~ s/\s+//g;
        if($AMB_actsite_hash{$actsites_num} eq 'A'){
          push @active_array, $actsites_tmp;
        }
        elsif($AMB_actsite_hash{$actsites_num} eq 'M'){
          push @metal_array, $actsites_tmp;
        }
        elsif($AMB_actsite_hash{$actsites_num} eq 'B'){
          push @binding_array, $actsites_tmp;
        }
        elsif($AMB_actsite_hash{$actsites_num} eq 'J'){
          push @active_array, $actsites_tmp;
          push @binding_array, $actsites_tmp;
        }
        elsif($AMB_actsite_hash{$actsites_num} eq 'O'){
          push @active_array, $actsites_tmp;
          push @metal_array, $actsites_tmp;
        }
        elsif($AMB_actsite_hash{$actsites_num} eq 'Z'){
          push @binding_array, $actsites_tmp;
          push @metal_array, $actsites_tmp;
        }
        elsif($AMB_actsite_hash{$actsites_num} eq 'X'){
          push @active_array, $actsites_tmp;
          push @metal_array, $actsites_tmp;
          push @binding_array, $actsites_tmp;
        }else{
          print STDERR "ERRPR: actsite changing is error at $mer_query_id.\nPlease redowload software Peptidy.\n";
          exit(1);
        }
      }
      if(@active_array){
        print OUT "AT   ACT_SITE: @active_array\n";
      }
      if(@metal_array){
        print OUT "AT   METAL: @metal_array\n";
      }
      if(@binding_array){
        print OUT "AT   BINDING: @binding_array\n";
      }

      if($domain){
        my @domain_array = ($domain =~ /\d+\-\d+/g);
        my %actsite_in_domain_hash;
        foreach my $domain_tmp (@domain_array){
          my @domain_tmp_array = ($domain_tmp =~ /\d+/g);
          my $domain_start = $domain_tmp_array[0];
          my $domain_end = $domain_tmp_array[1];
          foreach my $actsite_num_tmp (@AMB_actsite_num_array){
            if($actsite_num_tmp > $domain_start && $actsite_num_tmp < $domain_end){
              $actsite_in_domain_hash{$actsite_num_tmp} += 1;
            }
            if($actsite_num_tmp == $domain_start || $actsite_num_tmp == $domain_end){
              $actsite_in_domain_hash{$actsite_num_tmp} += 1;
            }
          }
        }
        if(%actsite_in_domain_hash){
          my @actsite_count_array = keys %actsite_in_domain_hash;
          my $actsite_count = scalar @actsite_count_array;
          my $actsite_count_sum = scalar @AMB_actsite_num_array;
          if($actsite_count == $actsite_count_sum){
            print OUT "DO   Domain: $domain\n";
          }else{
            print OUT "DO   Domain: NA\n";
          }
        }
        unless(%actsite_in_domain_hash){
          print OUT "DO   Domain: NA\n";
        }
      }
      unless($custom){
        unless($annotation_file){
          unless($domain){
            print OUT "DO   Domain: NA\n";
          }
          print OUT "FM   Family: $family\n";
        }
        print OUT "MD   MEROPS_id: $mer_id\n";
      }
      if($custom){
        print OUT "MD   MEROPS_id: $custom_database_seq_name_hash{$mer_id}\n";
      }
      print OUT "EV   E-Value: $evalue\n";
      print OUT "SC   Score: $score\n";
      print OUT "UP   Uniport_id: $uniport_id\n";
      print OUT "$uniport_id_DE_hash{$uniport_id}";
    }
    print OUT "end\n\n\n";
  }

}

sub fetch_annotated_sequence{
  my $in = $_[0];
  my $tag = $_[1];
  my $out;
  my $out1;
  my $out2;
  if($tag eq 'annotate'){
    $out = $workDir  . '/' . $query_file_name . '_annotated_query_seq';
  }elsif($tag eq 'build'){
    $out = $workDir  . '/' . $db_name . '_database_seq';
    $out1 = $workDir  . '/' . $db_name . '_seq_name';
    $out2  = $workDir  . '/' . $db_name . '_seq_length';
  }else{
    print STDERR "ERROR: The \'Peptidy.pl\' is modified. Please redownload it and replace the old one.\n";
    exit(1);
  }
  open IN, $in;
  open OUT, '>>', $out;
  if($tag eq 'build'){
    open OUT1, '>>', $out1;
    open OUT2, '>>', $out2;
  }
  my $print_key = 1;
  my $count_seq_length;
  my $query_id_num_old;
  while(<IN>){
    chomp;
    my $tmp_all = $_;
    if(/^>QER(\d+)$/){
      $print_key = 1;
      my $query_id_num = $1;
      if($query_id_annotation_seq_hash{$query_id_num}){
        if($tag eq 'annotate'){
          print OUT ">$query_all_name_hash{$query_id_num}\n";
        }
        if($tag eq 'build'){
          print OUT ">MER$query_id_num \n";
          print OUT1 "MER$query_id_num#$query_all_name_hash{$query_id_num}\n";
          if($query_id_num_old && $count_seq_length){
            print OUT2 "$query_id_num_old \{peptidase unit\: 1\-$count_seq_length\}\n";
          }
        }
        $print_key = 2;
        $query_id_num_old = $query_id_num;
        $count_seq_length = 0;
        next;
      }
    }
    if($print_key == 2){
      print OUT "$tmp_all\n";
      my $length_tmp = length($tmp_all);
      $count_seq_length += $length_tmp;
    }
  }
  close IN;
  close OUT;
  if($tag eq 'build'){
    print OUT2 "$query_id_num_old \{peptidase unit\: 1\-$count_seq_length\}\n";
    close OUT1;
    close OUT2;
  }
}


#create_actsites_annotation_all();
sub create_actsites_annotation_all{
  $actsites_annotation_all = $query_input_tmp . '_actsite_annotation_all';
  if($actsite_protease_key){
    if(! -e $actsites_annotation_protease){
      print STDERR "ERROR: no \'$actsites_annotation_protease\' file.\n";
      exit(1);
    }
    if($actsite_pepunit_key){
      if(! -e $actsites_annotation_pepunit){
        print STDERR "ERROR: no \'$actsites_annotation_pepunit\' file.\n";
        exit(1);
      }
      my $cmd_combine_all = "cat $actsites_annotation_protease $actsites_annotation_pepunit > $actsites_annotation_all";
      system("$cmd_combine_all")==0 or die("failed to execute command \'$cmd_combine_all\': $!\n");
    }
    unless($actsite_pepunit_key){
      $actsites_annotation_all = $actsites_annotation_protease;
    }
  }
  unless($actsite_protease_key){
    if($actsite_pepunit_key){
      if(! -e $actsites_annotation_pepunit){
        print STDERR "ERROR: no \'$actsites_annotation_pepunit\' file.\n";
        exit(1);
      }
      $actsites_annotation_all = $actsites_annotation_pepunit;
    }
  }

}

sub fetch_actsite_match_query_id{
  my $in = $_[0];
  open IN, $in;
  while(<IN>){
    chomp;
    if(/^(\d+)/){
      my $query_id = $1;
      $actsite_match_query_id_hash{$query_id} += 1;
    }
  }
  close IN;
}

# fetch_actsite_unmatch_query_blastp($blastp_database_out_seq);
sub fetch_actsite_unmatch_query_blastp{

  my $in = $_[0];
  my $out = $query_input_tmp . '_custom_database_blastp';
  $blastp_protease_out = $out;
  my $print_key = 2;
  open IN, $in;
  open OUT, '>>', $out;
  while(<IN>){
    my $tmp_all = $_;
    chomp;
    if(/^Query=\sQER(\d+)/){
      my $query_id = $1;
      $print_key = 1;
      unless($actsite_match_query_id_hash{$query_id}){
        $print_key = 2;
      }
    }
    if(/^\s+Database:\s/){
      $print_key = 2;
    }
    if($print_key == 2){
      print OUT "$tmp_all";
      $switch_next = 1;
    }
  }
  close IN;
  close OUT;
}



########################################## PepTidy peplines ##########################################

#blastp
# MEROPS proteases database: $database_MEROPS_protease
# MEROPS pepunit database: $database_MEROPS_pepunit

my $plus_next;
if($plus){

  $plus_next = 1;
  my $custom_database_protease = $database . '/' . $database_name . '_database_seq';
  Check_fasta_format($query_input);
  
  check_switch_next();
  if($switch_next == 1){
    Reformat_fasta_file($query_input);
  }else{
    print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
    exit(1);
  }

  check_switch_next();
  if($switch_next == 1){
    print STDOUT "Search peptidase homologues.\n";
    Blastp($query_input_tmp, $custom_database_protease, 'database');
  }else{
    print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
    exit(1);
  }

  check_switch_next();
  if($switch_next == 1){
    Check_blastp_input($blastp_database_out);
  }else{
    print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
    exit(1);
  }

  check_switch_next();
  if($switch_next == 1){
    #$blastp_pepunit_out_seq is the output of Fetch_hits_sequence(). it is used for proteases balstp at next actsites_matching step.
    #$count_hits_match is hits_count. it can be used for check the integrality of next actsites_matching step.
    Fetch_hits_sequence($blastp_database_out, $query_input_tmp, $blastp_database_out_seq);
    unless($count_hits_match){
      print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
      exit(1);
    }
  }elsif($switch_next == 2){
    print STDOUT "Ooops, no peptidase homologues is found. Please try to change Blastp searching parameters.\n$usage\n";
    exit(0);
  }else{
    print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
    exit(1);
  }

  delete_file($query_input_tmp);
  $query_input_tmp = $blastp_database_out_seq;


}

if(!defined($custom) || $plus){

  my $delete_file_name;
# Blastp
  
  unless($plus){
    Check_fasta_format($query_input);

    check_switch_next();
    if($switch_next == 1){
      Reformat_fasta_file($query_input);
    }else{
      print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
      exit(1);
    }
  }  

  check_switch_next();
  if($switch_next == 1){
    print STDOUT "Search peptidase homologues.\n";
    Blastp($query_input_tmp, $database_MEROPS_pepunit, 'pepunit');
  }else{
    print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
    exit(1);
  }

  check_switch_next();
  if($switch_next == 1){
    Check_blastp_input($blastp_pepunit_out);
  }else{
    print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
    exit(1);
  }


  check_switch_next();
  if($switch_next == 1){
    #$blastp_pepunit_out_seq is the output of Fetch_hits_sequence(). it is used for proteases balstp at next actsites_matching step.
    #$count_hits_match is hits_count. it can be used for check the integrality of next actsites_matching step.
    Fetch_hits_sequence($blastp_pepunit_out, $query_input_tmp, $blastp_pepunit_out_seq);
    unless($count_hits_match){
      print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
      exit(1);
    }
  }elsif($switch_next == 2){
    print STDOUT "Ooops, no peptidase homologues is found. Please try to change Blastp searching parameters.\n$usage\n";
    exit(0);
  }else{
    print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
    exit(1);
  }

  check_switch_next();
  if($switch_next == 1){
    print STDOUT "Search for Actsites matching.\n";
    Blastp($blastp_pepunit_out_seq, $database_MEROPS_protease, 'protease');
  }else{
    print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
    exit(1);
  }


# Actsites Annotation Peplines at protease lever
  my $key_protease = 1;
  print STDOUT "Actsites Annotation at protease lever:\n";
  check_switch_next();
  if($switch_next == 1){
    Check_blastp_input($blastp_protease_out);
  }else{
    print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
    exit(1);
  }

  check_switch_next();
  if($switch_next == 1){
    Filter_MEROPS_protease($blastp_protease_out, 'protease');
  }elsif($switch_next == 2){
    $switch_next = 1010;
  }else{
    print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
    exit(1);
  }



  if($switch_next == 1){
    Match_actsite($blastp_protease_out, 'protease');
  }
  unless($switch_next){
    unless($switch_ori){
      print STDOUT "Ooops, No matched actsites is found. Please try to change Blastp searching parameters.\n$usage\n";
      exit(0);
    }
    if($switch_ori == 1){
      #$protease_same_query (QER$query_id#MER$mer_id)
      unless($db_build){
        create_merops_annotation($protease_same_query);
        fetch_annotated_sequence($blastp_pepunit_out_seq, 'annotate');
        $delete_file_name = $query_input_tmp . '*';
        delete_file($delete_file_name);
        if($plus){
          delete_file($blastp_database_out);
        }
        if($switch_ori == 1 || $switch_ori == 2 || $switch_ori == 4){
          print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
          exit(1);
        }
        if($switch_ori == 3){
          print STDOUT "\nAnnotated putative peptidases existed in UEROPS: $annotation_ori_count.\nAnnotated putative peptidases not existed in UEROPS: NA.\nPeptidase is done.\n";
          exit(0);
        }
      }
      if($db_build && $merops_switch){
        create_merops_custom_database($protease_same_query);
        fetch_annotated_sequence($blastp_pepunit_out_seq, 'build');
        if($annotation_file){
          create_merops_annotation($protease_same_query);
          fetch_annotated_sequence($blastp_pepunit_out_seq, 'annotate');
          if($switch_ori == 1 || $switch_ori == 2 || $switch_ori == 4){
            print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
            exit(1);
          }
          if($switch_ori == 3){
            print STDOUT "\nAnnotated putative peptidases existed in UEROPS: $annotation_ori_count.\nAnnotated putative peptidases not existed in UEROPS: NA.\nPeptidase is done.\n";
          }
        }
        if($plus){
          delete_file($blastp_database_out);
        }
        $delete_file_name = $query_input_tmp . '*';
        delete_file($delete_file_name);
        exit(0);
      }
    }else{
      print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
      exit(1);
    }
  }

  if($switch_ori){
    if($switch_ori == 1){
      unless($db_build){
        create_merops_annotation($protease_same_query);
        if($switch_ori == 1 || $switch_ori == 2 || $switch_ori == 4){
          print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
          exit(1);
        }
      }
      if($db_build && $merops_switch){
        create_merops_custom_database($protease_same_query);
        if($annotation_file){
          create_merops_annotation($protease_same_query);
          if($switch_ori == 1 || $switch_ori == 2 || $switch_ori == 4){
            print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
            exit(1);
          }
        }
      } 
    }
  }

  check_switch_next();
  if($switch_next == 1){
    Filter_actsite_matching_result($actsite_out_protease, 'protease');
  }else{
    print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
    exit(1);
  }

  check_switch_next();
  if($switch_next == 1010){
    print STDOUT "No matched actsites found at protease lever.\n";
  }
  unless($switch_next == 1010){
    if($switch_next == 1){
      Rearrange_actsite_matching_result($actsite_out_protease_filter, 'protease');
    }else{
      print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
      exit(1);
    }

    check_switch_next();
    if($switch_next == 1){
      Delete_repeat_in_actsite_matching_result($actsite_out_protease_rearrange, 'protease');
    }else{
      print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
      exit(1);
    }

    check_switch_next();
    if($switch_next == 1){
      Separate_actsite_matching_result($actsite_out_protease_delete_repeat, 'protease');
      unless($switch_next || $actsite_result){
        print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
        exit(1);
      }
    }else{
      print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
      exit(1);
    }


    if($switch_next == 1){
      Separate_AA_actsite_matching_result($actsite_out_protease_left1, 'protease');
    }

    if($switch_next == 1){
      Delete_inclusion_actsite_matching_result($actsite_out_protease_left2, 'protease');
    }

    if($actsite_result){
      combine_actsites('protease');
      $actsite_protease_key = 1;
      $actsite_result = 0;
      print STDOUT "Actsites Annotation at protease lever is done. $count_actsites_match putative peptidases are found.\n";
    }

  }


# Actsites Annotation Peplines at pepunit lever

  if($count_actsites_match){
    if($count_hits_match > $count_actsites_match){
      $switch_next = 1010;
    }
  }

  if($switch_next == 1010){
    print STDOUT "Actsites Annotation at pepunit lever:\n";
    $key_protease = 0;
    if($plus){
      $plus_next = 2;
    }
    Fetch_actsites_sequences($protease_blastp_all_query_id ,$blastp_pepunit_out_seq);

    check_switch_next();
    if($switch_next == 1){
      Check_fasta_format($balstp_actsites_left_seq);
    }else{
      print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
      exit(1);
    }

    check_switch_next();
    if($switch_next == 1){
      print STDOUT "Search for Actsites matching.\n";
      Blastp($balstp_actsites_left_seq, $database_MEROPS_pepunit, 'pepunit');
    }else{
      print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
      exit(1);
    }

    check_switch_next();
    if($switch_next == 1){
      Check_blastp_input($blastp_pepunit_out);
    }else{
      print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
      exit(1);
    }

    check_switch_next();
    if($switch_next == 1){
      Match_actsite($blastp_pepunit_out, 'pepunit');
    }elsif($switch_next == 2){
      unless($actsite_protease_key){
        print STDOUT "Ooops, no peptidase homologues is found. Please try to change Blastp searching parameters.\n$usage\n";
        exit(0);
      }
      if($actsite_protease_key){
        print STDOUT "Actsites Annotation is done.\n";
      }
    }else{
      print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
      exit(1);
    }

    check_switch_next();
    unless($switch_next == 2){
      if($switch_next == 1){
        Filter_actsite_matching_result($actsite_out_pepunit, 'pepunit');
      }else{
        print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
        exit(1);
      }

      check_switch_next();
      if($switch_next == 1010){
        print STDOUT "No matched actsites found at pepunit lever.\n";
        unless($actsite_protease_key){
          print STDOUT "\nOoops, no peptidase homologues is found. Please try to change Blastp searching parameters.\n$usage\n";
          exit(0);
        }
        if($actsite_protease_key){
          print STDOUT "Actsites Annotation is done.\n";
        }
      }

      unless($switch_next == 1010){
        if($switch_next == 1){
          Rearrange_actsite_matching_result($actsite_out_pepunit_filter, 'pepunit');
        }

        check_switch_next();
        if($switch_next == 1){
          Delete_repeat_in_actsite_matching_result($actsite_out_pepunit_rearrange, 'pepunit');
        }else{
          print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
          exit(1);
        }

        check_switch_next();
        if($switch_next == 1){
          Separate_actsite_matching_result($actsite_out_pepunit_delete_repeat, 'pepunit');
          unless($switch_next || $actsite_result){
            print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
            exit(1);
          }
        }else{
          print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
          exit(1);
        }

        if($switch_next == 1){
          Separate_AA_actsite_matching_result($actsite_out_pepunit_left1, 'pepunit');
        }

        if($switch_next == 1){
          Delete_inclusion_actsite_matching_result($actsite_out_pepunit_left2, 'pepunit');
        }

        if($actsite_result){
          combine_actsites('pepunit');
          $actsite_pepunit_key = 1;
          $actsite_result = 0;
          print STDOUT "Actsites Annotation at pepunit lever is done. $count_actsites_match putative peptidases are found.\n";
          print STDOUT "Actsites Annotation is done.\n";
        }
      }
    }
  }

  if($db_build || $plus){
    create_actsites_annotation_all();
    Transform_actsites();
    fetch_annotated_sequence($blastp_pepunit_out_seq, 'build');
    if($annotation_file){
      create_query_annotation();
      fetch_annotated_sequence($blastp_pepunit_out_seq, 'annotate');
    }
    if($plus){
      if($plus_next == 2){
        if($actsite_protease_key){
          fetch_actsite_match_query_id($protease_blastp_all_query_id);
        }
        if($actsite_pepunit_key){
          fetch_actsite_match_query_id($pepunit_blastp_all_query_id);
        }
      }
    }
    unless($plus){
      $delete_file_name = $query_input_tmp . '*';
    }
    if($plus){
      if($plus_next == 2){
        $delete_file_name = $query_input_tmp . '_*';
      }
      if($plus_next == 1){
        $delete_file_name = $query_input_tmp . '*';
        delete_file($blastp_database_out);
      }
    }
    delete_file($delete_file_name);
    delete_file($protease_same_query);
    unless($plus){
      MakeBlastdb();
      exit(0);
    }
    if($plus){
      if($plus_next == 1){
        MakeBlastdb();
        exit(0);
      }
      if($plus_next == 2){
        fetch_actsite_unmatch_query_blastp($blastp_database_out);
        delete_file($blastp_database_out);
      }
    }
  }



# domain
  unless($plus){
    $switch_next = 0;
    unless($actsite_protease_key || $actsite_pepunit_key){
      print STDOUT "\nOoops, no peptidase homologues is found. Please try to change Blastp searching parameters.\n$usage\n\n";
      exit(0);
    }

    if($actsite_protease_key || $actsite_pepunit_key){
      $blastp_pepunit_out_first = $query_input_tmp . '_pepunit_blastp'; #add in
      Fetch_domain_blastp_result();
    }

    Transform_actsites();

    if($switch_option){
      if($switch_option == 2 || $switch_option == 3){
        Change_hits_format($blastp_result_for_miss_domain_analysis);
        #$domain_annotation_rank;
        Rank_all_MEROPS_id($blastp_result_for_miss_domain_analysis);
        Splited_domain_assembly($domain_annotation_change);
        $domain_annotation_combine = $domain_annotation_splited . '_combine';
        Domain_region_analysis($domain_annotation_splited, $domain_annotation_combine);
        #$domain_annotation_miss_match_input;
        Create_input_for_domain_analysis($domain_annotation_combine, $domain_annotation_rank);
        $domain_annotation_miss_match_result_tmp = $domain_annotation_miss_match_input . '_tmp_result';
        Domain_region_analysis($domain_annotation_miss_match_input, $domain_annotation_miss_match_result_tmp);
        #$domain_annotation_miss_match_input;
        #delete_file($domain_annotation_miss_match_input);
        create_input_for_miss_domain_analysis();
        $domain_annotation_miss_match_result = $domain_annotation_miss_match_input . '_result';
        Domain_region_analysis($domain_annotation_miss_match_input, $domain_annotation_miss_match_result);
        create_miss_domain_annotation();
      }
      if($switch_option == 1 || $switch_option == 3){
        my $i;
        for($i = 1; $i < $count_max_actsite_repeat; $i++){
          $ii = $i;
          my $in_Filter = "$balstp_result_for_domain_analysis" . "$i";
          Filter_hits_end($in_Filter);
          unless($switch_next){
            create_domain_annotation();
            next;
          }
          Change_hits_format($domain_annotation_hits_end);
          Splited_domain_assembly($domain_annotation_change);
          $domain_annotation_combine = $domain_annotation_splited . '_combine';
          Domain_region_analysis($domain_annotation_splited, $domain_annotation_combine);
          create_input_for_domain_match_analysis($domain_annotation_combine, $actsites_annotation_all);
          $domain_annotation_match_result = $domain_annotation_match_input . '_reslut';
          Domain_region_analysis($domain_annotation_match_input, $domain_annotation_match_result);
          create_domain_annotation();
        }

      }
      if($switch_option == 2){
        $switch_next = 0;
        create_domain_annotation();
      }
    }
    create_query_annotation();
    fetch_annotated_sequence($blastp_pepunit_out_seq, 'annotate');

    $delete_file_name = $query_input_tmp . '*';
    delete_file($delete_file_name);
    delete_file($protease_same_query);
    $delete_file_name = $balstp_result_for_domain_analysis . '*';
    delete_file($delete_file_name);
 
  }

}

if($plus){

  $database_MEROPS_protease = $database . '/' . $database_name . '_database_seq';
  $database_custom_database_seq_name = $database . '/' . $database_name . '_seq_name';
  $MEROPS_protease_length = $database . '/' . $database_name . '_seq_length';
  $MEROPS_protease_actsites_all = $database . '/' . $database_name . '_seq_actsites_all';
  $active_metal_binding_AMB_sum = $database . '/' . $database_name . '_seq_actsites_AMB_sum';
  $mer_uniport_id_sum = $database . '/' . $database_name . '_seq_uniport_annotation_sum';

  Check_peptidy_files($database_MEROPS_protease);
  my $database_MEROPS_protease_blastdb = "$database_MEROPS_protease.phr";
  Check_peptidy_files($database_MEROPS_protease_blastdb);
  $database_MEROPS_protease_blastdb = "$database_MEROPS_protease.pin";
  Check_peptidy_files($database_MEROPS_protease_blastdb);
  $database_MEROPS_protease_blastdb = "$database_MEROPS_protease.psq";
  Check_peptidy_files($database_MEROPS_protease_blastdb);
  Check_peptidy_files($MEROPS_protease_actsites_all);
  Check_peptidy_files($MEROPS_protease_length);
  Check_peptidy_files($active_metal_binding_AMB_sum);
  Check_peptidy_files($mer_uniport_id_sum);
  Check_peptidy_files($uniport_annotation_description);

  


}

if($custom){

  my $delete_file_name;
  # Blastp

  unless($plus){
    Check_fasta_format($query_input);

    check_switch_next();
    if($switch_next == 1){
      Reformat_fasta_file($query_input);
    }else{
      print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
      exit(1);
    }

    check_switch_next();
    if($switch_next == 1){
      print STDOUT "Search for Actsites matching.\n";
      Blastp($query_input_tmp, $database_MEROPS_protease, 'protease');
    }else{
      print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
      exit(1);
    }

  }

# Actsites Annotation Peplines at protease lever
  my $key_protease = 1;
  print STDOUT "Actsites Annotation at protease lever:\n";
  check_switch_next();
  if($switch_next == 1){
    Check_blastp_input($blastp_protease_out);
  }else{
    print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
    exit(1);
  }

  check_switch_next();
  if($switch_next == 1){
    Filter_MEROPS_protease($blastp_protease_out, 'protease');
  }elsif($switch_next == 2){
    print STDOUT "Ooops, no peptidase homologues is found. Please try to change Blastp searching parameters.\n$usage\n";
    exit(0);
  }else{
    print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
    exit(1);
  }

  if($switch_next == 1){
    Match_actsite($blastp_protease_out, 'protease');
  }

  unless($switch_next){
    unless($switch_ori){
      print STDOUT "Ooops, No matched actsites is found. Please try to change Blastp searching parameters.\n$usage\n";
      exit(0);
    }
    if($switch_ori == 1){
      unless($db_build){
        create_merops_annotation($protease_same_query);
        fetch_annotated_sequence($query_input_tmp, 'annotate');
        $delete_file_name = $query_input_tmp . '*';
        delete_file($delete_file_name);
        delete_file($protease_same_query);
        if($switch_ori == 1 || $switch_ori == 2 || $switch_ori == 4){
          print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
          exit(1);
        }
        if($switch_ori == 3){
          print STDOUT "\nAnnotated putative peptidases existed in custom database: $annotation_ori_count.\nAnnotated putative peptidases not existed in custom database: NA.\nPeptidase is done.\n";
          exit(0);
        }
      }
      if($db_build && $merops_switch){
        create_merops_custom_database($protease_same_query);
        fetch_annotated_sequence($query_input_tmp, 'build');
        if($annotation_file){
          create_merops_annotation($protease_same_query);
          fetch_annotated_sequence($query_input_tmp, 'annotate');
          if($switch_ori == 1 || $switch_ori == 2 || $switch_ori == 4){
            print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
            exit(1);
          }
          if($switch_ori == 3){
            print STDOUT "\nAnnotated putative peptidases existed in custom database: $annotation_ori_count.\nAnnotated putative peptidases not existed in custom database: NA.\nPeptidase is done.\n";
          }
        }
        $delete_file_name = $query_input_tmp . '*';
        delete_file($delete_file_name);
        exit(0);
      }
    }else{
      print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
      exit(1);
    }
  }
  if($switch_ori){
    if($switch_ori == 1){
      unless($db_build){
        create_merops_annotation($protease_same_query);
        if($switch_ori == 1 || $switch_ori == 2 || $switch_ori == 4){
          print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
          exit(1);
        }
      }
      if($db_build && $merops_switch){
        create_merops_custom_database($protease_same_query);
        if($annotation_file){
          create_merops_annotation($protease_same_query);
          if($switch_ori == 1 || $switch_ori == 2 || $switch_ori == 4){
            print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
            exit(1);
          }
        }
      }
    }
  }
  check_switch_next();
  if($switch_next == 1){
    Filter_actsite_matching_result($actsite_out_protease, 'protease');
  }else{
    print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
    exit(1);
  }

  check_switch_next();
  if($switch_next == 1010){
    print STDOUT "Ooops, No matched actsites is found. Please try to change Blastp searching parameters.\n$usage\n";
    exit(0);
  }

  unless($switch_next == 1010){
    if($switch_next == 1){
      Rearrange_actsite_matching_result($actsite_out_protease_filter, 'protease');
    }else{
      print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
      exit(1);
    }

    check_switch_next();
    if($switch_next == 1){
      Delete_repeat_in_actsite_matching_result($actsite_out_protease_rearrange, 'protease');
    }else{
      print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
      exit(1);
    }

    check_switch_next();
    if($switch_next == 1){
      Separate_actsite_matching_result($actsite_out_protease_delete_repeat, 'protease');
      unless($switch_next || $actsite_result){
        print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
        exit(1);
      }
    }else{
      print STDERR "Failed. Please redownload the Peptidy.pl and replace the old one.\n";
      exit(1);
    }

    if($switch_next == 1){
      Separate_AA_actsite_matching_result($actsite_out_protease_left1, 'protease');
    }

    if($switch_next == 1){
      Delete_inclusion_actsite_matching_result($actsite_out_protease_left2, 'protease');
    }

    if($actsite_result){
      combine_actsites('protease');
      $actsite_protease_key = 1;
      $actsite_result = 0;
      print STDOUT "Actsites Annotation by custom database is done. $count_actsites_match putative peptidases are found.\n";
    }

  }

  create_actsites_annotation_all();
  Transform_actsites();
  unless($db_build){
    create_query_annotation();  
    fetch_annotated_sequence($query_input_tmp, 'annotate');
  }
  if($db_build){
    fetch_annotated_sequence($query_input_tmp, 'build'); 
    if($annotation_file){
      create_query_annotation();
      fetch_annotated_sequence($query_input_tmp, 'annotate');
    }
  }
  
  $delete_file_name = $query_input_tmp . '*';
  delete_file($delete_file_name);
  delete_file($protease_same_query);

  if($db_build){
    MakeBlastdb();
  }

}





           ###################### Actsites Annotation Peplines ##########################
#Step_1: check the format of Blastp result file

#Step_1.1: if blastp at protease or database lever, filter sequences in database from blastp results

#Step_2: match Actsites

#Step_3: filter Actsite matching results to find out all correct results

#Step_4: rearrrange query actsite by numerical value from small to large

#step_5: delete repeated actsite type in one query_id

#step_6: separate one query_id with single actsite type from one query_id with poly actsite type

#step_7: separate one query_id with single amino acid type from one query_id with poly actsite amino acid type

#step_8: check the inclusion relationship of actsite types within one query_id and delete the inclued actsite type

#step_9: combine all actsites annotation results

#step_10: fetch unmatched sequences
#if($count_hits_match > $count_actsites_match){
  #Fetch_actsites_sequences($protease_blastp_all_query_id ,$blastp_pepunit_out_seq);
#}

#if($switch_next == 1){
#  Check_fasta_format($balstp_actsites_left_seq);
#}

#if($switch_next == 1){
  #Blastp($balstp_actsites_left_seq, $database_MEROPS_pepunit, 'pepunit');
#}



