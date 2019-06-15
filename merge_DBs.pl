#! /usr/bin/perl

eval 'exec perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shele;


###########################################################################
#                                                                         #
# merge_DBs.pl                                                            #
# custom databases concatenating                                          #
#                                                                         #
# Authors: Yi Jiang                                                       #
#                                                                         #
# Contact   : yijiang@peptidy.cn                                          #
# Bug report: bugs@peptidy.cn                                             # 
#                                                                         #
# Release date: April 15th, 2019                                          #
#                                                                         #
# This script is under the Artistic Licence                               #
# https://opensource.org/licenses/artistic-license-2.0                    #
#                                                                         #
# Usage:                                                                  #
# perl merge_DBs.pl [option] --db ... --db --out                          #
#                                                                         #
###########################################################################

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

merge_DBs.pl is to concatenate custom databases.

SYNOPSIS

perl merge_DBs.pl [option] --db ... --db --out

   --help    To print this help message.
             --help
   --db      The custom database name and the directory path to custom database files.
             --db <String>
             Example [--db /path/to/db1_dir/db1 --db /path/to/db2_dir/db2 ... --db /path/to/dbn_dir/dbn]
   --out     The concatenated custom database name and the path to the directory for saving the database.
             --out <String>
             Example [--out /path/to/out_dir/out_name]
   --nr      To delete redundant sequences in custom databases.
             --nr

DESCRIPTION

  Example:
    
    perl merge_DBs.pl [option] --db ... --db --out

ENDUSAGE

my $help;                                  # to print usage
my $currentDir = cwd();                    # working superdirectory where programme is called from
my $workDir;                               # in the working directory results and temporary files are stored
my $non_redundancy;                        # to delete redundant sequences
my @database_array;                        # containing all custom databases names and directorys;
my $database;                              # the directory containing all custom database files;
my $database_name;                         # custom database name
my $database_output;                       # the concatenated custom database name and directory
my $custom_database_seq_name;
my $custom_database_seq_actsite;
my $custom_database_seq_AMB;
my $custom_database_seq_length;
my $custom_database_seq_uniprot;
my $custom_database_seq_seq;
my $database_seq_name;
my $database_seq_actsite;
my $database_seq_AMB;
my $database_seq_length;
my $database_seq_uniprot;
my $database_seq_seq;
my %repeat_mer_id_hash;
my %seq_name_hash;
my %seq_actsite_hash;
my %seq_AMB_hash;
my %seq_length_hash;
my %seq_uniprot_hash;
my %seq_seq_hash;
my %output_mer_id_hash;
my $count_sum;
my $count;



if(@ARGV==0){
  print "$usage\n";
  exit(0);
}

GetOptions( 'help!'                         => \$help,
            'nr!'                           => \$non_redundancy,
            'db=s'                          => \@database_array,
            'out=s'                         => \$database_output) or die("$!");

if($help){
  print $usage;
  exit(0);
}

######################## make some regular checks ########################

my $cmd_check = "makeblastdb -h";
system("$cmd_check")==0 or die("Please install ncbi Blast+ first.\n");

if(@database_array){
  my @database_array_tmp = @database_array;
  foreach $database (@database_array_tmp){
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
      print STDERR "ERROR: Option \'--db\' needs a argument with correct format.\n$usage\n";
      exit(1);
    }
    print "Custom database directory: \'$database\'\nCustom database name: \'$database_name\'\n";
    $custom_database_seq_name = $database . '/' . $database_name . '_seq_name';
    Check_peptidy_files($custom_database_seq_name);
    $custom_database_seq_actsite = $database . '/' . $database_name . '_seq_actsites_all';
    Check_peptidy_files($custom_database_seq_actsite);
    $custom_database_seq_AMB = $database . '/' . $database_name . '_seq_actsites_AMB_sum';
    Check_peptidy_files($custom_database_seq_AMB);
    $custom_database_seq_length = $database . '/' . $database_name . '_seq_length';
    Check_peptidy_files($custom_database_seq_length);
    $custom_database_seq_uniprot = $database . '/' . $database_name . '_seq_uniport_annotation_sum';
    Check_peptidy_files($custom_database_seq_uniprot);
    $custom_database_seq_seq = $database . '/' . $database_name . '_database_seq';
    Check_peptidy_files($custom_database_seq_seq);
  }
}else{
  print STDERR "ERROR: Please set option \'--db\'\n$usage\n";
  exit(1);
}

if($database_output){
  if($database_output =~ /^(.+)\/([^\/]+)$/){
    $workDir = $1;
    $database_name = $2;
    $workDir = rel2abs($workDir);
  }elsif($database_output =~ /^([^\/]+)$/){
    $database_name = $1;
    $workDir = rel2abs($currentDir);
  }else{
    print STDERR "ERROR: Option \'--out\' needs a argument with correct format.\n$usage\n";
    exit(1);
  }
  if(! -e $workDir){
    print STDERR "ERROR: Output directory \'$workDir\' does not exist. Please check.\n";
    exit(1);
  }
  if(! -d $workDir){
    print STDERR "ERROR: Output directory \'$workDir\' does not exist. Please check.\n";
    exit(1);
  }
  if(! -w $workDir){
    print STDERR "ERROR: Output directory \'$workDir\' does not have the right to write. Please check.\n";
    exit(1);
  }
  $workDir = $workDir . '/' . $database_name . '_db';
  my $cmd_mkdir = "mkdir $workDir";
  system("$cmd_mkdir")==0 or die("failed to execute command \'$cmd_mkdir\': $!\n");;

  $database_seq_name = $workDir . '/' . $database_name . '_seq_name';
  if(-e $database_seq_name){
    print STDERR "ERROR: Output database name '$database_name' exist. Please change.\n";
    exit(1);
  }
  $database_seq_actsite = $workDir . '/' . $database_name . '_seq_actsites_all';
  if(-e $database_seq_actsite){
    print STDERR "ERROR: Output database name '$database_name' exist. Please change.\n";
    exit(1);
  }
  $database_seq_AMB = $workDir . '/' . $database_name . '_seq_actsites_AMB_sum';
  if(-e $database_seq_AMB){
    print STDERR "ERROR: Output database name '$database_name' exist. Please change.\n";
    exit(1);
  }
  $database_seq_length = $workDir . '/' . $database_name . '_seq_length';
  if(-e $database_seq_length){
    print STDERR "ERROR: Output database name '$database_name' exist. Please change.\n";
    exit(1);
  }
  $database_seq_uniprot = $workDir . '/' . $database_name . '_seq_uniport_annotation_sum';
  if(-e $database_seq_uniprot){
    print STDERR "ERROR: Output database name '$database_name' exist. Please change.\n";
    exit(1);
  }
  $database_seq_seq = $workDir . '/' . $database_name . '_database_seq';
  if(-e $database_seq_seq){
    print STDERR "ERROR: Output database name '$database_name' exist. Please change.\n";
    exit(1);
  }
}else{
  print STDERR "ERROR: Please set option \'--out\'\n$usage\n";
  exit(1);
}

#check whether files exist
sub Check_peptidy_files{
  my $file = $_[0];
  unless($file){
    print STDERR "FATAL ERROR: Please redownload \'merge_DBs.pl\'.\n";
    exit(1);
  }
  if(! -e $file || -d $file){
    print STDERR "ERROR: No file \'$file\'. Please check.\n";
    exit(1);
  }
}

#check fasta format
sub Check_fasta_format{
  print STDOUT "Check fasta format ...\n";
  my $in = $_[0];
  my $count_line;
  my $count_start;
  my %query_name_hash;
  my $count_add = 1;
  open IN, $in;
  while(<IN>){
    $_ =~ s/\r//g;
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
}

print STDOUT "All options setting are ok.\n";


############################ Main functon ##################################

concatenate_custom_databases();

sub concatenate_custom_databases{
 
  my $seq_key = 1;
 
  foreach $database (@database_array){
    if($database =~ /^(.+)\/(\w+)$/){
      $database = $1;
      $database_name = $2;
      $database = rel2abs($database);
    }
    $custom_database_seq_name = $database . '/' . $database_name . '_seq_name';
    $custom_database_seq_actsite = $database . '/' . $database_name . '_seq_actsites_all';
    $custom_database_seq_AMB = $database . '/' . $database_name . '_seq_actsites_AMB_sum';
    $custom_database_seq_length = $database . '/' . $database_name . '_seq_length';
    $custom_database_seq_uniprot = $database . '/' . $database_name . '_seq_uniport_annotation_sum';
    $custom_database_seq_seq = $database . '/' . $database_name . '_database_seq';

    %output_mer_id_hash = ();
    my $count_ckeck;
    #read sequnecs
    read_seq($custom_database_seq_seq);
    if(%output_mer_id_hash){
      $count_ckeck = scalar(keys %output_mer_id_hash);
    }else{
      print STDOUT "WARNING: No sequence in custom database of '$database/$database_name'. Please check.\n";
      next;
    }
    $seq_key = 2;
    #read actsites
    read_actsite($custom_database_seq_actsite);
    if($count){
      unless($count == $count_ckeck){
        print STDERR "ERROR: Some sequences are missing in file '$custom_database_seq_actsite'.\n Please make sure the custom database not be modified.\n";
        exit(1);
      }
    }else{
      print STDERR "ERROR: No sequence in file '$custom_database_seq_actsite'.\n Please make sure the custom database not be modified.\n";
      exit(1);
    }
    #read AMB
    read_AMB($custom_database_seq_AMB);
    if($count){
      unless($count == $count_ckeck){
        print STDERR "ERROR: Some sequences are missing in file '$custom_database_seq_AMB'.\n Please make sure the custom database not be modified.\n";
        exit(1);
      }
    }else{
      print STDERR "ERROR: No sequence in file '$custom_database_seq_AMB'.\n Please make sure the custom database not be modified.\n";
      exit(1);
    }
    #read uniprot
    read_uniprot($custom_database_seq_uniprot);
    if($count){
      unless($count == $count_ckeck){
        print STDERR "ERROR: Some sequences are missing in file '$custom_database_seq_uniprot'.\n Please make sure the custom database not be modified.\n";
        exit(1);
      }
    }else{
      print STDERR "ERROR: No sequence in file '$custom_database_seq_uniprot'.\n Please make sure the custom database not be modified.\n";
      exit(1);
    }
    #read length
    read_length($custom_database_seq_length);
    if($count){
      unless($count == $count_ckeck){
        print STDERR "ERROR: Some sequences are missing in file '$custom_database_seq_length'.\n Please make sure the custom database not be modified.\n";
        exit(1);
      }
    }else{
      print STDERR "ERROR: No sequence in file '$custom_database_seq_length'.\n Please make sure the custom database not be modified.\n";
      exit(1);
    }
    #read name
    read_name($custom_database_seq_name);
    if($count){
      unless($count == $count_ckeck){
        print STDERR "ERROR: Some sequences are missing in file '$custom_database_seq_name'.\n Please make sure the custom database not be modified.\n";
        exit(1);
      }
    }else{
      print STDERR "ERROR: No sequence in file '$custom_database_seq_name'.\n Please make sure the custom database not be modified.\n";
      exit(1);
    }
  }

  if($seq_key == 1){
    print "No sequence in custom databases. Please check.\n";
    exit(0);
  }

  if($non_redundancy){
    de_redundancy();
  }

  $count_sum = 0;
  #write concatenated custom database
  my @mer_id_array = sort(keys %seq_uniprot_hash);
  open OUT1, '>>', $database_seq_actsite;
  open OUT2, '>>', $database_seq_AMB;
  open OUT3, '>>', $database_seq_uniprot;
  open OUT4, '>>', $database_seq_length;
  open OUT5, '>>', $database_seq_name;
  open OUT6, '>>', $database_seq_seq;
  foreach my $mer_id (@mer_id_array){
    next if $repeat_mer_id_hash{$mer_id};
    $count_sum++;
    print OUT1 "$mer_id$seq_actsite_hash{$mer_id}\n";
    print OUT2 "$mer_id$seq_AMB_hash{$mer_id}\n";
    print OUT3 "$mer_id$seq_uniprot_hash{$mer_id}\n";
    print OUT4 "$mer_id$seq_length_hash{$mer_id}\n";
    print OUT5 "$mer_id$seq_name_hash{$mer_id}\n";
    print OUT6 ">$mer_id \n$seq_seq_hash{$mer_id}\n";
  }
  close OUT1;
  close OUT2;
  close OUT3;
  close OUT4;
  close OUT5;
  close OUT6;
  my $cmd_mkblastdb = "makeblastdb -in $database_seq_seq -dbtype prot -out $database_seq_seq";
  system("$cmd_mkblastdb")==0 or die("failed to execute command \'$cmd_mkblastdb\': $!\n");
  print "\n$count_sum sequences are in new custom database '$workDir'.\n";

}

sub read_uniprot{
  my $in = $_[0];
  $count = 0;
  open IN, $in;
  while(<IN>){
    chomp;
    if(/^(\d+)(#.+)$/){
      $count++;
      my $mer_id = $1;
      my $tmp = $2;
      unless($output_mer_id_hash{$mer_id}){
        print STDERR "ERROR: Pleas make sure the files in database '$database/$database_name' are not modified.\n";
        exit(1);
      }
      $seq_uniprot_hash{$output_mer_id_hash{$mer_id}} = $tmp;
    }else{
      print STDERR "ERROR: wrong format of '$in'.\nPleas make sure the file is not modified.\n";
      exit(1);
    }
  }
  close IN;
}


sub read_AMB{
  my $in = $_[0];
  $count = 0;
  open IN, $in;
  while(<IN>){
    chomp;
    if(/^(\d+)(#.+)$/){
      $count++;
      my $mer_id = $1;
      my $tmp = $2;
      unless($output_mer_id_hash{$mer_id}){
        print STDERR "ERROR: Pleas make sure the files in database '$database/$database_name' are not modified.\n";
        exit(1);
      }
      $seq_AMB_hash{$output_mer_id_hash{$mer_id}} = $tmp;
    }else{
      print STDERR "ERROR: wrong format of '$in'.\nPleas make sure the file is not modified.\n";
      exit(1);
    }
  }
  close IN;
}

sub read_actsite{
  my $in = $_[0];
  $count = 0;
  open IN, $in;
  while(<IN>){
    chomp;
    if(/^(\d+)(#.+)$/){
      $count++;
      my $mer_id = $1;
      my $tmp = $2;
      unless($output_mer_id_hash{$mer_id}){
        print STDERR "ERROR: Pleas make sure the files in database '$database/$database_name' are not modified.\n";
        exit(1);
      }
      $seq_actsite_hash{$output_mer_id_hash{$mer_id}} = $tmp;
    }else{
      print STDERR "ERROR: wrong format of '$in'.\nPleas make sure the file is not modified.\n";
      exit(1);
    }
  }
  close IN;
}

sub read_seq{
  my $in = $_[0];
  $count = 0;
  open IN, $in;
  while(<IN>){
    chomp;
    if(/^>MER(\d+)\s/){
      my $mer_id = $1;
      $count_sum++;
      $count++;
      $output_mer_id_hash{$mer_id} = $count_sum;
      next;
    }
    if(/^(\w+)/){
      my $seq_tmp = $1;
      $seq_seq_hash{$count_sum} .= $seq_tmp;
    }
  }
  close IN;
}

sub read_length{
  my $in = $_[0];
  $count = 0;
  open IN, $in;
  while(<IN>){
    chomp;
    if(/^(\d+)(\s.+)$/){
      $count++;
      my $mer_id = $1;
      my $tmp = $2;
      unless($output_mer_id_hash{$mer_id}){
        print STDERR "ERROR: Pleas make sure the files in database '$database/$database_name' are not modified.\n";
        exit(1);
      }
      $seq_length_hash{$output_mer_id_hash{$mer_id}} = $tmp;
    }else{
      print STDERR "ERROR: wrong format of '$in'.\nPleas make sure the file is not modified.\n";
      exit(1);
    }
  }
  close IN;
}

sub read_name{
  my $in = $_[0];
  $count = 0;
  open IN, $in;
  while(<IN>){
    chomp;
    if(/^MER(\d+)(#.+)$/){
      $count++;
      my $mer_id = $1;
      my $tmp = $2;
      unless($output_mer_id_hash{$mer_id}){
        print STDERR "ERROR: Pleas make sure the files in database '$database/$database_name' are not modified.\n";
        exit(1);
      }
      $seq_name_hash{$output_mer_id_hash{$mer_id}} = $tmp;
    }else{
      print STDERR "ERROR: wrong format of '$in'.\nPleas make sure the file is not modified.\n";
      exit(1);
    }
  }
  close IN;
}

sub de_redundancy{
  my %seq_hash;
  my @mer_id_seq_array = (keys %seq_seq_hash);
  foreach my $mer_id (@mer_id_seq_array){
    $seq_hash{$seq_seq_hash{$mer_id}} += 1;
    if($seq_hash{$seq_seq_hash{$mer_id}} > 1){
      $repeat_mer_id_hash{$mer_id} += 1;
    }
  }
}

