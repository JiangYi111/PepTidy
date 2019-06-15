#! /usr/bin/perl

eval 'exec perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shele;


###########################################################################
#                                                                         #
# make_DBs.pl                                                             #
# Proteins and residues databases creating                                #
#                                                                         #
# Authors: Yi Jiang                                                       #
#                                                                         #
# Contact: yijiang@peptidy.cn                                             #
# Bug report: bugs@peptidy.cn                                             # 
#                                                                         #
# Release date: April 15th, 2019                                          #
#                                                                         #
# This script is under the Artistic Licence                               #
# https://opensource.org/licenses/artistic-license-2.0                    #
#                                                                         #
# Usage:                                                                  #
# perl make_DBs.pl --seq --site -dbname --out                             #
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

make_DBs.pl is to build proteins and residues databases which are used 
for mining proteins by PepTidy.pl.

SYNOPSIS

perl make_DBs.pl --seq --site  --dbname --out

   --help     To print this help message.
              --help
   --seq      Name of protein sequences file and path to it.
              --seq <String>
              Example [--seq /path/to/sequences_file_name]
   --site     Name of residues file and path to it.
              Example [--site /path/to/residues_file_name]
   --dbname   The name of database.
              --dbname <String>
              Example [--dbname name]
              !NOTE: The name must start with letter, and only letter, 
                     integer and underline are allowed.
   --out      The directory path used to save output databases files.
              --out <String>
              Example [--out /path/to/out]
              Default = './'

DESCRIPTION

  Example:
    
    perl make_DBs.pl --seq --site  --dbname --out

ENDUSAGE

my $help;                                  # to print usage
my $currentDir = cwd();                    # working superdirectory where programme is called from
my $seq_file;                              # The file of protein sequences
my $site_file;                             # The file of residues
my $database;                              # the directory containing all custom database files
my $database_name;                         # custom database name
my $database_output;                       # in database output directory results and temporary files are stored
my $database_seq_name;
my $database_seq_actsite;
my $database_seq_AMB;
my $database_seq_length;
my $database_seq_uniprot;
my $database_seq_seq;
my %uni_AC_hash;
my %seq_hash;
my %seq_name_hash;




if(@ARGV==0){
  print "$usage\n";
  exit(0);
}

GetOptions( 'help!'                 => \$help,
            'seq=s'                 => \$seq_file,
            'site=s'                => \$site_file,
            'dbname=s'              => \$database_name,
            'out=s'                 => \$database_output) or die("$!");

if($help){
  print $usage;
  exit(0);
}

######################## make some regular checks ########################

my $cmd_check = "makeblastdb -h";
system("$cmd_check")==0 or die("Please install ncbi Blast+ first.\n");

if($database_name){
  unless($database_name =~ /^[A-Za-z][A-Za-z\d\_]*$/){
    print STDERR "ERROR: Please set option '--dbname' with right format.\n$usage\n";
    exit(1)
  }
}
unless($database_name){
  print STDERR "ERROR: Please set option '--dbname'.\n$usage\n";
  exit(1);
}

if($seq_file){
  if(! -e $seq_file || -d $seq_file){
    print STDERR "ERROR: No file \'$seq_file\'. Please reset option '--seq'.\n$usage\n";
    exit(1);
  }
}
unless($seq_file){
  print STDERR "ERROR: Please set option '--seq'.\n$usage\n";
  exit(1);
}
Check_fasta_format($seq_file);


if($site_file){
  if(! -e $site_file || -d $site_file){
    print STDERR "ERROR: No file \'$site_file\'. Please reset option '--seq'.\n$usage\n";
    exit(1);
  }
}
unless($site_file){
  print STDERR "ERROR: Please set option '--site'.\n$usage\n";
  exit(1);
}
Check_residues_format($site_file);

unless($database_output){
  $database_output = $currentDir;
}
if($database_output){
  if(! -d $database_output){
    print STDERR "ERROR: Directory \'$database_output\' does not exist. Please reset option '--out'.\n$usage\n";
    exit(1);
  }
  if(! -w $database_output){
    print STDERR "ERROR: Directory \'$database_output\' does not the right to write. Please check.\n$usage\n";
    exit(1);
  }
  $database_output = rel2abs($database_output);
  my $directory_tmp = $database_output;
  $database_output = $database_output . '/' . $database_name . '_peptidy_DB';
  if(! -e $database_output){
    my $cmd_mkdir = "mkdir $database_output";
    system("$cmd_mkdir")==0 or die("Creating a new folder \'$database_output\' for storing databases files failed. Please check.\n");
  }else{
    print STDERR "ERROR: Please make sure no file or directory named \'$database_output\' in directory \'$directory_tmp\'.\n";
    exit(1);
  }
  $database_seq_name = $database_output . '/' . $database_name . '_seq_name';
  if(-e $database_seq_name){
    print STDERR "ERROR: Output database name '$database_name' exist. Please change it.\n";
    exit(1);
  }
  $database_seq_actsite = $database_output . '/' . $database_name . '_seq_actsites_all';
  if(-e $database_seq_actsite){
    print STDERR "ERROR: Output database name '$database_name' exist. Please change it.\n";
    exit(1);
  }
  $database_seq_AMB = $database_output . '/' . $database_name . '_seq_actsites_AMB_sum';
  if(-e $database_seq_AMB){
    print STDERR "ERROR: Output database name '$database_name' exist. Please change it.\n";
    exit(1);
  }
  $database_seq_length = $database_output . '/' . $database_name . '_seq_length';
  if(-e $database_seq_length){
    print STDERR "ERROR: Output database name '$database_name' exist. Please change it.\n";
    exit(1);
  }
  $database_seq_uniprot = $database_output . '/' . $database_name . '_seq_uniport_annotation_sum';
  if(-e $database_seq_uniprot){
    print STDERR "ERROR: Output database name '$database_name' exist. Please change it.\n";
    exit(1);
  }
  $database_seq_seq = $database_output . '/' . $database_name . '_database_seq';
  if(-e $database_seq_seq){
    print STDERR "ERROR: Output database name '$database_name' exist. Please change it.\n";
    exit(1);
  }
}

#Check_residues_format
sub Check_residues_format{
  print STDOUT "Check residues format...\n";
  my $in = $_[0];
  my $count_line;

  open IN, $in;
  while(<IN>){
    $count_line++;
    if(/^Entry\t/){
      unless($count_line == 1){
        print STDERR "ERROR: Please read the section of 'Exploitation' in PepTidy_Userguide.pdf.\nAnd redownload the residues file.\n";
        exit(1);
      }
    }
    if(/^([^\t]+)\t+ACT_SITE \d+ \d+/){
      my $uni_AC = $1;
      $uni_AC_hash{$uni_AC} += 1;
    }
    if(/^([^\t]+)\t+METAL \d+ \d+/){
      my $uni_AC = $1;
      $uni_AC_hash{$uni_AC} += 1;
    }
    if(/^([^\t]+)\t+BINDING \d+ \d+/){
      my $uni_AC = $1;
      $uni_AC_hash{$uni_AC} += 1;
    }
  }
  close IN;
  unless(%uni_AC_hash){
    print STDERR "ERROR: Please read the section of 'Exploitation' in PepTidy_Userguide.pdf.\nAnd redownload the residues file.\n";
    exit(1);
  }
  print STDOUT "Done.\n";
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
          print STDERR "ERROR: Error at $count_start_tmp line of $in. No sequences.\nPlease read the section of 'Exploitation' in PepTidy_Userguide.pdf.\nAnd redownload the sequences file in fasta format.\n";
          exit(1);
        }
      }
      $count_start = $count_line;
      my $name_tmp = $1;
      $query_name_hash{$name_tmp} += 1;
      if($query_name_hash{$name_tmp} > 1){
        print STDERR "WARNING: query name: \'$name_tmp\' repeat in $in.\n";
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


############################ Main function ##################################
read_seq_write_seq_name();
read_uni_AC();
my $cmd_mkbldb = "makeblastdb -in $database_seq_seq -dbtype prot -out $database_seq_seq";
system("$cmd_mkbldb")==0 or die("Error occured when running makeblastdb.\n");
print "\nDONE! PepTidy databases of '$database_name' is builded.\nFiles are stored in '$database_output'\n";


sub read_uni_AC{
  my $in = $site_file;

  open IN, $in;
  open OUT1, '>>', $database_seq_length;
  open OUT2, '>>', $database_seq_actsite;
  open OUT3, '>>', $database_seq_AMB;
  open OUT4, '>>', $database_seq_uniprot;

  while(<IN>){
    chomp;
    next if /^Entry\t/;
    my $tmp = $_;
    $tmp =~ s/\r//g;
    my %site_hash;
    if($tmp =~ /^([^\t]+)\t\s*ACT_SITE \d+ \d+/){
      my $Uni_AC = $1;
      my $seq_length = length $seq_hash{$Uni_AC};
      print OUT1 "$seq_name_hash{$Uni_AC} \{peptidase unit: 1\-$seq_length\}\n";
      print OUT4 "$seq_name_hash{$Uni_AC}#$Uni_AC\n";
      my @active_array = ($tmp =~ /ACT_SITE \d+ \d+/g);
      my @metal_array = ($tmp =~ /METAL \d+ \d+/g);
      my @binding_array = ($tmp =~ /BINDING \d+ \d+/g);
      if(@active_array){
        foreach my $site_tmp (@active_array){
          my @site_tmp_array = ($site_tmp =~ /\d+/g);
          if($site_tmp_array[0] == $site_tmp_array[1]){
            $site_hash{$site_tmp_array[0]} .= 'A';
          }
          elsif($site_tmp_array[0] < $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[0]; $i <= $site_tmp_array[1]; $i++){
               $site_hash{$i} .= 'A';
            }
          }elsif($site_tmp_array[0] > $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[1]; $i <= $site_tmp_array[0]; $i++){
               $site_hash{$i} .= 'A';
            }
          }else{
            print STDERR "ERROR1: Residue information has bug at \'$Uni_AC\'.\nIf the residue file is download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
            exit(1);
          }
        }
      }
      if(@metal_array){
        foreach my $site_tmp (@metal_array){
          my @site_tmp_array = ($site_tmp =~ /\d+/g);
          if($site_tmp_array[0] == $site_tmp_array[1]){
            $site_hash{$site_tmp_array[0]} .= 'M';
          }
          elsif($site_tmp_array[0] < $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[0]; $i <= $site_tmp_array[1]; $i++){
               $site_hash{$i} .= 'M';
            }
          }elsif($site_tmp_array[0] > $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[1]; $i <= $site_tmp_array[0]; $i++){
               $site_hash{$i} .= 'M';
            }
          }else{
            print STDERR "ERROR2: Residue information has bug at \'$Uni_AC\'.\nIf the residue file is download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
            exit(1);
          }
        }
      }
      if(@binding_array){
        foreach my $site_tmp (@binding_array){
          my @site_tmp_array = ($site_tmp =~ /\d+/g);
          if($site_tmp_array[0] == $site_tmp_array[1]){
            $site_hash{$site_tmp_array[0]} .= 'B';
          }
          elsif($site_tmp_array[0] < $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[0]; $i <= $site_tmp_array[1]; $i++){
               $site_hash{$i} .= 'B';
            }
          }elsif($site_tmp_array[0] > $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[1]; $i <= $site_tmp_array[0]; $i++){
               $site_hash{$i} .= 'B';
            }
          }else{
            print STDERR "ERROR3: Residue information has bug at \'$Uni_AC\'.\nIf the residue file is download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
            exit(1);
          }
        }

      }
      if(%site_hash){
        my @site_sort_array = sort{$a <=> $b} (keys %site_hash);
        my $count_sum = scalar @site_sort_array;
        my $count_tmp;
        print OUT2 "$seq_name_hash{$Uni_AC}#\"";
        print OUT3 "$seq_name_hash{$Uni_AC}#\"";
        foreach my $site_tmp (@site_sort_array){
          $count_tmp++;
          if($site_hash{$site_tmp} =~ /^[A]+$/){
            print OUT3 "A$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[M]+$/){
            print OUT3 "M$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[B]+$/){
            print OUT3 "B$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[A]+[M]+$/){
            print OUT3 "O$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[A]+[B]+$/){
            print OUT3 "J$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[A]+[M]+[B]+$/){
            print OUT3 "X$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[M]+[B]+$/){
            print OUT3 "Z$site_tmp";
          }
          else{
            print STDERR "ERROR4: Residue information has bug at \'$Uni_AC\'.\nIf the residue file is download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
            exit(1);
          }
          my $site_substr = $site_tmp -1;
          my $site_residue = substr($seq_hash{$Uni_AC},$site_substr,1);
          unless($site_residue){
            print STDERR "ERROR: Residue and sequence file has bug at \'$Uni_AC\'.\nIf files are download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
            exit(1);
          }
          print OUT2 "$site_residue$site_tmp";
 
          if($count_tmp < $count_sum){
            print OUT3 "\, ";
            print OUT2 "\, ";
          }
          if($count_tmp == $count_sum){
            print OUT3 "\"\,\n";
            print OUT2 "\"\,\n";
          }
        }
      }else{
        print STDERR "ERROR5: Residue information has bug at \'$Uni_AC\'.\nIf the residue file is download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
        exit(1);
      }
      next;
    }
    if($tmp =~ /^([^\t]+)\t\s*METAL \d+ \d+/){
      my $Uni_AC = $1;
      my $seq_length = length $seq_hash{$Uni_AC};
      print OUT1 "$seq_name_hash{$Uni_AC} \{peptidase unit: 1\-$seq_length\}\n";
      print OUT4 "$seq_name_hash{$Uni_AC}#$Uni_AC\n";
      my @active_array = ($tmp =~ /ACT_SITE \d+ \d+/g);
      my @metal_array = ($tmp =~ /METAL \d+ \d+/g);
      my @binding_array = ($tmp =~ /BINDING \d+ \d+/g);
      if(@active_array){
        foreach my $site_tmp (@active_array){
          my @site_tmp_array = ($site_tmp =~ /\d+/g);
          if($site_tmp_array[0] == $site_tmp_array[1]){
            $site_hash{$site_tmp_array[0]} .= 'A';
          }
          elsif($site_tmp_array[0] < $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[0]; $i <= $site_tmp_array[1]; $i++){
               $site_hash{$i} .= 'A';
            }
          }elsif($site_tmp_array[0] > $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[1]; $i <= $site_tmp_array[0]; $i++){
               $site_hash{$i} .= 'A';
            }
          }else{
            print STDERR "ERROR6: Residue information has bug at \'$Uni_AC\'.\nIf the residue file is download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
            exit(1);
          }
        }
      }
      if(@metal_array){
        foreach my $site_tmp (@metal_array){
          my @site_tmp_array = ($site_tmp =~ /\d+/g);
          if($site_tmp_array[0] == $site_tmp_array[1]){
            $site_hash{$site_tmp_array[0]} .= 'M';
          }
          elsif($site_tmp_array[0] < $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[0]; $i <= $site_tmp_array[1]; $i++){
               $site_hash{$i} .= 'M';
            }
          }elsif($site_tmp_array[0] > $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[1]; $i <= $site_tmp_array[0]; $i++){
               $site_hash{$i} .= 'M';
            }
          }else{
            print STDERR "ERROR7: Residue information has bug at \'$Uni_AC\'.\nIf the residue file is download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
            exit(1);
          }
        }
      }
      if(@binding_array){
        foreach my $site_tmp (@binding_array){
          my @site_tmp_array = ($site_tmp =~ /\d+/g);
          if($site_tmp_array[0] == $site_tmp_array[1]){
            $site_hash{$site_tmp_array[0]} .= 'B';
          }
          elsif($site_tmp_array[0] < $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[0]; $i <= $site_tmp_array[1]; $i++){
               $site_hash{$i} .= 'B';
            }
          }elsif($site_tmp_array[0] > $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[1]; $i <= $site_tmp_array[0]; $i++){
               $site_hash{$i} .= 'B';
            }
          }else{
            print STDERR "ERROR8: Residue information has bug at \'$Uni_AC\'.\nIf the residue file is download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
            exit(1);
          }
        }

      }
      if(%site_hash){
        my @site_sort_array = sort{$a <=> $b} (keys %site_hash);
        my $count_sum = scalar @site_sort_array;
        my $count_tmp;
        print OUT2 "$seq_name_hash{$Uni_AC}#\"";
        print OUT3 "$seq_name_hash{$Uni_AC}#\"";
        foreach my $site_tmp (@site_sort_array){
          $count_tmp++;
          if($site_hash{$site_tmp} =~ /^[A]+$/){
            print OUT3 "A$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[M]+$/){
            print OUT3 "M$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[B]+$/){
            print OUT3 "B$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[A]+[M]+$/){
            print OUT3 "O$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[A]+[B]+$/){
            print OUT3 "J$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[A]+[M]+[B]+$/){
            print OUT3 "X$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[M]+[B]+$/){
            print OUT3 "Z$site_tmp";
          }
          else{
            print STDERR "ERROR9: Residue information has bug at \'$Uni_AC\'.\nIf the residue file is download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
            exit(1);
          }
          my $site_substr = $site_tmp -1;
          my $site_residue = substr($seq_hash{$Uni_AC},$site_substr,1);
          unless($site_residue){
            print STDERR "ERROR: Residue and sequence file has bug at \'$Uni_AC\'.\nIf files are download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
            exit(1);
          }
          print OUT2 "$site_residue$site_tmp";

          if($count_tmp < $count_sum){
            print OUT3 "\, ";
            print OUT2 "\, ";
          }
          if($count_tmp == $count_sum){
            print OUT3 "\"\,\n";
            print OUT2 "\"\,\n";
          }
        }
      }else{
        print STDERR "ERROR10: Residue information has bug at \'$Uni_AC\'.\nIf the residue file is download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
        exit(1);
      }
      next;
    }
    if($tmp =~ /^([^\t]+)\t\s*BINDING \d+ \d+/){
      my $Uni_AC = $1;
      my $seq_length = length $seq_hash{$Uni_AC};
      print OUT1 "$seq_name_hash{$Uni_AC} \{peptidase unit: 1\-$seq_length\}\n";
      print OUT4 "$seq_name_hash{$Uni_AC}#$Uni_AC\n";
      my @active_array = ($tmp =~ /ACT_SITE \d+ \d+/g);
      my @metal_array = ($tmp =~ /METAL \d+ \d+/g);
      my @binding_array = ($tmp =~ /BINDING \d+ \d+/g);
      if(@active_array){
        foreach my $site_tmp (@active_array){
          my @site_tmp_array = ($site_tmp =~ /\d+/g);
          if($site_tmp_array[0] == $site_tmp_array[1]){
            $site_hash{$site_tmp_array[0]} .= 'A';
          }
          elsif($site_tmp_array[0] < $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[0]; $i <= $site_tmp_array[1]; $i++){
               $site_hash{$i} .= 'A';
            }
          }elsif($site_tmp_array[0] > $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[1]; $i <= $site_tmp_array[0]; $i++){
               $site_hash{$i} .= 'A';
            }
          }else{
            print STDERR "ERROR11: Residue information has bug at \'$Uni_AC\'.\nIf the residue file is download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
            exit(1);
          }
        }
      }
      if(@metal_array){
        foreach my $site_tmp (@metal_array){
          my @site_tmp_array = ($site_tmp =~ /\d+/g);
          if($site_tmp_array[0] == $site_tmp_array[1]){
            $site_hash{$site_tmp_array[0]} .= 'M';
          }
          elsif($site_tmp_array[0] < $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[0]; $i <= $site_tmp_array[1]; $i++){
               $site_hash{$i} .= 'M';
            }
          }elsif($site_tmp_array[0] > $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[1]; $i <= $site_tmp_array[0]; $i++){
               $site_hash{$i} .= 'M';
            }
          }else{
            print STDERR "ERROR12: Residue information has bug at \'$Uni_AC\'.\nIf the residue file is download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
            exit(1);
          }
        }
      }
      if(@binding_array){
        foreach my $site_tmp (@binding_array){
          my @site_tmp_array = ($site_tmp =~ /\d+/g);
          if($site_tmp_array[0] == $site_tmp_array[1]){
            $site_hash{$site_tmp_array[0]} .= 'B';
          }
          elsif($site_tmp_array[0] < $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[0]; $i <= $site_tmp_array[1]; $i++){
               $site_hash{$i} .= 'B';
            }
          }elsif($site_tmp_array[0] > $site_tmp_array[1]){
            my $i;
            for($i = $site_tmp_array[1]; $i <= $site_tmp_array[0]; $i++){
               $site_hash{$i} .= 'B';
            }
          }else{
            print STDERR "ERROR13: Residue information has bug at \'$Uni_AC\'.\nIf the residue file is download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
            exit(1);
          }
        }

      }
      if(%site_hash){
        my @site_sort_array = sort{$a <=> $b} (keys %site_hash);
        my $count_sum = scalar @site_sort_array;
        my $count_tmp;
        print OUT2 "$seq_name_hash{$Uni_AC}#\"";
        print OUT3 "$seq_name_hash{$Uni_AC}#\"";
        foreach my $site_tmp (@site_sort_array){
          $count_tmp++;
          if($site_hash{$site_tmp} =~ /^[A]+$/){
            print OUT3 "A$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[M]+$/){
            print OUT3 "M$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[B]+$/){
            print OUT3 "B$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[A]+[M]+$/){
            print OUT3 "O$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[A]+[B]+$/){
            print OUT3 "J$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[A]+[M]+[B]+$/){
            print OUT3 "X$site_tmp";
          }
          elsif($site_hash{$site_tmp} =~ /^[M]+[B]+$/){
            print OUT3 "Z$site_tmp";
          }
          else{
            print STDERR "ERROR14: Residue information has bug at \'$Uni_AC\'.\nIf the residue file is download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
            exit(1);
          }
          my $site_substr = $site_tmp -1;
          my $site_residue = substr($seq_hash{$Uni_AC},$site_substr,1);
          unless($site_residue){
            print STDERR "ERROR: Residue and sequence file has bug at \'$Uni_AC\'.\nIf files are download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
            exit(1);
          }
          print OUT2 "$site_residue$site_tmp";

          if($count_tmp < $count_sum){
            print OUT3 "\, ";
            print OUT2 "\, ";
          }
          if($count_tmp == $count_sum){
            print OUT3 "\"\,\n";
            print OUT2 "\"\,\n";
          }
        }
      }else{
        print STDERR "ERROR15: Residue information has bug at \'$Uni_AC\'.\nIf the residue file is download from Uniprot with out edition, please contact PepTidy team.\(yijiang\@peptidy.cn\)\n";
        exit(1);
      }
    }
  }
  close IN;
  close OUT1;
  close OUT2;
  close OUT3;
  close OUT4;
}


sub read_seq_write_seq_name{
  my $in = $seq_file;
  my $Uni_AC;
  my %reapt_hash;
  my $mer_count;
  my $name_tmp;
  open IN, $in;
  open OUT1, '>>', $database_seq_seq;
  open OUT2, '>>', $database_seq_name;
  while(<IN>){
    chomp;
    my $tmp = $_;
    $tmp =~ s/\r//g;
    $name_tmp = $tmp;
    $name_tmp =~ s/\>//;
    if(/^>sp\|([^\|]+)\|/){
      $Uni_AC = $1;
      $reapt_hash{$Uni_AC} += 1;
      if($reapt_hash{$Uni_AC} == 1){
        if($uni_AC_hash{$Uni_AC}){
          $mer_count++;
          my $mer_id = $mer_count;
          print OUT1 ">MER$mer_id \n";
          $seq_name_hash{$Uni_AC} = $mer_id;
          print OUT2 "MER$mer_id#$name_tmp\n";
        }
      }
      next;
    }
    next if $reapt_hash{$Uni_AC} > 1;
    next unless $uni_AC_hash{$Uni_AC};
    print OUT1 "$tmp\n";
    unless($Uni_AC){
      print STDERR "The format of sequences files is wrong. If it is downloaded from Uniprot without edition, please email to yijiang\@peptidy.cn\n";
      exit(1);
    }
    $seq_hash{$Uni_AC} .= $tmp;
  }
  close IN;
  close OUT1;
  close OUT2;
}


