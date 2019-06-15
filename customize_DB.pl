#! /usr/bin/perl

eval 'exec perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shele

###################################################################
#                                                                 #
# customize_DB.pl                                                 #
# creating custom database from MEROPS database                   #
#                                                                 #
# Authors: Yi Jiang                                               #
#                                                                 #
# Contact   : yijiang@peptidy.cn                                  #
# Bug report: bugs@peptidy.cn                                     # 
#                                                                 #
# Release date: April 15th, 2019                                  #
#                                                                 #
# This script is under the Artistic Licence                       #
# https://opensource.org/licenses/artistic-license-2.0            #
#                                                                 #
# Usage:                                                          #
# perl customize_DB.pl [option] --out                             #
#                                                                 #
###################################################################

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

customize_DB.pl is to create custom database from MEROPS database.

SYNOPSIS

perl customize_DB.pl [option] --out

   --help          To print this help message.
                   --help
   --in            path to directory of Peptidy MEROPS database.
                   --in <String>
                   Example [--in /path/to/Peptidy_MEROPS_database_dir]
                   Default = './'
   --organism      To include organism's proteases in custom database.
                   --organism <String>
                   Example [--organism organism1 --organism organism2 ... --organism organismn]
                   Default = 'Fungi'
                   Please set the organism name from the following options:
                   ===================================================================================
                   Archaea(Euryarchaeota, Crenarchaeota); 
                   Bacteria(Proteobacteria, Actinobacteria, Bacteroidetes, Cyanobacteria, Firmicutes);
                   Eukaryota(Fungi, Viridiplantae, Metazoa);
                   Virus.
                   ===================================================================================
                   such as [--organism Viridiplantae --Bacteria]
                   '--organism Archaea Euryarchaeota' is invalid. Because Archaea contains Euryarchaeota.
                   '--organism archaea' is invalid. Because the fist letter should be uppercase.

   --de_organism   To exclude organism's proteases in custom database.
                   --de_organism <String>
                   Example [--de_organism organism1 --de_organism organism2 ... --de_organism organismn]
   --family        To include MEROPS family's proteases in custom database.
                   --family <String>
                   Example [--family family1 --family family2 ... --family familyn]
                   such as [--family A --family M]
                   Please set the family name from nine options "A, S, M, C, N, G, U, T, P".
   --de_family     To exclude MEROPS family's proteases in custom database.
                   --de_family <String>
                   Example [--de_family family1 --de_family family2 ... --de_family familyn]
   --code          To include MEROPS code's proteases in custom database.
                   --code <String>
                   Example [--code code1 --code code2 ... --code coden]
                   such as [--code A11 --code S08 --code M24]
                   Code name like 'A11.01' will not be accepted.
   --de_code       To exclude MEROPS code's proteases in custom database.
                   --de_code <String>
                   Example [--de_code code1 --de_code code2 ... --de_code coden]
   --out           custom database name and path to it.
                   --out <String>
                   Example [--out /path/to/out_dir/out_name]

   !!!NOTE: The priority of '--de_organism', '--de_family' and '--de_code' are higher than '--organism', '--family' and '--code'.

perl customize_DB.pl [options] --out

DESCRIPTION

  Example:
    
    perl customize_DB.pl [options] --out

ENDUSAGE

my $help;                                  # to print usage
my $currentDir = cwd();                    # working superdirectory where programme is called from
my $workDir;                               # in the working directory results and temporary files are stored
my @organism_array;                        # to include organism's proteases in custom database
my %organism_hash;                         # organism name
my @de_organism_array;                     # to exclude organism's proteases in custom database
my %de_organism_hash;                      # organism name
my @family_array;                          # to include MEROPS family's proteases in custom database
my %family_hash;                           # family name
my @de_family_array;                       # to exclude MEROPS family's proteases in custom database
my %de_family_hash;                        # family name
my @code_array;                            # to include MEROPS code's proteases in custom database
my %code_hash;                             # code name
my @de_code_array;                         # to exclude MEROPS code's proteases in custom database
my %de_code_hash;                          # code name
my $output;                                # custom database name and path to it
my $input;                                 # protein sequences input in fasta fromat
my $database_name;                         # custom database name
my %mer_id_output_hash;
my $count_check;
my $MEROPS_protease_actsites_all;
my $MEROPS_protease_AMB_sum;
my $MEROPS_protease_length;
my $MEROPS_protease_uniport_sum;
my $MEROPS_database_protease;
my $MEROPS_protease_tax_code;
my $database_seq_name;
my $database_seq_actsite;
my $database_seq_AMB;
my $database_seq_length;
my $database_seq_uniprot;
my $database_seq_seq;
my %organism_sum_hash = ('Archaea'        => 'SUM1',
                         'Euryarchaeota'  => 'SUB1',
                         'Crenarchaeota'  => 'SUB1',
                         'Bacteria'       => 'SUM2',
                         'Proteobacteria' => 'SUB2',
                         'Actinobacteria' => 'SUB2',
                         'Bacteroidetes'  => 'SUB2',
                         'Cyanobacteria'  => 'SUB2',
                         'Firmicutes'     => 'SUB2',
                         'Eukaryota'      => 'SUM3',
                         'Fungi'          => 'SUB2',
                         'Viridiplantae'  => 'SUB2',
                         'Metazoa'        => 'SUB2',
                         'Virus'          => 'SUM4');

my %family_sum_hash = ('A' => 1,
                       'S' => 1,
                       'M' => 1,
                       'C' => 1,
                       'N' => 1,
                       'G' => 1,
                       'U' => 1,
                       'T' => 1,
                       'P' => 1);

if(@ARGV==0){
  print "$usage\n";
  exit(0);
}

GetOptions( 'help!'                         => \$help,
            'in=s'                          => \$input,
            'organism=s'                    => \@organism_array,
            'de_organism=s'                 => \@de_organism_array,
            'family=s'                      => \@family_array,
            'de_family=s'                   => \@de_family_array,
            'code=s'                        => \@code_array,
            'de_code=s'                     => \@de_code_array,
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
  $MEROPS_protease_actsites_all = $input . '/MEROPS_protease_actsites_all';
  $MEROPS_protease_AMB_sum = $input . '/MEROPS_protease_actsites_AMB_sum';
  $MEROPS_protease_length = $input . '/MEROPS_protease_length';
  $MEROPS_protease_uniport_sum = $input . '/MEROPS_uniport_annotation_sum';
  $MEROPS_protease_tax_code = $input . '/MEROPS_protease_tax_code';
  $MEROPS_database_protease = $input . '/MEROPS_database/protease_lib';
}else{
  $MEROPS_protease_actsites_all = $directory . '/MEROPS_protease_actsites_all';
  $MEROPS_protease_AMB_sum = $directory . '/MEROPS_protease_actsites_AMB_sum';
  $MEROPS_protease_length = $directory . '/MEROPS_protease_length';
  $MEROPS_protease_uniport_sum = $directory . '/MEROPS_uniport_annotation_sum';
  $MEROPS_protease_tax_code = $directory . '/MEROPS_protease_tax_code';
  $MEROPS_database_protease = $directory . '/MEROPS_database/protease_lib';
}

Check_peptidy_files($MEROPS_protease_actsites_all);
Check_peptidy_files($MEROPS_protease_AMB_sum);
Check_peptidy_files($MEROPS_protease_length);
Check_peptidy_files($MEROPS_protease_uniport_sum);
Check_peptidy_files($MEROPS_database_protease);
Check_peptidy_files($MEROPS_protease_tax_code);

if($output){
  if($output =~ /^(.+)\/([^\/]+)$/){
    $workDir = $1;
    $database_name = $2;
    $workDir = rel2abs($workDir);
  }elsif($output =~ /^([^\/]+)$/){
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
  my $workDir_tmp = "$workDir/$database_name" . '_db';
  my $cmd_mkdir_tmp = "mkdir $workDir_tmp";
  system("$cmd_mkdir_tmp")==0 or die("failed to execute command \'$cmd_mkdir_tmp\': $!\n");
  $workDir = $workDir_tmp;
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


if(@organism_array){
  my %sum1;
  my %sum2;
  my %sum3;
  foreach my $organism (@organism_array){
    if($organism_sum_hash{$organism}){
      $organism_hash{$organism} += 1;
      if($organism_sum_hash{$organism} eq 'SUM1' or $organism_sum_hash{$organism} eq 'SUB1'){
        $sum1{$organism_hash{$organism}} += 1;
      }elsif($organism_sum_hash{$organism} eq 'SUM2' or $organism_sum_hash{$organism} eq 'SUB2'){
        $sum2{$organism_hash{$organism}} += 1;
      }elsif($organism_sum_hash{$organism} eq 'SUM3' or $organism_sum_hash{$organism} eq 'SUB3'){
        $sum3{$organism_hash{$organism}} += 1;
      }elsif($organism_sum_hash{$organism} eq 'SUM4'){
      }else{
        print STDERR "FATAL ERROR: Please redownload \'customize_DB.pl\'.\n";
        exit(1);
      }
    }else{
      print STDERR "ERROR: The vaild arguments of '--organism' are Archaea\nEuryarchaeota\nCrenarchaeota\nBacteria\nProteobacteria\nActinobacteria\nBacteroidetes\nCyanobacteria\nFirmicutes\nEukaryota\nFungi\nViridiplantae\nMetazoa\nVirus\n";
      exit(1);
    }
  }
  if(%sum1){
    my $count = scalar(keys %sum1);
    if($count > 1){
      print STDERR "ERROR: if set 'Archaea', there is no need to set 'Euryarchaeota' or 'Crenarchaeota' in option '--organism'.\n";
      exit(1);
    }
  }
  if(%sum2){
    my $count = scalar(keys %sum2);
    if($count > 1){
      print STDERR "ERROR: if set 'Bacteria', there is no need to set 'Proteobacteria', 'Actinobacteria', 'Bacteroidetes', 'Cyanobacteria' or 'Firmicutes' in option '--organism'.\n";
      exit(1);
    }
  }
  if(%sum3){
    my $count = scalar(keys %sum3);
    if($count > 1){
      print STDERR "ERROR: if set 'Eukaryota', there is no need to set 'Fungi', 'Viridiplantae' or 'Metazoa' in option '--organism'.\n";
      exit(1);
    }
  }
}
  
if(@de_organism_array){
  my %sum1;
  my %sum2;
  my %sum3;
  foreach my $organism (@de_organism_array){
    if($organism_sum_hash{$organism}){
      $de_organism_hash{$organism} += 1;
      if($organism_sum_hash{$organism} eq 'SUM1' or $organism_sum_hash{$organism} eq 'SUB1'){
        $sum1{$de_organism_hash{$organism}} += 1;
      }elsif($organism_sum_hash{$organism} eq 'SUM2' or $organism_sum_hash{$organism} eq 'SUB2'){
        $sum2{$de_organism_hash{$organism}} += 1;
      }elsif($organism_sum_hash{$organism} eq 'SUM3' or $organism_sum_hash{$organism} eq 'SUB3'){
        $sum3{$de_organism_hash{$organism}} += 1;
      }elsif($organism_sum_hash{$organism} eq 'SUM4'){
      }else{
        print STDERR "FATAL ERROR: Please redownload \'customize_DB.pl\'.\n";
        exit(1);
      }
    }else{
      print STDERR "ERROR: The vaild arguments of '--de_organism' are Archaea\nEuryarchaeota\nCrenarchaeota\nBacteria\nProteobacteria\nActinobacteria\nBacteroidetes\nCyanobacteria\nFirmicutes\nEukaryota\nFungi\nViridiplantae\nMetazoa\nVirus\n";
      exit(1);
    }
  }
  if(%sum1){
    my $count = scalar(keys %sum1);
    if($count > 1){
      print STDERR "ERROR: if set 'Archaea', there is no need to set 'Euryarchaeota' or 'Crenarchaeota' in option '--de_organism'.\n";
      exit(1);
    }
  }
  if(%sum2){
    my $count = scalar(keys %sum2);
    if($count > 1){
      print STDERR "ERROR: if set 'Bacteria', there is no need to 'Proteobacteria', 'Actinobacteria', 'Bacteroidetes', 'Cyanobacteria' or 'Firmicutes' in option '--de_organism'.\n";
      exit(1);
    }
  }
  if(%sum3){
    my $count = scalar(keys %sum3);
    if($count > 1){
      print STDERR "ERROR: if set 'Eukaryota', there is no need to 'Fungi', 'Viridiplantae' or 'Metazoa' in option '--de_organism'.\n";
      exit(1);
    }
  }
}

if(@organism_array && @de_organism_array){
  my %organism_mix_hash;
  my %sum1;
  my %sum2;
  my %sum3;
  my %sub1;
  my %sub2;
  my %sub3;
  foreach my $organism (keys %organism_hash){
    $organism_mix_hash{$organism} = 1;
    if($organism_sum_hash{$organism} eq 'SUM1'){
       $sum1{$organism_hash{$organism}} = 1;
    }
    if($organism_sum_hash{$organism} eq 'SUB1'){
      $sub1{$organism_hash{$organism}} += 1;
    }
    if($organism_sum_hash{$organism} eq 'SUM2'){
      $sum2{$organism_hash{$organism}} = 1;
    }
    if($organism_sum_hash{$organism} eq 'SUB2'){
      $sub2{$organism_hash{$organism}} += 1;
    }
    if($organism_sum_hash{$organism} eq 'SUM3'){
      $sum3{$organism_hash{$organism}} = 1;
    }
    if($organism_sum_hash{$organism} eq 'SUB3'){
      $sub3{$organism_hash{$organism}} += 1;
    }
  }
  foreach my $organism (keys %de_organism_hash){
    $organism_mix_hash{$organism} += 1;
    if($organism_mix_hash{$organism} > 1){
      print STDERR "ERROR: '$organism' can only be set in '--organism' or '--de_organism'.\n";
      exit(1);
    }
    if($organism_sum_hash{$organism} eq 'SUM1'){
       $sum1{$organism_hash{$organism}} += 1;
       $sub1{$organism_hash{$organism}} += 1;
    }
    if($organism_sum_hash{$organism} eq 'SUM2'){       
      $sum2{$organism_hash{$organism}} += 1;
      $sub2{$organism_hash{$organism}} += 1;
    }
    if($organism_sum_hash{$organism} eq 'SUM3'){
      $sum3{$organism_hash{$organism}} += 1;
      $sub3{$organism_hash{$organism}} += 1;
    }
  }
  if(%sum1){
    if($sum1{'SUM1'} > 1){
      print STDERR "ERROR: if set 'Archaea' in '--de_organism', there is no need to 'Archaea' in option '--organism'.\n";
      exit(1);
    }
  }
  if(%sub1){
    my $count = scalar(keys %sub1);
    if($count > 1){
      print STDERR "ERROR: if set 'Archaea' in '--de_organism', there is no need to 'Euryarchaeota' or 'Crenarchaeota' in option '--organism'.\n";
      exit(1);
    }
  }
  if(%sum2){
    if($sum2{'SUM2'} > 1){
      print STDERR "ERROR: if set 'Bacteria' in '--de_organism', there is no need to 'Bacteria' in option '--organism'.\n";
      exit(1);
    }
  }
  if(%sub2){
    my $count = scalar(keys %sub2);
    if($count > 1){
      print STDERR "ERROR: if set 'Bacteria' in '--de_organism', there is no need to set 'Proteobacteria', 'Actinobacteria', 'Bacteroidetes', 'Cyanobacteria' or 'Firmicutes' in option '--organism'.\n";

      exit(1);
    }
  }
  if(%sum3){
    if($sum3{'SUM3'} > 1){
      print STDERR "ERROR: if set 'Eukaryota' in '--de_organism', there is no need to set 'Eukaryota' in option '--organism'.\n";
      exit(1);
    }
  }
  if(%sub3){
    my $count = scalar(keys %sub3);
    if($count > 1){
      print STDERR "ERROR: if set 'Eukaryota' in '--de_organism', there is no need to 'Fungi', 'Viridiplantae' or 'Metazoa' in option '--organism'.\n";
      exit(1);
    }
  }

}

unless(@organism_array || @de_organism_array || @family_array || @de_family_array || @code_array || @de_code_array){
  %organism_hash = ('Fungi' => 1);
}


my %family_tmp_hash;
my %code_tmp_hash;
if(@family_array){
  foreach my $family (@family_array){
    unless($family_sum_hash{$family}){
      print STDERR "ERROR: The vaild arguments of '--family' are A\nS\nM\nC\nN\nG\nU\nT\nP\n";
      exit(1);
    }
    $family_tmp_hash{$family} += 1;
    $family_hash{$family} += 1;
  }
}

if(@code_array){
  foreach my $code (@code_array){
    if($code =~ /^([ASMCNGUTP])(\d+)$/){
      my $family = $1;
      my $number = $2;
      my $code;
      next if $family_hash{$family};
      if($number){
        $code = $family . '0' . $number  if length($number) == 1;
        $code = $family . $number if length($number) > 1;
        $code_hash{$code} += 1;
        $family_tmp_hash{$family} += 1;
      }else{
        print STDERR "ERROR: The vaild format of argument of '--code' is one capital with numbers such as A1 or A12.\nThe number is a non-zero integer.\n";
        exit(1);
      }
      
    }else{
      print STDERR "ERROR: The vaild format of argument of '--code' is one capital with numbers such as A1 or A12.\nThe capital is one of these 'ASMCNGUTP'. The number is a non-zero integer.\n";
      exit(1); 
    }
  }
}

if(@de_family_array){
  foreach my $family (@de_family_array){
    unless($family_sum_hash{$family}){
      print STDERR "ERROR: The vaild arguments of '--de_family' are A\nS\nM\nC\nN\nG\nU\nT\nP\n";
      exit(1);
    }
    if($family_tmp_hash{$family}){
      print STDERR "ERROR: if set '$family' in '--de_family', there is no need to set '$family' in '--family' or '--code'.\n";
      exit(1);
    }
    $de_family_hash{$family} += 1;
  }
}

if(@de_code_array){
  foreach my $code (@de_code_array){
    if($code =~ /^([ASMCNGUTP])(\d+)$/){
      my $family = $1;
      my $number = $2;
      my $code;
      next if $de_family_hash{$family};
      if($number){
        $code = $family . '0' . $number  if length($number) == 1;
        $code = $family . $number if length($number) > 1;
        if($code_hash{$code}){
          print STDERR "ERROR: Arguments of '--code' and '--de_code' should be different.\n";
          exit(1);
        }
        $de_code_hash{$code} += 1;
      }else{
        print STDERR "ERROR: The vaild format of argument of '--de_code' is one capital with numbers such as A1 or A12.\nThe number is a non-zero integer.\n";
        exit(1);
      }

    }else{
      print STDERR "ERROR: The vaild format of argument of '--de_code' is one capital with numbers such as A1 or A12.\nThe capital is one of these 'ASMCNGUTP'. The number is a non-zero integer.\n";
      exit(1);
    }
  }
}


#check whether files exist
sub Check_peptidy_files{
  my $file = $_[0];
  unless($file){
    print STDERR "FATAL ERROR: Please redownload \'customize_DB.pl\'.\n";
    exit(1);
  }
  if(! -e $file || -d $file){
    print STDERR "ERROR: No file \'$file\'. Please check.\n";
    exit(1);
  }
}

print STDOUT "All options setting are ok.\n";

############################ Main functon ##################################

read_mer_id($MEROPS_protease_tax_code);
if(%mer_id_output_hash){
  my $count_mer_id = scalar(keys %mer_id_output_hash);
  write_actsite_like($MEROPS_protease_actsites_all, $database_seq_actsite);
  if($count_check){
    unless($count_check == $count_mer_id){
      print STDERR "ERROR: File '$MEROPS_protease_actsites_all' or file '$MEROPS_protease_tax_code' is modified. Please redownload them.\n";
      exit(1);
    }
  }else{
    print STDERR "ERROR: File '$MEROPS_protease_actsites_all' or file is modified. Please redownload them.\n";
    exit(1);
  }
  write_actsite_like($MEROPS_protease_AMB_sum, $database_seq_AMB);
  if($count_check){
    unless($count_check == $count_mer_id){
      print STDERR "ERROR: File '$MEROPS_protease_AMB_sum' or file '$MEROPS_protease_tax_code' is modified. Please redownload them.\n";
      exit(1);
    }
  }else{
    print STDERR "ERROR: File '$MEROPS_protease_actsites_all' or file is modified. Please redownload them.\n";
    exit(1);
  }
  write_actsite_like($MEROPS_protease_uniport_sum, $database_seq_uniprot);
  if($count_check){
    unless($count_check == $count_mer_id){
      print STDERR "ERROR: File '$MEROPS_protease_uniport_sum' or file '$MEROPS_protease_tax_code' is modified. Please redownload them.\n";
      exit(1);
    }
  }else{
    print STDERR "ERROR: File '$MEROPS_protease_actsites_all' or file is modified. Please redownload them.\n";
    exit(1);
  }
  write_length($MEROPS_protease_length, $database_seq_length);
  if($count_check){
    unless($count_check == $count_mer_id){
      print STDERR "ERROR: File '$MEROPS_protease_length' or file '$MEROPS_protease_tax_code' is modified. Please redownload them.\n";
      exit(1);
    }
  }else{
    print STDERR "ERROR: File '$MEROPS_protease_actsites_all' or file is modified. Please redownload them.\n";
    exit(1);
  }
  write_seq_name($MEROPS_database_protease);
  if($count_check){
    unless($count_check == $count_mer_id){
      print STDERR "ERROR: File '$MEROPS_database_protease' or file '$MEROPS_protease_tax_code' is modified. Please redownload them.\n";
      exit(1);
    }
    my $cmd_makeblastdb = "makeblastdb -in $database_seq_seq -dbtype prot -out $database_seq_seq";
    system("$cmd_makeblastdb")==0 or die("failed to execute command \'$cmd_makeblastdb\': $!\n");
  }else{
    print STDERR "ERROR: File '$MEROPS_protease_actsites_all' or file is modified. Please redownload them.\n";
    exit(1);
  }
}else{
  print STDOUT "No matching protease. Please change your search criteria.\n";
  exit(0);
}

sub read_mer_id{
  my $in = $_[0];
  open IN, $in;
  while(<IN>){
    chomp;
    if(/^(\d+)\scode:(([\w])\d+)\..+\srank1:(.+)\srank2:(.+)$/){
      my $mer_id = $1;
      my $code = $2;
      my $family = $3;
      my $oganism1 = $4;
      my $oganism2 = $5;
      next if $de_organism_hash{$oganism1};
      next if $de_organism_hash{$oganism2};
      next if $de_family_hash{$family};
      next if $de_code_hash{$code};
      $mer_id_output_hash{$mer_id} += 1 if $organism_hash{$oganism1};
      $mer_id_output_hash{$mer_id} += 1 if $organism_hash{$oganism2};
      $mer_id_output_hash{$mer_id} += 1 if $family_hash{$family};
      $mer_id_output_hash{$mer_id} += 1 if $code_hash{$code};
    }
  }
  close IN;
}

sub write_actsite_like{
  my $in = $_[0];
  my $out = $_[1];
  $count_check = 0;
  open IN, $in;
  open OUT, '>>', $out;
  while(<IN>){
    chomp;
    if(/^(\d+)#/){
      my $mer_id = $1;
      if($mer_id_output_hash{$mer_id}){
        print OUT "$_\n";
        $count_check++;
      }
    }else{
      print STDERR "ERROR: File '$in' is modified. Please redownload it.\n";
      exit(1);
    }
  }
  close IN;
  close OUT;
}

sub write_length{
  my $in = $_[0];
  my $out = $_[1];
  $count_check = 0;
  open IN, $in;
  open OUT, '>>', $out;
  while(<IN>){
    chomp;
    if(/^(\d+)\s/){
      my $mer_id = $1;
      if($mer_id_output_hash{$mer_id}){
        print OUT "$_\n";
        $count_check++;
      }
    }else{
      print STDERR "ERROR: File '$in' is modified. Please redownload it.\n";
      exit(1);
    }
  }
  close IN;
  close OUT;
}

sub write_seq_name{
  my $in = $_[0];
  my $print_key = 1;
  $count_check = 0;
  open IN, $in;
  open OUT1, '>>', $database_seq_name;
  open OUT2, '>>', $database_seq_seq;
  while(<IN>){
    chomp;
    if(/^>/){
      if(/^>(MER(\d+)\s.+)$/){
        my $tmp = $1;
        my $mer_id = $2;
        $print_key = 1;
        if($mer_id_output_hash{$mer_id}){
          $count_check++;
          print OUT1 "MER$mer_id#$tmp\n";
          $print_key = 2;
        }
      }else{
        print STDERR "ERROR: File '$in' is modified. Please redownload it.\n";
        exit(1);
      }
    }
    if($print_key == 2){
      print OUT2 "$_\n";
    }
  }
  close IN;
  close OUT1;
  close OUT2;
}
