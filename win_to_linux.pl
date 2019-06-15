#! /usr/bin/perl

eval 'exec perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shele;

use strict;
use warnings;
use Getopt::Long;
use Cwd;
use File::Spec::Functions qw(rel2abs);
use File::Basename qw(dirname basename);

#########################################################
#                                                       #
# Authors: Yi Jiang                                     #
#                                                       #
# Contact: yijiang@peptidy.cn                           #
#                                                       #
# This script is under the Artistic Licence             #
# https://opensource.org/licenses/artistic-license-2.0  #
#                                                       #
# Usage:                                                #
# perl linux_to_win.pl -in -out                         #
#                                                       #
# #######################################################

my $usage = <<'ENDUSAGE';

win_to_linux.pl is used to convert windows-type file to linux-type file.

SYNOPSIS

perl win_to_linux.pl -in -out

   --help    To print this help message.
             --help
   --in      input file with win-type
             --in <String>
             Example [--in \path\to\input_file]
   --out     name of output file
             --out <String>
             Example [--in output_file_name]

ENDUSAGE

if(@ARGV==0){
  print "$usage\n";
  exit(0);
}

my $in;
my $out;
my $help;

GetOptions( 'help!'                         => \$help,
            'in=s'                          => \$in,
            'out=s'                         => \$out) or die("$!");

if($help){
  print $usage;
  exit(0);
}

unless($out){
  print STDERR "ERROR: No output file name.\n$usage\n";
  exit(1);
}

unless($in){
  print STDERR "ERROR: Please specify the input file.\n$usage\n";
  exit(1);
}

unless(-e $in){
  print STDERR "ERROR: File \'$in\' does not exist. Please check.\n$usage\n";
  exit(1);
}

if(-d $in){
  print STDERR "ERROR: \'$in\' is a directory. Please specify a file as input.\n$usage\n";
  exit(1);
}

if(-d $out){
  print STDERR "ERROR: \'$out\' is a directory. Please also specify a file name for output.\n$usage\n";
  exit(1);
}

if(-e $out){
  print STDERR "ERROR: \'$out\' exists. Please check.\n$usage\n";
  exit(1);
}

my $currentDir = cwd();

unless(-w $currentDir){
  print STDERR "ERROR: No right to write in current directory. Please check.\n$usage\n";
  exit(1);
}

if($out){
  if($out =~ /\\/ || $out =~ /\//){
    print STDERR "ERROR: option '--out' only accept file name but not a path. Please reset the option \'--out\'.\n$usage\n"; 
    exit(1);
  }
  unless($out =~ /^\w+$/){
    print STDERR "ERROR: Output file name only allow letter, number, and underline. Please reset the option \'--out\'.\n$usage\n";
    exit(1);
  }
}

convert_W2L();

sub convert_W2L{
  open IN, $in;
  open OUT, '>>', $out;
  while(<IN>){
    chomp;
    my $tmp = $_;
    if(/\r$/){
      $tmp =~ s/\r$//;
    }
    print OUT "$tmp\n";
  }
  close IN;
  close OUT;
}




