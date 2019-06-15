#! /usr/bin/perl
use warnings;

my $in1 = shift;
my $in2 = shift;
my $out = shift;

open IN, $in1;

while(<IN>){
  chomp;
  if(/^Query=(.+)$/){
    $query_id = $1;
  }
  if(/^(\d+)\s(.+)$/){
    my $mer_id = $1;
    my $domain = $2;
    $$query_id{$mer_id} = $domain;
  }
}
close IN;

$query_id = '';

open IN, $in2;
open OUT, '>>', $out;

while(<IN>){
  chomp;
  my $tmp = $_;
  if(/^Query=(.+)$/){
    $query_id = $1;
    print OUT "Query=$query_id\n";
  }else{
    if(/^(.+?)\s/){
      $code = $1;
    }
    $tmp =~ s/$code\s//;
    my @mer_id_array = split /\s/, $tmp;
    print OUT "code=$code\n";
    foreach(@mer_id_array){
      my $mer_id = $_;
      if($$query_id{$mer_id}){
        print OUT "$mer_id $$query_id{$mer_id}\n";
        #$switch_next++;
      }else{
        print "ERROR $query_id $mer_id\n";
      }
    }
  }
}
close IN;
close OUT;


