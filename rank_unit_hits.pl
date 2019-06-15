#! perl
use warnings;


my $in1 = shift;
my $in2 = shift;
my $out = shift;
my ($id,$code);
open IN1, $in1;
our %code_hash;
while(<IN1>){
  chomp;
  if(/^>MER(\d+)\s.+\[(.+?\..+?)\]/){
    $id = $1;
    $code = $2;
    $code_hash{$id} = $code;
  }
}
close IN1;


open IN2, $in2;
open OUT, '>>', $out;
our %query_hash;
our @query_array;
my $count = 0;
while(<IN2>){
  chomp;
  if(/^Query=\s(.+)$/){
    %query_hash = ();
    $code_id = '';
    print OUT "Query=$query_id\n" if $query_id;
    if($count == 1){
      foreach(@query_array){
        $code_id = $_;
        print OUT "$code_id\t@{$code_id}\n";
        @{$code_id} = ();
      }
    }
    $query_id = $1;
    $count = 0;
    @query_array = ();
  }
  if(/^>\sMER(\d+)\s/){
    $count = 1;
    $id1 = $1;
    $query_hash{"$code_hash{$id1}"} += 1;
    if($query_hash{"$code_hash{$id1}"} == 1){
      push @query_array, $code_hash{$id1};
      push @{$code_hash{$id1}}, $id1;
    }
    if($query_hash{"$code_hash{$id1}"} > 1){
      push @{$code_hash{$id1}}, $id1;
    }
  }
}
close IN2;

print OUT "Query=end\n";



