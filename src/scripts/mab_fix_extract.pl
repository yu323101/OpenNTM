#!usr/bin/perl
use warnings;

while(<>){
chomp;
@a=split "\t",$_;
$a[4]=~s/%//;
if($a[4]>=75){
print "$_\n";
}
}
