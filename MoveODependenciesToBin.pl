#!/usr/bin/env perl
$in = shift @ARGV;
$out = shift @ARGV;
open(IN, "$in") or die "cannot open $in\n";
open(OUT, ">$out") or die "cannot open $out\n";
while(<IN>) {
		if ($_ =~ /(.*\.o.*)/) {
				print OUT "bin/$1\n";
		}
		else {
				print OUT $_;
		}
}
