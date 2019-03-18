use strict;

my $infile = shift;
my $outfile = shift;

open (INFILE, $infile) or die "input error\n";
open (OUTFILE, '>'.$outfile) or die "output error\n";

my %name2seq;

while (my $line = <INFILE>)	{

	chomp $line;

	if ($line =~ /^(\S+)\s+(\S+)$/g)	{
		my $name = $1;
		my $seq = $2;
		$name2seq{$name} = $name2seq{$name} . $seq;
	}
}

foreach my $name (keys %name2seq)	{
	print OUTFILE ">$name\n$name2seq{$name}\n";
}

