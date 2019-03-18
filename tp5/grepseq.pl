use strict;

my $gbfile = shift;
my $selfile = shift;
my $outfile = shift;

open (SELFILE, $selfile) or die "input error\n";

my %idlist;

my $count = 0;

foreach my $line (<SELFILE>)	{
	chomp $line;
	$line =~ s/\s//g;
	$idlist{$line} = 1;
	$count++;
	#if ($line =~ /^((\w+)\_(\w+))\s/g)	{
	#	$idlist{$1} = 1;
	#	$count++;
	#}
	#else	{
	#	die "did not recognize $line\n";
	#}
}

print "$count items in selection\n";

open (GBFILE, $gbfile) or die "input error: gbfile\n";
open (OUTFILE, '>'.$outfile) or die "output error\n";

while (my $line = <GBFILE>)	{

	chomp $line;
	if ($line =~ /^ID\s+(\S+)\s/g)	{
		my $id = $1;
		if (exists $idlist{$id})	{
			$idlist{$id} = 0;
			while ($line !~ /^\/\/$/g)	{
				print OUTFILE "$line\n";
				$line = <GBFILE>;
				chomp $line;
			}
			print OUTFILE "\/\/\n";
		}
		else	{
			while ($line !~ /^\/\/$/g)	{
				$line = <GBFILE>;
				chomp $line;
			}
		}
	}
	else	{
		die "did not recognize $line\n";
	}
}

foreach my $id (keys %idlist)	{
	if ($idlist{$id} == 1)	{
		print "did not find $id\n";
	}
}

