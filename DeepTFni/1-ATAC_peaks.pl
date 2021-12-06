use strict;


my $dir = "./Input_data/";
my $outdir = "./ATAC_peaks/";
mkdir $outdir unless -e $outdir;

my $c = $ARGV[0];
my $reference = $ARGV[1];
{
	print "$c\n";
	open (f1,"$dir$c.csv") || die "Error";
	readline(f1);
	open (o1,">$outdir$c.clean.ATAC_peak.temp") || die "Error";
	open (o2,">$outdir$c.clean.ATAC_peak.matrix.temp") || die "Error";
	while (<f1>)
	{
		$_=~s/\s+$//;
		my @a = split /,/,$_;
		my ($n,$total);
		my ($chr, $start, $end) = split /_|\t|:|-/,$a[0],3;
		for (my $i = 1;$i<scalar(@a);$i++)
		{
			$n++ if $a[$i] > 0;
			$total++;
		}
		my $per = $n/$total;
		print o1 "$chr\t$start\t$end\t$per\n" if $per > 0.1;
		if ($per > 0.1){
			print o2 "$chr\t$start\t$end";
			for (my $i = 1;$i<scalar(@a);$i++){print o2 "\t$a[$i]";}
			print o2 "\n";
		}
	}
	close f1;
	close o1;
	close o2;
	system ("sort-bed $outdir$c.clean.ATAC_peak.temp > $outdir$c.clean.ATAC_peak.bed");
	system ("sort-bed $outdir$c.clean.ATAC_peak.matrix.temp > $outdir$c.clean.ATAC_peak.matrix.txt");
	system ("rm $outdir$c.clean.ATAC_peak.temp");
	system ("rm $outdir$c.clean.ATAC_peak.matrix.temp");
}

system ("perl 2-SequenceGet.pl $outdir$c.clean.ATAC_peak.bed $reference");

