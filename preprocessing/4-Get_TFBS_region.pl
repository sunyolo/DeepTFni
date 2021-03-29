use strict;

my $outdir = "./TFBS_region/";
mkdir $outdir unless -e $outdir;

my $dataset = $ARGV[0];
{
	print "dataset : $dataset\n";
	my $dir = "./fimo-res/$dataset.clean.ATAC_peak_1e4/";
	opendir (DIR,"$dir") || die "Error";
	my @file = grep {/txt/}readdir (DIR);
	closedir (DIR);

	open (o1,">$outdir$dataset.temp") || die "Error";
	my $fn;
	foreach my $m(@file)
	{
		my $name = (split /\.txt/,$m)[0];
		open (f1,"$dir$m") || die "Error";
		open (o2,">$dir$name.temp") || die "Error";
		while (<f1>)
		{
			$_=~s/\s+$//;
			last if($_ eq '');
			my @a = split /\t/,$_;
			next if $a[0] =~ /motif_id/;
			if ($a[7] > 1e-6){next;}
			my $tf = (split /_/,$a[0])[0];
			my ($chr, $start, $end) = split /-/,$a[2],3;
			($start, $end) = ($start+$a[3]-1, $start+$a[4] - 1);
			print o1 "$chr\t$start\t$end\t$tf\n";
			print o2 "$chr\t$start\t$end\t$tf\n";
		}
		close f1;
		close o2;
		system ("sort-bed $dir$name.temp > $dir$name.bed");
		system ("rm $dir$name.temp");
	}
	close o1;

	system ("sort-bed $outdir$dataset.temp > $outdir$dataset.txt");
	system ("rm $outdir$dataset.temp");
	print "run adjance matrix construction\n";
	system ("perl 5-GetMatrix.pl $dataset");
}
