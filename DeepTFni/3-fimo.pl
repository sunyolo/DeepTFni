use strict;
use File::Copy;

my $fimo = "/Your_path/fimo";
my $dir = "/Your_path/data_resource/human_HOCOMO/";

opendir (DIR,"$dir") || die "Error";
my @file = grep{/meme/}readdir (DIR);
closedir (DIR);

my $fimoout = "./fimo-res/";
mkdir $fimoout unless -e $fimoout;

my $cell = $ARGV[0];
{
	print "$cell\n";
	my $input = "./Fasta/$cell";
	my $outdir = "./fimo-res/$cell\_1e4";
	mkdir $outdir unless -e $outdir;
	foreach my $m(@file)
	{
		system "$fimo --oc $outdir --no-qvalue $dir$m $input.fasta";
		#system "$fimo --oc $outdir --thresh 1e-5 --no-qvalue $dir$m $input.fasta";
		my @a = split /_HUMAN.H10MO./,$m;
		my $b = (split /\.meme/,$a[1])[0];
		my $name = join ".",$a[0],$b;
		print "\t$name\n";
		system ("mv $outdir/fimo.tsv $outdir/$name.txt");		
	}
}
my $c = (split /\./,$cell)[0];

system ("perl 4-Get_TFBS_region.pl $c");
