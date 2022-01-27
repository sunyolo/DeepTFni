use strict;

my $motifDir = "/Your_path/data_resource/human_HOCOMO/";
opendir (DIR,"$motifDir") || die "Error";
my @file = grep {/meme/}readdir(DIR);
closedir (DIR);
my %motif;
my $n;
foreach my $m(@file){
	my $name = (split /_/,$m)[0];
	$motif{$name} = 1;
	$n++;
}
print "total TF motif: $n\n";

open (f1,"./gencode.v19.ProteinCoding_gene_promoter.txt") || die "Error";
open (o1,">./gencode.v19.TF.promoter.temp") || die "Error";
my $n;
while (<f1>)
{
	$_=~s/\s+$//;
	my @a = split /\t/,$_;
	next if $motif{$a[-1]} eq "";
	$n++;
	print o1 "$_\n";
}
close f1;
close o1;
print "TF motif in gencode19: $n\n";
system ("sort-bed gencode.v19.TF.promoter.temp > gencode.v19.TF.promoter.txt");
system ("rm gencode.v19.TF.promoter.temp");

my $input = $ARGV[0];
my $dir = "./TFBS_region/";
my $outtemp = "./bedops_temp/";
mkdir $outtemp unless -e $outtemp;
system ("bedops -e -100% $dir$input.txt gencode.v19.TF.promoter.txt > $outtemp$input.txt");# unless -e "$outtemp$input.txt"
my $n = count("$outtemp$input.txt");
print "TFBS in TF promoter: $n\n";

my %hash;
my %val;
open (f1,"gencode.v19.TF.promoter.txt") || die "Error";
my $tempCount;
while (<f1>)
{
	$_=~s/\s+$//;
	print "read $tempCount\n" if $tempCount%100 == 0;
	$tempCount++;
	my @a = split /\t/,$_;
	my ($chr,$start,$end,$TF) = ($a[0],$a[1],$a[2],$a[4]);
	open (f2,"$outtemp$input.txt") || die "Error";
	my $sum;
	while (<f2>)
	{
		$_=~s/\s+$//;
		my @a = split /\t/,$_;
		next if $a[0] ne $chr;
		if ($a[1] >= $start && $a[2] <= $end){
			my $p = join "\t",$TF,$a[3];
			my $q = join "\t",$a[3],$TF;
			$hash{$p} = 1;$hash{$q} = 1;
			$sum++;
		}
	}
	close f2;
	$val{$TF} = 1 if $sum > 0;
}
close f1;
my $n = scalar(keys %val);
print "val TF number: $n\n";
my $outdir = "./Adjacency_matrix/";
mkdir $outdir unless -e $outdir;
open (o1,">./$outdir/$input.txt") || die "Error";
open (o2,">./$outdir/$input.TF.list.txt") || die "Error";
print o1 "TF";
foreach my $m(keys %val){print o1 "\t$m";}
print o1 "\n";
my $sum;
foreach my $m(keys %val){
	print o1 "$m";
	print o2 "$m\n";
	foreach my $n(keys %val){
		my $p = join "\t",$m,$n;
		$hash{$p} = 0 if $hash{$p} eq "";
		$sum+=$hash{$p};
		print o1 "\t$hash{$p}";
	}
	print o1 "\n";
}
close o1;
my $density = sprintf("%.2f",$sum/scalar(keys %val)/scalar(keys %val));
print "$input density:$sum\t$density\n";


my $out_dir = "./$input/";
mkdir $out_dir unless -e $out_dir;
my $input_file = "./Input_data/$input.csv";
print ("input_file = $input_file\n");
system ("python 6_dense_matrix_2_h5_for_MAESTRO_input.py -sample $input -input_file $input_file -output_dir $out_dir");

sub count{
	my @a = @_;
	my $n;
	open (f1,"$a[0]") || die "Error $a[0]";
	while (<f1>){$n++;}
	close f1;
	return $n;
}
