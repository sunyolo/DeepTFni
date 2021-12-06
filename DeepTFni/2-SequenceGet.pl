use strict;

use lib qw(BioPerl-1.7.4); 
$| = 1;
use Getopt::Long;
use File::Spec;
use File::Basename;
use Bio::DB::Fasta;

my $seqs;

my $input = $ARGV[0];
my $reference = $ARGV[1];
print STDERR "reading genome file $reference\n";
my $GENOMEFILE;
if ($reference eq "hg19") {
	$GENOMEFILE = "/Your_path/hg19/index";
}
else {system ("unsupported reference version\n");}

my $seqDb =  Bio::DB::Fasta->new($GENOMEFILE);
print "seqDb info : $seqDb\n";
my @ids      = $seqDb->get_all_primary_ids;
print "id info : @ids\n";

my $outdir = "./Fasta/";
mkdir $outdir unless -e $outdir;

my $name = (split /\/|\.bed/,$input)[2];
print "name : $name\n";
open (f1,"$input") || die "Error $input";
open (o1,">$outdir$name.fasta") || die "Error";
while (<f1>)
{
	$_=~s/\s+$//;
	my @a = split /\t/,$_;
	my ($chr,$start,$end,@rest)=split /\t/,$_;
	next if $chr eq "chrY";
	print o1 ">$a[0]\-$a[1]\-$a[2]\-$a[4]\n";
	my $Seq = uc($seqDb->get_Seq_by_id($chr)->subseq($start=>$end));
	print o1 "$Seq\n";
}
close f1;
close o1;

system ("perl 3-fimo.pl $name");
