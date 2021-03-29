#!/usr/bin/perl
#use warnings; 
#use strict;
use File::Basename;
use File::Copy;
use Parallel::ForkManager;

$dir  = @ARGV[0];
$reference  = @ARGV[1];	
@suffixlist=qw(.csv);

sub run_DeepTFni{				
    my $matrix_file	= $_[0];   
	my @a = split/\.csv/,$matrix_file;
	my $sample_name = $a[0];
    print STDOUT "Running: DeepTFni.v1 --file $matrix_file\n"; 
	
	open (f1,"./run_template_csv.sh") || die "Error";
	open (o1,">./run_$sample_name.sh") || die "Error";
	while (<f1>)
		{
			my $s1 = substr($_,0,6);
			my $s2 = substr($_,0,9);
			if ($s1 eq "sample") {
				print o1 ('sample_type="');
				print o1 ("$sample_name");
				print o1 ('"');
				print o1 ("\n");
			}
			elsif ($s2 eq "reference") {
				print o1 ('reference="');
				print o1 ("$reference");
				print o1 ('"');
				print o1 ("\n");
			}
			else{
				print o1 $_;
			}	
		}

	close f1;
	close o1;	
	system("sh ./run_$sample_name.sh") == 0 or die "Error in running run_$sample_name.sh\n";	
}

opendir (DIR,"$dir") || die "Error : cannot open input directory";
@file_list = grep{/csv/}readdir(DIR);
closedir (DIR);

$nprocs = 5; 
$MAX_PROCESSES = $nprocs;
$pm = new Parallel::ForkManager($MAX_PROCESSES);
print STDOUT "Start: An $nprocs parallel manager!\n==================\n";

DATA_LOOP:
foreach my $file (@file_list) {
      
	my $file_name = basename($file,@suffixlist);
	my @a = split/\./,$file;
	my $type = $a[0];
	print $type;
	my $pid = $pm->start and next DATA_LOOP;
	#print STDOUT "pid	:	$pid\n";
    	run_DeepTFni($file);
    	$pm->finish; 
}
$pm->wait_all_children;
undef @file_list;


