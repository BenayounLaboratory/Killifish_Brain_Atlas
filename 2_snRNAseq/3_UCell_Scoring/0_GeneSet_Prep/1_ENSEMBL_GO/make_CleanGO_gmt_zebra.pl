#! /usr/bin/perl

use warnings;
use strict;
use List::MoreUtils qw(uniq);

# a script to extract a gene list from big files

unless (scalar @ARGV == 1) {
	die "\n make_CleanGO_gmt_zebra.pl <ENSEMBL GO Ann file> .\n";

}


# database file
my $database = shift @ARGV;

unless (open(DATA, $database)) {
	die "cannot open $database file. \n";
}

# read in list of genes in array
my %database_BP  = ();
my %database_CC  = ();
my %database_MF  = ();
my %database_ALL = ();

while (my $line = <DATA>) {
	next if ($. == 1); # skip header
	
	my @tmp = get_line_data ($line);
	
	next if (scalar @tmp < 6   ); # skip lines without GO terms
	next if ($tmp[4] eq "NAS");  # Do Not use GO with just author statements
	next if ($tmp[4] eq "TAS");  # Do Not use GO with just author statements
	next if ($tmp[4] eq "ND");   # Do Not use GO without evidence

	
	my $GOterm = join("_",($tmp[2],$tmp[3]));
	#print $GOterm."\n";
	
	push (@{$database_ALL{$GOterm}}, $tmp[1]);
	
	if ($tmp[5] eq "molecular_function") {
		push (@{$database_MF{$GOterm}}, $tmp[1]);
		
	} elsif ($tmp[5] eq "cellular_component") {
		push (@{$database_CC{$GOterm}}, $tmp[1]);
	
	} elsif ($tmp[5] eq "biological_process") {
  		push (@{$database_BP{$GOterm}}, $tmp[1]);
    
	} else {
		print "GO domain unclear\n";
	}
}

close DATA;

#print keys %database;
open(GMT, '>', "2024-03-20_Zebrafish_Ens111_GO_ALL.gmt");
foreach my $GOterm (keys %database_ALL) {
	print GMT $GOterm."\tEns111\t".join("\t", uniq @{$database_ALL{$GOterm}} )."\n";
}
close GMT;

open(GMT, '>', "2024-03-20_Zebrafish_Ens111_GO_BP.gmt");
foreach my $GOterm (keys %database_BP) {
	print GMT $GOterm."\tEns111\t".join("\t", uniq @{$database_BP{$GOterm}} )."\n";
}
close GMT;

open(GMT, '>', "2024-03-20_Zebrafish_Ens111_GO_MF.gmt");
foreach my $GOterm (keys %database_MF) {
	print GMT $GOterm."\tEns111\t".join("\t", uniq @{$database_MF{$GOterm}} )."\n";
}
close GMT;

open(GMT, '>', "2024-03-20_Zebrafish_Ens111_GO_CC.gmt");
foreach my $GOterm (keys %database_CC) {
	print GMT $GOterm."\tEns111\t".join("\t", uniq @{$database_CC{$GOterm}} )."\n";
}
close GMT;

exit;

###########################################################
# SUBROUTINES
###########################################################

###########################################################
# a subroutine that separates fields from a data line and
# returns them in an array

sub get_line_data {

    my $line = $_[0];
    
    chomp $line;  

    my @linedata = split(/\t/, $line);
        
    return @linedata;
}