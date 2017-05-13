#!/usr/bin/perl
#
# This is the main script that calls the mor modules responsible for downloading
# new sequences, screening and aligning. No further analysis is done at that
# point, since there need not be an explicit sequence whereby testmor downloads new
# sequences and then immediately processes them.
use strict;

chdir '/misc/mor/';
use vars qw($log %settings $maindb);

use mor::Acquisition;
use mor::DatabaseIO;
use mor::Log;
use mor::Screening;
use mor::Settings;
use mor::Utilities;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Seq::RichSeq;


*log = \$mor::Log::log;
$log = new mor::Log("logs/.log");
*maindb = \$mor::DatabaseIO::maindb;
*settings = \%mor::Settings::settings;

# Load the settings from the settings file
mor::Settings::LoadSettings("settings/settings.dat");

my $lsu = $settings{'LSUtable'};

mor::Utilities::lockProcessOrDie("download");

#Connect the database
mor::DatabaseIO::LoadMainDB();

# Query GenBank up to three times. Program aborts if nothing found.
my @seqs = mor::Acquisition::tryQueryGB("","","",$lsu,"i");

# should check to see if genbank graced us with any sequences instead of throwing an error
if (scalar(@seqs) ne 0 ) {
    
# Process the sequences downloaded from GenBank (accept or reject) and also stores them in the passed table
    $log->print("Sorting and depositing the sequences into [". $lsu ."] MySQL table.",1);
    mor::Screening::sortSequences(@seqs,$lsu,$settings{'primerA'},$settings{'primerB'},$settings{'HMMfile'});
    
    my $FASTAin = mor::Utilities::getUniqueFile($settings{'alignIn'});
    my $FASTAout = mor::Utilities::getUniqueFile($settings{'alignOut'});
    
# Gets arrays of the loosely accepted sequences and accession numbers from the table

    $log->print("Getting the accepted sequences from ". $lsu,1);

    my $accnoRef = $maindb->getArray("SELECT accno FROM ". $lsu ." WHERE accepted = 1 AND HRI != 0 ORDER BY date ASC");
    my @accnos = @{$accnoRef};
    my $seqRef = $maindb->getArray("SELECT accepted_seq FROM ". $lsu ." WHERE accepted = 1 AND HRI != 0 ORDER BY date ASC");
    my @sequences = @{$seqRef};
    my $nameRef = $maindb->getArray("SELECT species FROM ". $lsu ." WHERE accepted = 1 AND HRI != 0 ORDER BY date ASC");
    my @names = @{$nameRef};
    my @looseAccepted;

    for(my $i = 0; $i < scalar(@accnos); $i++){
	if($accnos[$i] eq ""){ 
	    $log->print("Found the accno at $i");
	    sleep 60;
	}
	push(@looseAccepted, (Bio::Seq::RichSeq->new(-seq => $sequences[$i], -alphabet => "dna", -accession_number => $accnos[$i], -display_id => $accnos[$i], -primary_id=>$names[$i])));
    }

# Check for redundant species
    $log->print("Checking for redundant sequences in $lsu",1);
    for ( my $i = 0; $i < scalar(@looseAccepted); $i++) {
	$log->print("Creating group $i of a potential ". scalar(@looseAccepted) ."", 1);
	my @group;
	$group[0] = shift(@looseAccepted);
	for (my $j = 0; $j < scalar(@looseAccepted); $j++){
	    if($group[0]->primary_id eq $looseAccepted[$j]->primary_id ){
		push(@group, splice(@looseAccepted,$j,1));
	    }
	}
	if(scalar(@group > 1)){
	    $log->print("Screening group ". $group[0]->primary_id ." of size ". scalar(@group) ."", 1);
	    mor::Screening::sequenceSimilarity(@group, $lsu);
	}
    }

# Get all of the non-redundant sequences

    my $alnObj = $maindb->getAlignObject($lsu, "accepted_seq", "accepted = 1");

# Writes the array of sequence objects to a file
    $log->print("Writing the sequence objects for processing.",1);
    my $out = Bio::AlignIO->new(-file => ">$FASTAin", -format => "fasta");
    $out->write_aln($alnObj);
    
# Align the sequences using MAFFT
# This could maybe use some better control and wrapping but it works for right now
    $log->print("Aligning the sequences.",1);
    my $cmd = "mafft --auto ". $FASTAin ." > " . $FASTAout;
    system($cmd);
    
# get the alignment and store it mysql table
    my $in  = Bio::SeqIO->new(-file   => $FASTAout , -format => 'fasta');
    
    while ( my $obj = $in->next_seq() ) {
	my $name = $obj->display_id;
	$name =~ s/\/\d*\-\d*$//;
	my $query = "UPDATE $lsu SET align_seq = \"". $obj->seq ."\" WHERE accno=\"". $name ."\"";
	$maindb->do("$query");
    }
    
    $log->print("Sequence acquisition/alignment complete completed.",1);
} else {

    $log->print("There were no new records from genbank to add.");

}

mor::Utilities::clearProcessLocks("download");
