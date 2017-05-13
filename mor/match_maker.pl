#!/usr/bin/perl
#
# "Matchmaker, Matchmaker, Make me a match, find me a find, catch me a catch."  Anways, this scripts goes through two tables and tries to make
# a lot of connections between the two.  It then uses some mysql conditions to say that something is a match or not and then marks it in
# a third table which has an option for human interaction to break the match or not. (It works better if you sing the song while it's running.)

use strict;

chdir '/misc/mor/';
use vars qw($log %settings $maindb);

use mor::DatabaseIO;
use mor::Log;
use mor::Settings;
use mor::Utilities;
use Bio::AlignIO;
use Bio::SeqIO;

*log = \$mor::Log::log;
$log = new mor::Log("logs/.log");
*maindb = \$mor::DatabaseIO::maindb;
*settings = \%mor::Settings::settings;

# Load the settings from the settings file
mor::Settings::LoadSettings("settings/settings.dat");

my $its = $settings{'ITStable'};
my $lsu = $settings{'LSUtable'};
my $clusterTable = "clusters";
my $match = "matches";

mor::Utilities::lockProcessOrDie("match_maker");

#Connect the database
mor::DatabaseIO::LoadMainDB();

# get an array of all the fields from lsu
my $lsuQuery = "SELECT accno, clade, species, aftol_id, db_xref, isolate, culture_collection, specimen_voucher, strain FROM ". $lsu ." WHERE accepted = 1 ORDER BY date";

my $dbArray = $maindb->query($lsuQuery);
my @entries;
while(my @row = $dbArray->fetchrow_array()){
    push(@entries, \@row);
}

my $count = $maindb->getValue("SELECT COUNT(accno) FROM $lsu where accepted = 1");
my @fields = ("species", "aftol_id", "db_xref", "isolate", "culture_collection", "specimen_voucher", "strain");
# for each otu in lsu
for(my $x = 0; $x < $count; $x++) { # This would be a perfect place to multithread this script because it's going to take a LONG time
    if($x % 25 == 0) {
	$log->print("On LSU sequence number $x out of $count", 1);
    }
    my $lsuRowRef = @entries[$x];
    my @lsuRow = @{$lsuRowRef};
    # check each field after accno (that's why we start at 1) 
    for( my $i = 2; $i < scalar(@lsuRow); $i++){
	# if the field is not null
	if($lsuRow[$i] ne "NULL"){
	    # get the its info if its present
	    my $query = "SELECT accno, $fields[$i-2], cluster FROM $its where $fields[$i-2] != \"NULL\" ORDER BY date";
	    my $itsArray = $maindb->query($query);
	    my @itsMatches;

	    while(my @row = $itsArray->fetchrow_array()){
		push(@itsMatches, \@row);
	    }

	    # then for each its sequence
	    my $itsCount = scalar(@itsMatches);
	    for(my $y = 0; $y < $itsCount; $y++){
		my $itsRowRef = @itsMatches[$y];
		my @itsRow= @{$itsRowRef};
		my $value = 0;
		if($itsRow[1] ne "NULL"){
		    #clean up the its field
		    my $itsTrimmed = $itsRow[1];
		    $itsTrimmed =~ s/\s//g;
		    $itsTrimmed =~ s/\W//g;
		    #clean up the lsu field
		    my $lsuTrimmed = $lsuRow[$i];
		    $lsuTrimmed =~s/\s//g;
		    $lsuTrimmed =~s/\W//g;
		    
		    # check to see if it has a loose match
		    if($lsuTrimmed=~ m/^$itsTrimmed$/i) {
			# note if it has a loose match
			$value = 1;
			# if it strictly matches note that instead
			if($itsRow[1] =~ m/^$lsuRow[$i]$/i){
			    $value = 2;
			}
			# figure out if there's already a relation between the two accnos
			my $matchNumber = $maindb->getValue("SELECT match_num FROM $match WHERE lsu_accno = \"$lsuRow[0]\" AND its_accno = \"$itsRow[0]\"");
			# if there's not
			if (!$matchNumber) { 
			    # we create a new match
			    $matchNumber = $maindb->getValue("SELECT max(match_num) + 1 FROM $match");
			}
			# Then insert the appropriate data into table
			my $insert = "INSERT INTO $match (match_num, match_name, lsu_accno, its_accno, cluster_number, ".$fields[$i-2].", clade) ";
			my $values = "VALUES ($matchNumber, \"".$lsuRow[2]."\", \"".$lsuRow[0]."\", \"".$itsRow[0]."\", \"".$itsRow[2]."\", \'$value\', \'".$lsuRow[1]."\') ";
			my $update = "ON DUPLICATE KEY UPDATE ".$fields[$i-2]." = $value WHERE lsu_accno = \"$lsuRow[0]\" AND its_accno = \"$itsRow[0]\"";
			$query = $insert . $values . $update;
			$maindb->do($query);   
		    }
		}
	    }
	}
    }
}

# using the matches mark the solid links between the tables
# this is working with only the matches table
$maindb->do("UPDATE $match SET machine_matched = 1 WHERE species > 0 AND db_xref > 0 ");
$maindb->do("UPDATE $match SET machine_matched = 2 WHERE (aftol_id > 0 OR isolate > 0 OR strain > 0 OR culture_collection > 0 OR specimen_voucher > 0) AND (species > 0 OR db_xref > 0)");

#This is for the clusters table
$log->print("Using the match information to update the tables displayed to the website");
my $maxCluster = $maindb->getValue("SELECT MAX(cluster) FROM $its");
for(my $i = 1; $i <= $maxCluster; $i++){
    my $resultRef = $maindb->getArray("SELECT lsu_accno, clade, machine_matched FROM $match WHERE cluster_number = $i AND machine_matched > 0 AND human_broken != 1");
    my @results = @{$resultRef};
    foreach my $ref (@results) {
	my @row = @{$ref};
	$maindb->do("");
    }
}

mor::Utilities::clearProcessLocks("match_maker");
