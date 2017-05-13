#!/usr/bin/perl
#
# tansley.pl is a rewrite of the previous tansley code with mor code; this was done to minimize the code base to maintain.  It works pretty much the same as mor_download.pl, but instead 
# of reducing the number of sequences based on HRI it clusters them together for further analysis.  Make sure when you're running this, it has your oldest sequence otherwise if you would
# like to include it I would recommend that you truncate the table and start again.

use strict;

chdir '/misc/mor/';
use vars qw($log %settings $maindb);

use mor::Acquisition;
use mor::DatabaseIO;
use mor::Log;
use mor::Screening;
use mor::Settings;
use mor::Utilities;
use Bio::Align::AlignI;
use Bio::AlignIO;
use Bio::SeqIO;

*log = \$mor::Log::log;
$log = new mor::Log("logs/.log");
*maindb = \$mor::DatabaseIO::maindb;
*settings = \%mor::Settings::settings;

# Load the settings from the settings file
mor::Settings::LoadSettings("settings/settings.dat");

my $its = $settings{'ITStable'};
my $clusterTable = "clusters";
my $minimum_similarity_value = 90;
mor::Utilities::lockProcessOrDie("tansley");

#Connect the database
mor::DatabaseIO::LoadMainDB();

my $queryPhrase = "Agaricomycetes[ORGN] AND 500[SLEN]:1000[SLEN] AND (\"ITS1\" OR \"ITS 1\" OR \"internal transcribed spacer 1\") AND (\"ITS2\" OR \"ITS 2\" OR \"internal transcribed spacer 2\") NOT mitochondrial";

# Query GenBank up to three times. Program aborts if nothing found.
my @seqs = mor::Acquisition::tryQueryGB("","","$queryPhrase",$its,"i");

# should check to see if genbank graced us with any sequences instead of throwing an error
if (scalar(@seqs) ne 0 ) {
    
# Process the sequences downloaded from GenBank (accept or reject) and also stores them in the passed table
    $log->print("Sorting and depositing the sequences into [". $its ."] MySQL table.",1);
    mor::Screening::sortSequences(@seqs,$its,"","","");
    
# Gets arrays of the loosely accepted sequences and accession numbers from the table
    $log->print("Getting the accepted sequences from [". $its ."]",1);


    my $cladeFile = "clade.txt";
    open(FILE, $cladeFile);
    my @cladeList = <FILE>;
    close(FILE);
    
    foreach my $clade (@cladeList) {
	chomp($clade);
# if that's done it could only check the keys against the new sequences
# pull the unclustered sequences that match the clade
	my $uncObj = $maindb->getAlignObject("its", "original_seq", "accepted = 1 AND cluster_key = 1 AND cluster = 0 AND taxonomy LIKE \"\%\:$clade\:\%\"");
	my $keyObj = $maindb->getAlignObject("its", "original_seq", "accepted = 1 AND cluster_key = 1 AND cluster != 0 AND taxonomy LIKE \"\%\:$clade\:\%\"");
	my $numSeqs = $maindb->getValue("select count(accno) from its where accepted = 1 and cluster_key = 1 and taxonomy like \"\%\:$clade\:\%\"");
	my @unclustered = ();
	$log->print("Clustering the $clade portion of the [$its] table with $numSeqs",1);
	my @keySeqs = ();
	if($numSeqs ne 0){
	    @unclustered = $uncObj->each_seq;
	    @keySeqs = $keyObj->each_seq;
	} 
        # Cluster the its data
	while(@unclustered) {
	    
	    # This creates clusters that are keyed to one sequence and removes all the associated sequences so there are no overlaps between clusters
	    # it might make more sense to have overlapping clusters for supertree analysis
	    my $seq;
	    if(scalar(@keySeqs) == 0){
		$seq = shift(@unclustered);
	    } else {
		$seq = shift(@keySeqs);
	    }
	    
	    my ($uncRef, $clusteredSeqsRef, $percentsRef) = &cluster($minimum_similarity_value, $seq, \@unclustered);

	    @unclustered = @{$uncRef};
	    my @clusteredSeqs = @{$clusteredSeqsRef};
	    my @percents = @{$percentsRef};
	    
	    # we now have a "new" cluster, we still have to include the older version of that cluster if there is one
	    my $cluster = "";
	    my $num = 0;
	    my $accnoList = "";
	    
	    # pull the first sequence and check to see if it already has a cluster
	    my $clusterKey = shift(@clusteredSeqs);
	    if ($maindb->getValue("SELECT cluster FROM $its WHERE accno = \'".$clusterKey->display_id()."\'")) {
		# pull its cluster number
		$cluster = $maindb->getValue("SELECT cluster FROM $its WHERE accno = \'".$clusterKey->display_id."\'");
		$log->print("Cluster is present, getting the associated information for cluster number: $cluster", 1);
		# pull all the associated accnos
		$accnoList = $maindb->getValue("SELECT similar_accnos FROM $its WHERE accno = \'".$clusterKey->display_id()."\'");
	    } else {
		$cluster = $maindb->getValue("SELECT cluster + 1 from $its ORDER BY cluster DESC limit 0,1");
		$log->print("Creating cluster number: ". $cluster ." with sequence ". $clusterKey->display_id);
		$accnoList = $clusterKey->display_id() ."; ";
	    }
	    
	    # rework this so the table doesn't contain a lot of redundent information, the cluster info should be stored in the "clusters" table
	    # Do a little bit of prep work to create information for the cluster
	    foreach my $element (@clusteredSeqs) {
		my $display = $element->display_id();
		#make sure we have a catch on the old data that could be in there
		if($accnoList !~ /$display/ ){
		    $accnoList = $accnoList . $element->display_id()."; ";
		}
	    }
	    # remove the trailing "; "
	    $accnoList =~ s/\W*$//;
	    
	    my @clusteredAccnos = split('; ', $accnoList);
	    $num = scalar(@clusteredAccnos);
	    
	    # Store the information
	    # because the first one needs to still be flagged as important to clustering
	    $log->print("Storing cluster: $cluster of size $num into the database", 1);
	    my $key = shift(@clusteredAccnos);
	    $maindb->do("UPDATE $its SET num_associated = \'$num\', cluster = \'$cluster\', similar_accnos = \'".$accnoList."\', percent = \'".$percents[0]."\', cluster_key=1 WHERE accno = \'".$key."\'") or die "Updating the clusterKey failed.\n";
	    # the rest can be now marked as already clustered
	    my $place = 1; 
	    foreach my $element (@clusteredAccnos) {
		my $updateNonKeys = "UPDATE $its SET num_associated=\'$num\', cluster=\'$cluster\', similar_accnos=\'".$accnoList."\', percent = \'$percents[$place]\', cluster_key = 0 WHERE accno = \'".$element."\'";
		$maindb->do($updateNonKeys) or die "Updating the rest of the cluster failed.\n";
		$place++;
	    }
	}
    }
    
#Now to check for mistakes and finish the unclustered sequences
    $log->print("Starting error checking in clusters", 1);
    my $alnObj = $maindb->getAlignObject($its, "original_seq", "accepted = 1 AND cluster_key = 1");
    my @accepted = $alnObj->each_seq;
    my $count = scalar(@accepted);
    #loop through doing pairwise alignments
    while(@accepted) {
	
	my $seq = shift(@accepted);
	my ($acceptedRef, $clusteredSeqsRef, $percentsRef) = &cluster($minimum_similarity_value, $seq, \@accepted);
	@accepted = @{$acceptedRef};
	my @clusteredSeqs = @{$clusteredSeqsRef};
	my @percents = @{$percentsRef};
	
	if(scalar(@clusteredSeqs > 1)){
	    my $goodSeq = shift(@clusteredSeqs);
	    foreach my $badSeq (@clusteredSeqs){
		my $badCluster = $maindb->getValue("SELECT cluster FROM ".$its." WHERE accno = \"".$badSeq->display_id()."\"");
		if($badCluster ne 0){
		    my $resetQuery = "UPDATE ".$its." SET cluster_key = \'1\', num_associated = \'0\', similar_accnos = \"\" WHERE cluster = ".$badCluster;
		    $maindb->do($resetQuery) or die "Couldn't remove cluster $badCluster.";
		    #reset the number for everything about it down one
		    $resetQuery = "UPDATE ".$its." SET cluster = cluster-1 WHERE cluster > ".$badCluster;
		    $maindb->do($resetQuery) or die "Couldn't decrement the clusters.";
		}
	    }
	    my $date = "SELECT date FROM ".$its." WHERE accno = \"".$goodSeq->display_id()."\"";
	    # This won't work.... try limiting to the last/first N sequences ?
	    my $progress = $count - scalar(@accepted);
	    $alnObj = $maindb->getAlignObject($its, "original_seq", "accepted = 1 AND cluster_key = 1");
	    @accepted = $alnObj->each_seq;
	    splice(@accepted, 0, $progress);
	}	
	
	# we now have a "new" cluster, we still have to include the older version of that cluster if there is one
	my $cluster = "";
	my $num = 0;
	my $accnoList = "";
	
	# pull the first sequence and check to see if it already has a cluster
	my $clusterKey = shift(@clusteredSeqs);
	if ($maindb->getValue("SELECT cluster FROM $its WHERE accno = \'".$clusterKey->display_id()."\'")) {
	    # pull its cluster number
	    $cluster = $maindb->getValue("SELECT cluster FROM $its WHERE accno = \'".$clusterKey->display_id."\'");
	    $log->print("Cluster is present, getting the associated information for cluster number: $cluster", 1);
	    # pull all the associated accnos
	    $accnoList = $maindb->getValue("SELECT similar_accnos FROM $its WHERE accno = \'".$clusterKey->display_id()."\'");
	} else {
	    $cluster = $maindb->getValue("SELECT cluster + 1 from $its ORDER BY cluster DESC limit 0,1");
	    $log->print("Creating cluster number: ". $cluster ." with sequence ". $clusterKey->display_id);
	    $accnoList = $clusterKey->display_id() ."; ";
	}
	
	# Do a little bit of prep work to create information for the cluster
	foreach my $element (@clusteredSeqs) {
	    my $display = $element->display_id();
	    #make sure we have a catch on the old data that could be in there
	    if($accnoList !~ /$display/ ){
		$accnoList = $accnoList . $element->display_id()."; ";
	    }
	}
	# remove the trailing "; "
	$accnoList =~ s/\W*$//;
	
	my @clusteredAccnos = split('; ', $accnoList);
	$num = scalar(@clusteredAccnos);
	
	# Store the information
	# because the first one needs to still be flagged as important to clustering
	$log->print("Storing the cluster into the database", 1);
	my $key = shift(@clusteredAccnos);
	$maindb->do("UPDATE $its SET num_associated=\'$num\', cluster=\'$cluster\', similar_accnos=\"$accnoList\", percent = \"$percents[0]\", cluster_key=1 WHERE accno = \"$key\"") or die "Updating the clusterKey failed.\n";
	# the rest can be now marked as already clustered
	my $i = 1; # 0 because of cluster key wasn't on the percent thing
	foreach my $element (@clusteredAccnos) {
	    my $updateNonKeys = "UPDATE $its SET num_associated=\'$num\', cluster=\'$cluster\', similar_accnos=\'".$accnoList."\', percent = \'$percents[$i]\', cluster_key = 0 WHERE accno = \'".$element."\'";
	    $maindb->do($updateNonKeys) or die "Updating the rest of the cluster failed.\n";
	    $i++;
	}
    }
    $log->print("Analyzing the data collected and storing it into $clusterTable");
    my $maxCluster = $maindb->getValue("SELECT MAX(cluster) FROM $its");
    for ( my $i = 1; $i < ($maxCluster+1); $i++){
	my $newDate = $maindb->getValue("SELECT MAX(date) FROM $its WHERE cluster = $i");
	my $refAccno = $maindb->getValue("SELECT accno FROM $its WHERE cluster = $i AND cluster_key = 1");
	my $refSpecies = $maindb->getValue("SELECT species FROM $its WHERE cluster = $i AND cluster_key = 1");
	my $size = $maindb->getValue("SELECT COUNT(accno) FROM $its WHERE cluster = $i");
	my $avg = $maindb->getValue("SELECT AVERAGE(percent) FROM $its WHERE cluster = $i");
	my $lineage = $maindb->getValue("SELECT taxonomy FROM $its WHERE cluster = $1 AND cluster_key = 1");
	my $precluster = "None";
	pop(@cladeList);
	foreach my $class (@cladeList) {
	    chomp($class);
	    if($class =~ $lineage){
		$precluster = $class;
	    }
	}
	$maindb->do("INSERT INTO $clusterTable (cluster,date,ref_accno,ref_species,size,avg_sim,precluster) VALUE (\"$newDate\", \"$refAccno\", \"$refSpecies\", \"$size\", \"$avg\", \"$precluster\") ON DUPLICATE KEY UPDATE date = \"$newDate\", size = \"$size\", avg_sim = \"$avg\", precluster = \"$precluster\"");
    }
    
    $log->print("Now aligning clusters from tansley");
    #clusters start counting at 1, sorry
    my $clusterCount = $maindb->getValue("SELECT cluster + 1 from $its ORDER BY cluster DESC limit 0,1");
    for(my $i = 1; $i <= $clusterCount; $i++){
	my $clusterAlign = $maindb->getAlignObject($its, "original_seq", "cluster = $i");
	# write the alignment
	my @alignedCluster = $clusterAlign->each_seq();
	if (scalar(@alignedCluster) != 1){
	    my $clusterFile = "tempfiles/clusters" . $i .".fasta";
	    my $clusterOut = "tempfiles/clusters"  . $i .".fasta";
	    my $outCluster = Bio::AlignIO->new(-file => ">$clusterFile", -format => 'fasta');
	    $outCluster->write_align($clusterAlign);
	    
            # Align the sequences using MAFFT
            # This could maybe use some better control and wrapping but it works for right now
	    $log->print("Aligning the sequences.",1);
	    my $cmd = "mafft --auto ". $clusterFile ." > " . $clusterOut;
	    system($cmd);

            # get the alignment and store it mysql table
	    my $inCluster  = Bio::SeqIO->new(-file  => $clusterOut , -format => 'fasta');
	    
	    while ( my $obj = $inCluster->next_seq() ) {
		my $query = "UPDATE $its SET align_seq = \"". $obj->seq ."\" WHERE accno=\"". $obj->display_id ."\"";
		$maindb->do("$query");
	    }
	} else {
	    my $query = "UPDATE $its SET align_seq = \"". $alignedCluster[0]->seq ."\" WHERE accno=\"". $alignedCluster[0]->display_id ."\"";
	    $maindb->do("$query");
	}

    }

    $log->print("Tansley storing and sorting finished.",1);

} else {

    $log->print("There were no new records from genbank to add.");

}

mor::Utilities::clearProcessLocks("tansley");


sub cluster {
    my $percentSimilarity = shift;
    my $refSeq = shift;
    my $clusterRef = shift;
    my @percents = (100);
    my @accepted = @{$clusterRef};
    my @clusteredSeqs = ($refSeq);

    
    for(my $i = 0; $i < scalar(@accepted); $i++) { # For each record after this record all the previously compared ones                                                                                                              
	if($i % 500 == 0) { # If the number of sequences left is divisible by 500                                                                                                                                                    
	    $log->print("Still clustering.  We're on $i of " . scalar(@accepted) . " with ". scalar(@clusteredSeqs) ." in the cluster.", 1);
	}

	my $outFile = "cluster/cluster.fasta";
	my $mafftOut = "cluster/clusterALN.fasta";

	my $out = Bio::SeqIO->new(-file => ">$outFile", -format => 'fasta');
	$out->write_seq($clusteredSeqs[0]);
	$out->write_seq($accepted[$i]);

	system("mafft --quiet $outFile > $mafftOut"); # Run MAFFT on the FASTA file                                                                                                                                                  
	my $in = Bio::AlignIO->new(-file => $mafftOut , -format => "fasta"); # Create a Bio::AlignIO object
	my $aln = $in->next_aln(); # Get the alignment info                                                                                                                                                                          
	my $percent = $aln->percentage_identity(); # Get the percentage similarity                         
	if($percent > $minimum_similarity_value) { # If the sequences are above the minimum similarity value
	    my $goodSeq = splice(@accepted, $i, 1); # Pull this sequence (arrayref) out of the seqs array
	    push @clusteredSeqs, $goodSeq; # Add this sequence (arrayref) to the array of clustered sequences
	    push @percents, $percent; # Push this percentage value onto the percents array
	    $i--; # We need to check this element position again, since we removed what was there and now there's something new it that position
	}
    }

    return (\@accepted, \@clusteredSeqs, \@percents);

}
