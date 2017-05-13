package mor::Screening;

use strict;

use DBI;
use Bio::Factory::EMBOSS;
use Bio::AlignIO;
use Bio::SeqIO;

use mor::Log;
use mor::Settings;
use mor::DatabaseIO;

use vars qw( $log *settings $maindb);

# The global log object and settings object are created dynamically at the start
# of the program and assumed to be valid by the time this module is used. I.e.,
# this module is not valid outside of the context of the full program.
*log = \$mor::Log::log;
*settings = \%mor::Settings::settings;
*maindb = \$mor::DatabaseIO::maindb;


sub reasonForRejection {
    my $seq = shift;
    my $HMMprofile = shift;
    
    #assume HMMER failed	
    my $reason = 2;
    my $Ncount = 0;
    my $currentSeq = $seq->seq;	
    $currentSeq =~ s/(N)/$Ncount++;'N'/eg;
    $seq->seq($currentSeq);
    my $percent = $Ncount/$seq->length();
    
    #Sequence has too many N's? (> 3 percent)
    if($percent > 3){
	$reason = 1;
    } elsif ($seq->length() < 600) {
	$reason = 3; 
    } elsif ($seq->length() > 1450) {
	$reason = 4; 
    }
    if($HMMprofile and $reason == 2){
	if (HMMERtester($seq, $HMMprofile) == 1) {
	    $reason = ""; 
	}
    } elsif (!$HMMprofile) {
	$reason = "";
    }
    return $reason;
}

sub sortSequences(\@$$$$) {
    my (@alignI) = @{(shift)};
    my $table = shift;
    my $primerA = shift;
    my $primerB = shift;
    my $HMMprofile = shift;
    
    foreach my $seq (@alignI){
	
	my $currentSeq = $seq->seq();
	# all that is not CcTtAaGgNn gets turned into an N
	# most programs can deal with ambiguity codes
	# $currentSeq =~ tr/AaGgCcTtNn/N/cd; 
	
	# force the sequence to be uppercase
	$currentSeq = uc($currentSeq);
	$seq->seq($currentSeq);
	
	# remove all initial and trailing N's, such as Langer's sequences
	$currentSeq =~ s/^N+//m;
	$currentSeq =~ s/N+$//m;

	my $reason = 1;
	if(!$seq->seq){
	    $log->print("The sequence is empty!", 1);
	    next;
	}

	my @info = getInfo($seq);

	if ($reason ne ""){
	    my $normSeq = $seq;
	    if($primerA && $primerB){
		my $trimmed = vectorFilter($normSeq->seq(),$primerA,$primerB);
		$normSeq->seq($trimmed);
	    }	    
	    $log->print("trying screening for [".$seq->accession_number ."]",1);
	    $reason = reasonForRejection($normSeq, $HMMprofile);
	    if($reason eq ""){
		$log->print("Sequence of [". $seq->accession ."] worked!",1);
		$info[12]=0;
		$info[13]=$seq->seq;
		$info[3] = 1;
	    }
	}
	if ($reason ne "") {
		$log->print("Sequence [" . $seq->accession_number . "] did not pass screening; trying reverse-complement", 1);
		
		my $revSeq = $seq->revcom();
		if($primerA && $primerB){
		    my $trimmed = vectorFilter($revSeq->seq(),$primerA,$primerB);
		    $revSeq->seq($trimmed);
		}
		$reason = reasonForRejection($revSeq, $HMMprofile);
		if ($reason eq "") {
			$log->print("Reverse-complement of [" . $seq->accession_number . "] worked!", 1);
			$info[12] = 1;
			$info[13] = $revSeq->seq;
			$info[3] = 1;
		}
	}
	if ($reason ne ""){
		$log->print("Sequence [" . $seq->accession_number . "] did not pass reverse-complement; trying compliment", 1);
		#reversing the reverse compliment gives compliment
		my $complimentSeq = $seq->revcom;
		my $dummy = $complimentSeq->seq;
		$complimentSeq->seq(reverse($dummy));
                if($primerA && $primerB){
                    my $trimmed = vectorFilter($complimentSeq->seq(),$primerA,$primerB);
                    $complimentSeq->seq($trimmed);
                }
		my $reason = reasonForRejection($complimentSeq, $HMMprofile);
		if($reason eq ""){
			$log->print("Compliment of [". $seq->accession_number ."] worked!", 1);
			$info[12] = 2;
			$info[13] = $complimentSeq->seq;
			$info[3] = 1;
		}
	}
	if ($reason ne ""){	
		$log->print("Sequence [" . $seq->accession_number . "] did not pass compliment; trying reverse", 1);
		my $reverseSeq = $seq;
		$reverseSeq->seq(reverse($reverseSeq->seq()));
                if($primerA && $primerB){
                    my $trimmed = vectorFilter($reverseSeq->seq(),$primerA,$primerB);
                    $reverseSeq->seq($trimmed);
                }
		$reason = reasonForRejection($reverseSeq, $HMMprofile);
		if($reason eq ""){
			$log->print("Reverse of [". $seq->accession_number ."] worked!", 1);
			$info[12] = 3;
			$info[13] = $reverseSeq->seq;
			$info[3] = 1;
		}
	}
	
	if ($reason eq "") {
	    # passed HMMER successfully
	    # HMMER said: this is a homobasidiomycete
	    $log->print("Sequence for [".$seq->accession_number."] passed HMMER test.",1);
	}
	else{
	    @info = getInfo($seq);
	    $log->print("Sequence for [".$seq->accession_number."] failed HMMER test.",1);
	}
	my $success;
	# This could potentially be expanded on to take into account the two other possibilities of hmmer+!primers and !hmmer+primers
	if($primerA && $primerB){
	    my $insert = "INSERT INTO ". $table." (accno,species,date,accepted,version,clade,HRI,aftol_id,db_xref,isolate,culture_collection,specimen_voucher,strain,taxonomy,environmental,GI,seq_manip,accepted_seq,original_seq, reason) ";
	    my $value = "VALUES (\'$info[0]\',\'$info[1]\',\'$info[2]\',\'$info[3]\',\'$info[4]\',\'$info[5]\',\'$info[6]\',\'$info[7]\',\'$info[8]\',\'$info[9]\',\'$info[16]\',\'$info[17]\',\'$info[18]\',\'$info[15]\',\'$info[10]\',\'$info[11]\',\'$info[12]\',\'$info[13]\',\'$info[14]\',\'$reason\') ";
	    my $update = "ON DUPLICATE KEY UPDATE date = \'$info[2]\', version = \'$info[4]\', aftol_id = \'$info[7]\', isolate = \'$info[9]\', culture_collection = \'$info[16]\', specimen_voucher = \'$info[17]\', strain = \'$info[18]\', taxonomy = \'$info[15]\', GI = \'$info[11]\'";
	    my $query = $insert . $value . $update;
	    $success = $maindb->do($query);
	} else {
	    my $insert = "INSERT INTO ". $table." (accno,species,date,accepted,version,clade,cluster_key,num_associated,aftol_id,db_xref,isolate,culture_collection,specimen_voucher,strain,taxonomy,environmental,GI,original_seq,reason) ";
	    my $value = "VALUES (\'$info[0]\',\'$info[1]\',\'$info[2]\',\'$info[3]\',\'$info[4]\',\'$info[5]\',\'1\',\'0\',\'$info[7]\',\'$info[8]\',\'$info[9]\',\'$info[16]\',\'$info[17]\',\'$info[18]\',\'$info[15]\',\'$info[10]\',\'$info[11]\',\'$info[14]\',\'$reason\') ";
	    my $update = "ON DUPLICATE KEY UPDATE date = \'$info[2]\', version = \'$info[4]\', aftol_id = \'$info[7]\', isolate = \'$info[9]\', culture_collection = \'$info[16]\', specimen_voucher = \'$info[17]\', strain = \'$info[18]\', taxonomy = \'$info[15]\', GI = \'$info[11]\'";
	    my $query = $insert . $value . $update;
	    $success = $maindb->do($query);
	}
	if($success){
	    $log->print("Accession number [". $seq->accession_number ."] successfully added to the database.",1);
	}
	else{
	    $log->print("Insertion for [". $seq->accession_number ."] failed.",1);
	}	
    }
    
}

sub getInfo{
    my $seq = shift;
    
    my @info;
    my $speciesObj = $seq->species;
    my @dates = $seq->get_dates;
    my @datesArray = split(/-/,$dates[0]);
    
    if ($datesArray[1] eq "JAN") { $datesArray[1] = "01"; }
    elsif ($datesArray[1] eq "FEB") { $datesArray[1] = "02"; }
    elsif ($datesArray[1] eq "MAR") { $datesArray[1] = "03"; }
    elsif ($datesArray[1] eq "APR") { $datesArray[1] = "04"; }
    elsif ($datesArray[1] eq "MAY") { $datesArray[1] = "05"; }
    elsif ($datesArray[1] eq "JUN") { $datesArray[1] = "06"; }
    elsif ($datesArray[1] eq "JUL") { $datesArray[1] = "07"; }
    elsif ($datesArray[1] eq "AUG") { $datesArray[1] = "08"; }
    elsif ($datesArray[1] eq "SEP") { $datesArray[1] = "09"; }
    elsif ($datesArray[1] eq "OCT") { $datesArray[1] = "10"; }
    elsif ($datesArray[1] eq "NOV") { $datesArray[1] = "11"; }
    elsif ($datesArray[1] eq "DEC") { $datesArray[1] = "12"; }
    
    my $date = "$datesArray[2]$datesArray[1]$datesArray[0]";
    $info[0] = $seq->accession_number;
# The following line didn't work for screening for redundency
#    $info[1] = $speciesObj->binomial;
    my @organism = getTagFrom("organism", $seq);
    $info[1] = $organism[0];
    $info[1] =~ s/[\'|\"|\`]//g;
    $info[2] = $date;
    $info[3] = "0"; # "accepted" field
    $info[4] = $seq->seq_version;
    $info[5] = "incertae sedis"; # "clade" field
    $info[6] = "1"; # "HRI" or "num_associated" field
    my @tags = findWord("aftol",$seq);
    my @values = getTagFrom($tags[0], $seq);
    $values[0] =~ s/\D//g;
    if(!$values[0]){
	$values[0] = "NULL";
    }
    $info[7] = $values[0];
    $info[8] = $speciesObj->ncbi_taxid;
    
    my @isolate = getTagFrom("isolate", $seq);
    if(!$isolate[0]){
	$isolate[0] = "NULL";
    }
    $info[9] = $isolate[0];
    $info[10] = getTagFrom("environmental_sample",$seq);
    $info[11] = $seq->primary_id;
    $info[12] = "0"; # "seq_manip" field
    $info[13] = "0"; # "accepted_seq" field
    $info[14] = $seq->seq;
    
    my @class = $speciesObj->classification;
#    $info[1] = $class[0];
    @class = reverse(@class);
    my $lineage;
    foreach my $element (@class) {
	$lineage .= $element . ":";
    }
    $lineage = substr($lineage, 0, -2);
    $info[15] = $lineage;	
    $info[15] =~ s/[\'|\"|\`]//g;
    my @cc = getTagFrom("culture_collection", $seq);
    if(!$cc[0]){
	$cc[0] = "NULL";
    }
    $info[16] = $cc[0];
    
    my @sv = getTagFrom("speciment_voucher", $seq);
    if(!$sv[0]){
	$sv[0] = "NULL";
    }
    $info[17] = $sv[0];
    
    my @strain = getTagFrom("strain", $seq);
    if(!$strain[0]){
	$strain[0] = "NULL";
    }
    $info[18] = $strain[0];
    
    
    return @info;
    
}

sub sequenceSimilarity(\@$) {
# This function could be implemented more accurately if we only work with the genus from the classifcation ?
# This could also be be rewritten so that it is passed an array of seq objects and a cut-off value, and returns two arrays
# one grouped above that cut off, another grouped below
    my (@seq) = @{(shift)};
    my $table = shift;
    
    # temporary storage file
    my $outfile= mor::Utilities::getUniqueFile('tempfiles/pwout');
    my $tempfile= mor::Utilities::getUniqueFile('tempfiles/pwin');    

    
    for (my $i=1; $i< scalar(@seq); $i++) {
	open(SEQS, ">$tempfile" );
	print SEQS ">".$seq[0]->accession_number." \n ".$seq[0]->seq." \n";
	print SEQS ">".$seq[$i]->accession_number." \n ".$seq[$i]->seq." \n";
	close(SEQS);

	system("mafft-linsi --quiet $tempfile > $outfile");
	
	# BioPerl AlignIO
	my $in = new Bio::AlignIO ( -file => $outfile, -format => 'fasta' );
	my $aln = $in->next_aln();
	
	# these pertains to the pairwise alignment
	my $percent  = $aln->percentage_identity();
	my $length   = $aln->length();
	my $residues = $aln->no_residues();
	
	# some statistics
	$log->print("Length: $length", 1);
	$log->print("Similarity: $percent (some $residues in-alignment mismatches)", 1);
	
	if ($percent > 99.5) {
	    # too similar. reject it.
	    $log->print("The incoming sequence is >99.5% identical to an in-db sequence and will be marked as redundant.", 1);
	    # UPDATE THE DATABASE
	    # Decrease the HRI of the redundant sequence to 0
	    my $query = "UPDATE $table SET HRI = 0 WHERE accno=\"".$seq[$i]->accession_number."\" ";
	    $maindb->do($query);
	    # Increase the HRI of the passed sequence by 1
	    $query = "UPDATE $table SET HRI=HRI+1 WHERE accno=\"". $seq[0]->accession_number ."\"";
	    $maindb->do($query);
	    # Add the accno of the redundant to the redundant field of the passed seq
	    $query = "UPDATE $table SET redundant_list=\"". $seq[0]->accession_number ."\"  WHERE accno=\"". $seq[$i]->accession_number ."\" ";
	    $maindb->do($query);
	    # Add the accno of the passed to the redundant field of the redundant seq
	    $query = "UPDATE $table SET redundant_list = CONCAT_WS(\',\',redundant_list,\'". $seq[$i]->accession_number ."\')  WHERE accno=\'". $seq[0]->accession_number ."\' ";
	    $maindb->do($query);
	}
    }
}


sub vectorFilter {
    my $sequence = shift;
    my $primerA = shift;
    my $primerB = shift;
    my $vector_mismatch_percentage = $settings{"vector_mismatch_percentage"};
    
    my $temp_in = "vectorFilter_in";
    my $temp_out = "vectorFilter_out";
    open (IN, ">$temp_in") or die "error in vectorFilter, Screening.pm: $!";
    print IN $sequence;
    close IN;
    
    open (OUT, ">$temp_out") or die "error in vectorFilter, Screening.pm: $!";
    
    my $factory = new Bio::Factory::EMBOSS;
    my $vectorstripapp = $factory->program('vectorstrip');
    
    
    # Variables to store the sequences while we work on them
    my @out_file;
    my $trimmed;
    
    # Set up the input for vectorstrip. This apparently needs to use a file for both input and output of sequences
    my %vectorstrip_object = (   -sequence       => $temp_in,
				 -linkera        => $primerA,
				 -linkerb        => $primerB,
				 -vectorfile     => "no",
				 -mismatch       => $vector_mismatch_percentage,
				 -outseq         => $temp_out);
    
    # Run vectorstrip
    $vectorstripapp->run(\%vectorstrip_object);
    close OUT;
    # Get the trimmed sequence back
    open (OUT1, "<$temp_out") or die "another error in vectorFilter(): $!";
    @out_file = <OUT1>;
    close OUT1;
    $trimmed = join("",@out_file);
    
    # sometimes vectorstrip generates blank files. assuming for the moment
    # that this means that it had some kind of problem processing the file,
    # so we might as well ignore it and conclude that it couldn't strip
    # the sequence
    
    # Get rid of all the info
    $trimmed =~ s/\>_from_\d+_to_\d+//g;
    
    # Clean up any spaces or newlines that may be left.
    $trimmed =~ s/\s//g;
    
    if ($trimmed eq "") {
	return $sequence;
    }
    # Hand back the trimmed (or maybe not?) sequence
    return $trimmed;
}

sub HMMERtester {
    my $seq = shift;
    my $profile = shift;
    
    # name of the (systemwide) HMMER executable
    my $HMM = "hmmpfam";
    my $passedTest = 1;
    
    # HMMER needs the seuqence as a FASTA file
    my $tempfile="HMM.fasta";
    
    my $out = Bio::SeqIO->new(-file => ">$tempfile", -format => 'fasta');
    $out->write_seq($seq);
    
    my $HMMERoutfile="tempfiles/hmmer_results.txt";
    
    # for info about parameters, see HMMER docs. T = score, E = e value
    
    my $cmd = "$HMM -n -A 0 -E 10 -T 0 " .  $profile . " " . $tempfile . " > " . $HMMERoutfile;
    system($cmd);
    #-A 0 -E 0 -T 900
    open (HMMEROUT, "<$HMMERoutfile") or die "error: $!";
    my @file = <HMMEROUT>;
    
    # if the line contains this, it didn't pass HMMER
    my $noHitsString = "thresholds";
    
    if ($file[16] =~ "$noHitsString"){
	return 0;
    } elsif ($file[21] =~ "$noHitsString"){
	return 0;
    }
    $log->print("Yes, sequence [" . $seq->accession_number . "] passed HMMER test.", 1);
    return 1;

}

sub findWord{

    my $word = shift;
    my $seq = shift;
    my @location;

    foreach my $featureObj ( $seq->get_SeqFeatures) {
        foreach my $tag ( $featureObj->get_all_tags ){
            foreach my $value ($featureObj->get_tag_values($tag)){
                if ($value =~ m/$word/i){
                    push(@location, $tag);
                }
            }
        }
    }
    return @location;
}

sub getTagFrom{
    my $tagWanted = shift;
    my $seq = shift;
    my @value;

    foreach my $featureObj ( $seq->get_SeqFeatures) {
        if($featureObj->has_tag($tagWanted)) {
            foreach my $tag ( $featureObj->get_tag_values($tagWanted)){
                push(@value, $tag);
            }
        }
    }
    return @value;
}

1;


__END__


=head1 NAME
                                                                                                                    
mor::Screening - Analysis the incoming sequences from genbank

=head1 SYNOPSIS

	use DBI;
	use Bio::Factory::EMBOSS;
	use Bio::AlignIO;
	use Bio::SeqIO;
	use mor::Log;
	use mor::Settings;
	use mor::DatabaseIO;
	use vars qw( $log *settings $maindb);

=head1 DESCRIPTION

This module handles all the sequences from genbank retrieves the information and stores it into a MySQL table

=head2 Methods

=over 6

=item B<reasonForRejection>

 Title   : reasonForRejection
 Usage   : my $reason = mor::Screening::reasonForRejection($seq, $HMMprofile)
 Function: Determines if a sequence is valid and if not why it is not valid
 Returns : The reason for the sequence to be rejected
 Args    : a Bio::Seq object
	   the location and name of a hmm profile

=item B<sortSequences>
                                                                                                                 
 Title   : sortSequences
 Usage   : mor::Screening(@seqs, $table, $primerA, $primerB, $HMMprofile);
 Function: Sorts through an array of Bio::Seq::RichSeq objects and stores information in the object.
	   It also analysis the sequence to see if it has been deposited into genbank as 
	   reverse-complemented (1), complemented (2) or reversed (3), (0 if unmanipulated).
	   It does this by using a hmm profile and vector stripping the sequences (if primers are provided).
 Returns : Nothing
 Args    : an array of Bio::Seq::RichSeq objects (preferably from a genbank file)
	   the MySQL table that the sequences will be stored in.
	   two of the primers used for vector stripping (pass an empty string if you do not want this to happen)
	   the file name for a hmm profile (pass an emptry string if you do not want this to happen)

=item B<getInfo>

 Title   : getInfo
 Usage   : my @info = mor::Screening::getInfo($seq)
 Function: Retrieves all the information a Bio::Seq::RichSeq object and returns it in an array
 Returns : an array of a lot of information from the passed object
 Args    : a Bio::Seq::RichSeq object

=item B<sequenceSimilarity>

 Title   : sequenceSimilarity
 Usage   : my @info = mor::Screening::seqeuenceSimilarity(@seq, $table)
 Function: Retrieves all the information a Bio::Seq::RichSeq object and returns it in an array after
           marking the sequences as redundant in the passed table.
 Returns : nothing
 Args    : an array Bio::Seq::RichSeq object
           the table in wich the seqeuences will be marked as redundant 

=item B<vectorFilter>

 Title   : vectorFilter
 Usage   : my $trimmed = mor::Screening::seqeuenceSimilarity($sequence, $primerA, $primerB)
 Function: vector strips a sequence in string from (ATGC) in $sequence, look at Bio::Factory:EMBOSS
 Returns : the trimmed sequence
 Args    : a string representation of a sequence
           two string representations of the primers

=item B<HMMERtester>

 Title   : HMMERtester
 Usage   : my $result = mor::Screening::HMMERtester($seq, $profile)
 Function: Tests to see if the passed Bio::Seq object matches the hmm profile
 Returns : 1 if the sequence matches 0 if not.
 Args    : a Bio::Seq object
           the file name of a hmm profile

=item B<findWord>

 Title   : findWord
 Usage   : my @tag = mor::Screening::findWord($word, $seq)
 Function: Finds where a word is located in a Bio::Seq::RichSeq object
 Returns : an array of the tag names where it can be found
 Args    : a string
           a Bio::Seq::RichSeq object

=item B<getTagFrom>

 Title   : getTagFrom
 Usage   : my @values = mor::Screening::getTagFrom($tag, $seq)
 Function: Gets the values from a bioperl feature object (see the bioperl documentation)
 Returns : an array of the values from that tag (generally only the first element)
 Args    : the name of the tag
           a Bio::Seq::RichSeq object

=back

=head1 AUTHOR

The I<mor> Team - L<http://mor.clarku.edu/morTeam.php>

=cut
