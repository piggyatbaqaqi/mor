package mor::Acquisition;

use strict;
use DBI;

use Bio::DB::GenBank;
use Bio::DB::Query::GenBank;
use Bio::SeqIO;

use mor::Log;
use mor::Settings;
use mor::DatabaseIO;
use mor::Utilities;

use vars qw( $log *settings $maindb );

# The global log object and settings object are created dynamically at the start 
# of the program and assumed to be valid by the time this module is used. I.e.,
# this module is not valid outside of the context of the full program.
*log = \$mor::Log::log;
*maindb = \$mor::DatabaseIO::maindb;
*settings = \%mor::Settings::settings;



sub queryGB {
    
    my ($table,$saveFilename,$query_string,$mindate,$maxdate) = @_;
    
    
    if ($query_string eq "") { 
	$query_string = $settings{'query_string'}; 
    }
    
    $log->print("Using query string: $query_string", 1);
    
    # if we haven't specified a minimum date to use...
    if ($mindate eq "") {
	
	# the last run file stores the last day mor was run, for the sake of
	# limiting the GenBank query to only what is presumably new for mor.
	# not critical if missing or tampered with, but it directly influences
	# the query
	
	# query mysql to see if there are any entries in the table
	my $count = $maindb->getValue("SELECT COUNT(accno) from $table");
	if($count == 0) {
	    # set lastrun to 6 months ago and inform the user of a strategy to get older sequences
	    $log->print("Setting last run to 6 months ago, replace the download file with a genbank file to get older sequences", 1);
	    $mindate = mor::Utilities::getDateToday();
	    my @tempDate = split("\/",$mindate);
	    $tempDate[1] = (($tempDate[1] + 6) %12);
	    if(($tempDate[1] + 6) > 12) {
		$tempDate[0] = $tempDate[0] - 1;
	    }
	    $mindate = "$tempDate[0]\/$tempDate[1]\/$tempDate[2]";
	} else {
	    # get date from table ordered by desc/asec and limit to 1
	    $mindate = $maindb->getValue("SELECT date FROM $table ORDER BY date DESC LIMIT 0,1");
	    $mindate =~ s/-/\//g;
	}
	$log->print("Mindate is now $mindate",1);
    }
    if ($maxdate eq "") {
	$maxdate = mor::Utilities::getDateToday();
    }
    $log->print("Maxdate is now $maxdate",1);

    $log->print("Last fetched sequences on [" . $mindate . "]. Searching until [" . $maxdate . "]", 1);
    
    # use BioPerl to do a Genbank query. please refer to BioPerl docs.
    # create a BioPerl GB query object, feeding in the query limiters
    # Genbank is a bit of a diva so we have to work around it a bit.
    # It only gives us a small subset of the number of sequences that we request because of various timeout issues.
    
    my $gb = new Bio::DB::GenBank;
    my @sequences;
    
    # prime all the variables in case we need to loop
    my $queryObj = Bio::DB::Query::GenBank->new(-db => 'nucleotide',
						-query => $query_string,
						-mindate => $mindate,
						-maxdate => $maxdate,
						# -maxids=>5000,
						# -retrievaltype => 'tempfile');
						);
    
    
    my $out= Bio::SeqIO->new( -file=>">" . $saveFilename, -format=>'GenBank');
    my $totalSaved = 0;
    my $count = $queryObj->count;
    $log->print("Downloading $count new sequences.", 1);
    
    $log->print("Downloading could take a while, please be patient.", 1);
    my $seqio = $gb->get_Stream_by_query($queryObj);
    
    while( my $seq =  $seqio->next_seq ) {
	$out->write_seq($seq);
	$log->print("Saving ". $seq->accession_number . " to $saveFilename");
	push(@sequences, $seq);
    }
    # this next step takes into account how finicky genbank can be
    while( scalar(@sequences) < $count ) {
	$log->print("Waiting 30 seconds before downloading again", 1);
	sleep 30;
	my $date;
	if(scalar(@sequences != 0)) {

	    # pull off all the ones that were deposited on the last day
	    my @older = $sequences[-1]->get_dates;
	    my @younger = $sequences[-2]->get_dates;
	    
	    while($older[0] eq $younger[0]) {
		pop(@sequences);
		@older = $sequences[-1]->get_dates;
		@younger = $sequences[-2]->get_dates;
	    }
	    
	    pop(@sequences);
	    
	    my @dates = $sequences[-1]->get_dates;
	    
	    #convert to a genbank recognizable string
	    my $size = scalar(@dates);
	    $size--;
	    my @datesArray = split(/-/,$dates[$size]);
	    
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
	    $date = "$datesArray[2]\/$datesArray[1]\/$datesArray[0]";

	} else {
	    $date = $mindate;
	}
	$log->print("Going to try downloading again with new min date of $date", 1);
	
	my $renewQuery = Bio::DB::Query::GenBank->new(-db => 'nucleotide',
						      -query => $query_string,
						      -mindate => $date,
						      -maxdate => $maxdate,
						      );
	
	my $outNew= Bio::SeqIO->new( -file=>">>" . $saveFilename, -format=>'GenBank');
	my $newCount = $renewQuery->count;
	$log->print("Attempting to download $newCount more new sequences.", 1);
	
	$log->print("Downloading could take a while, please be patient.", 1);
	$seqio = $gb->get_Stream_by_query($queryObj);
	
	while( my $seq =  $seqio->next_seq ) {
	    $out->write_seq($seq);
	    $log->print("Saving ". $seq->accession_number . " to $saveFilename");
	    push(@sequences, $seq);
	}
    }

    $totalSaved = scalar @sequences;
    $log->print("Total number of records found: $count", 1);
    $log->print("Total number of records saved: $totalSaved", 1);
    if($totalSaved != $count) {
	die "The total number of saved sequences did not match how many genbank told us there should be.";
    }
    $log->print("Recording today's date as the date of the last run.", 1);
    return @sequences;
}

sub tryQueryGB {
    my $mindate = shift;
    my $maxdate = shift;
    my $queryString = shift;
    my $table = shift;
    my $answer = shift;
    
    my $dataStorage = "download/". $table ."_download.txt";
    my @sequences;	
    my $answerIsInvalid = 1;
    
    if(!$answer && -e$dataStorage){
	$log->print("Going to use previously acquired data.", 1);
	$answer = "u";
    }
    if(!$answer && !(-e$dataStorage)){
	$log->print("Going to download sequences from genbank", 1);
	$answer = "i";
    }
    if($answer eq "u" or $answer eq "i"){
	$answerIsInvalid = 0;
    }
    while($answerIsInvalid && -e$dataStorage){
	$log->print("Previous data exists, what do you want to do?\n",1);
	$log->print("(u)se or (i)gnore?\nEnter your choice:",1);
	$answer = <>;
	chomp($answer);
	if($answer ne "u" and $answer ne "i"){
	    $log->print("That is not a valid answer.",1);
	}
	if($answer eq "u" or $answer eq "i"){
	    $answerIsInvalid = 0;
	}
    }
    
    if($answer eq "u"){
	my $stream = Bio::SeqIO->new(-file => "$dataStorage", -format => 'GenBank');
	while ( my $seq = $stream->next_seq() ) {
	    push(@sequences, $seq);
	}
	$log->print("Sequences have been loaded.", 1);
	return @sequences;
    }
    if($answer eq "i" or !$answer){
	open (DATA, ">$dataStorage") or die "error: $!";
	@sequences = queryGB($table,$dataStorage, $queryString, $mindate, $maxdate);
	if ($sequences[0]) { 
	    $log->print("Sequences retrieved successfully. Continuing to screening stage.", 1);
	    return @sequences;
	}
	else{
	    $log->print("There has been an error writing the downloaded sequences.", 1);
	    exit;
	}
    }
}

1;


__END__

=head1 NAME

mor::Acquisition - Genbank downloader

=head1 SYNOPSIS

	use DBI;
	use Bio::DB::GenBank;
	use Bio::DB::Query::GenBank;
	use Bio::SeqIO;
	use mor::Log;
	use mor::Settings;
	use mor::DatabaseIO;
	use mor::Utilities;
	use mor::Acquisition;

=head1 DESCRIPTION
                                                                                                                    
This package deals with querying GenBank for (new) sequences matching
the query string and saving the data for later processing.
                                                                                                                    
=head2 Methods
                                                                                                                    
=over 6
                                                                                                                    
=item B<queryGB>
                                                                                                                 
 Title   : queryGB 
 Usage   : my @sequences = mor::Acquisition::queryGB($table, $saveFile, $queryString, $mindate, $maxdate)
 Function: Downloads the sequences using $queryString that were deposited in genbank between 
	   the dates: $mindate and $maxdate.  It then saves the files to $table in a mysql database
	   and to the file $saveFile.  It returns an array of the sequences for further processing.
 Returns : An array of Bio::Seq::RichSeq objects.
 Args    : name of the table for the sequences to be stored
	   name of the file that the sequences will be saved to
	   genbank query string
	   earliest date to retrieve from genbank
	   latest date to retrieve from genbank

=item B<tryQueryGB>

 Title   : tryQueryGB
 Usage   : my @seqs = mor::Acquisition::tryQueryGB($mindate,$maxdate,$queryString,$table,$answer)
 Function: Acts as a wrapper around the queryGB function, allows for the program to
	   choose between the previous downloaded file or to download more sequences
	   from genbank to be saved to the database
 Returns : An array of Bio::Seq::RichSeq objects.
 Args    : name of the table for the sequences to be stored
	   name of the file that the sequences will be saved to
	   genbank query string
	   earliest date to retrieve from genbank
	   latest date to retrieve from genbank
	   'u' or 'i' for (u)se the file or (i)gnore the file and download from genbank.

=back

=head1 AUTHOR

The I<mor> Team - L<http://mor.clarku.edu/morTeam.php>

=cut
