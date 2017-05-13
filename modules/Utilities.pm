package mor::Utilities;

use strict;

# we use this date class for intelligent handling of date conversions.
# no need for fancy format checking when this will take care of it,
# and in a modular fashion
use DBI;
use Bio::TreeIO;
use Bio::AlignIO;
use Class::Date qw(:errors date localdate gmdate now -DateParse -EnvC);
use mor::Settings;
use mor::DatabaseIO;
use mor::Log;

use vars qw(%settings $maindb $log);

*log = \$mor::Log::log;		
*settings = \%mor::Settings::settings;
*maindb = \$mor::DatabaseIO::maindb;

sub getDateToday {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime time;

	# month ranges from 0 to 11
	$mon = $mon + 1;

	$year +=1900;

	# e.g., 1 -> 01 for months
	if ($mon < 10) { $mon = '0' . $mon; }

	# same but for day. genbank compliance.
	if ($mday < 10) { $mday = '0' . $mday; }

	return "$year/$mon/$mday";
}

sub getUniqueFile {
	my $fileTemplate = shift;

	my @components = split(/\./, $fileTemplate, 2);
	my $fileBase = $components[0];
	my $ext = "";

	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime time;
	$year += 1900;
	$mon += 1;

	if ($mday < 10) { $mday = "0$mday"; }

	if ($fileBase ne "" && $fileBase !~ /\/$/) { $fileBase .= "-"; }

	$fileBase .= "$year-$mon-$mday";

	$components[1] = ".".$components[1];

	my $filename = $fileBase . $ext . $components[1];
	while(-e $filename) {
	    $ext = $ext + 1;
	    $filename = $fileBase . $ext . $components[1];
	}
	return $filename;
}

sub lockProcessOrDie {
	my $lockName = shift;
	my $lockPath = $settings{'concurrency'};
	my $lockFile = $lockPath . $lockName;

	# This function might fail in race conditions the chance of two instances of a script
	# starting split seconds apart from each other is minimal
	if (-e $lockFile) { 
	    die "Lock file is present [statefiles/mor_lock_".$lockFile."]. Another process is potentially running, manually check to see if mor is running and then remove the file to continue.\n";
	}

	# if we're here, we own the lock
	open(LOCKFILE, ">$lockFile");
	print LOCKFILE localtime time;
	close(LOCKFILE);

}

sub clearProcessLocks {

    my $filename = shift;
    $filename = $settings{'concurrency'} . $filename;
    if (-e $filename) {
	unlink($filename) or die "Unable to delete lock $filename: $!";
    }

}

sub storeSingleTree { 
  my $tree_string = shift;
  my $clade_name = shift;
  
  my $tree_query = "UPDATE clades SET asciitree = '$tree_string' WHERE cladename='$clade_name'";
  $maindb->query($tree_query);
}

sub introStatsHTML {
	my $treename = shift;
	my $treeDBName = shift;
	my $pageTitle = shift;
    my $html = "";

	# $$ means the current process PID number
	my $pid = $$;
	#The next two lines don't exactly give the sequences added in the last week (tends to be an over estimation) but it's accurate enough
	my $lastDate = $maindb->getValue("SELECT max(date) FROM ".$settings{'LSUtable'}."");
	my $addedSequences = $maindb->getValue("SELECT COUNT(accno) FROM ".$settings{'LSUtable'}." WHERE date>DATE_SUB(\"".$lastDate."\", INTERVAL 7 DAY)");
	my $numSpecies = $maindb->getValue("SELECT COUNT(accno) from ".$settings{'LSUtable'});

	my $date = $maindb->getValue('SELECT DATE_FORMAT(MAX(date), "%W %M %e, %Y") FROM trees WHERE id="'.$treeDBName.'"');
	$log->print("Found date: $date for tree: $treeDBName");
    $html = <<HTML;
<html>
  <head>
    <link rel="stylesheet" href="/css/mainstyle.css" type="text/css" />
    <title>$pageTitle</title>
  </head>
  <body class="mainTree">
  <div id="header"><a href="/index.php?mor=home"><img id="logo" src="/images/mor_logo.gif" alt="Mor" /></a></div>
  <ul>
    <li>Date Created: $date</li>
    <li>Tree Type: $treename</li>
    <li>Process ID: $pid </li>
    <li>New Sequences: $addedSequences</li>
    <li>Number of Sequences: $numSpecies</li>
    <ul>
      <li><h4 class="linkKey">What the links mean:</h4></li>
      <li>Species Name links you to the placement of the species in our accepted sequences list</li>
      <li>Accession Number links you to the NCBI entry</li>
      <li>Image links you to a google image search of the species name</li>
      <li>nj/mp/ml are links between the trees.</li>
      <ul>
        <li>nj = Neighbor-Joining Jacknife Tree</li>
        <li>mp = Maximum Parsimony Tree</li>
        <li>ml = Maximum Liklihood Tree</li>
      </ul>
    </ul>
  </ul>
  <pre>
HTML
# These lines are not being used in the tree, there is no more HRI...nor are there any *s. If it's needed just add it to the above code.
#	$html .= "* = Species present in constraints tree<br />";
#	$html .= "HRI = <em>n</em>; <span style=\"font-size: smaller;\">the H Rejection Index; <em>n</em> distinct Genbank sequences with the same species name<br /> have been tested and found to be over 99.5% identical to this sequence, and hence discarded.</span></p>";

	return $html;
}

sub HTMLformat{
	my $htmlName = shift;
	my $treeDBID = shift;
	my $pageTitle = shift;
	my $file = shift;
	my $linksRef = shift;
	my $table = shift;
	my @links = @{$linksRef};

	open(LOG, "<$file");
	my @lines = <LOG>;
	my @returnTree;
	my $NCBI = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nuccore&id=";
	my $google = "http://images.google.com/images?q=";
	my $seqLink = "http://mor.clarku.edu/index.php?mor=accepted&sort=s#";
#saving for CSS editing, will delete after css has been created
#	push(@returnTree, "<body><pre style=\"font-size: 6pt; font-family: Courier New\">");
	my $trim;
	until ($trim =~ m/\s*\/\-\s*/){
	        $trim = shift(@lines);
		if(!@lines){
			die "Not a valid paup log fed into HTMLformat";
		}
	}
	unshift(@lines, $trim);
	until ($trim =~ m/\s*\\\-\s*/){
		$trim = pop(@lines);
		if(!@lines){
			die "The top was OK the bottom was not, What happened?";
		}
	}
	push(@lines, $trim);
	my $decoration;
	foreach my $line (@lines){
	        chomp($line);
	        if ($line =~ s/([A-Z]+\w+)//){
	                my $name = $maindb->getValue("SELECT species FROM ".$table." WHERE accno = \"$1\"");
			$decoration = "$line<a name=\"$1\" href=\"$seqLink$1\">$name</a> <a href =\"$NCBI$1\">($1)</a> <a href=\"$google$name\">images</a>";
			my $accno = $1;
			foreach my $url (@links){
			    my $text = $url;
			    $text =~ s/\.html.*$//;
			    $text =~ s/^.*\///;
			    $decoration .= " <a href=\"$url$accno\">$text</a>";
			}
	                push(@returnTree, $decoration);
	        }else{
	                push(@returnTree, $line);
	        }
	}
	push(@returnTree, "</pre></body></html>");
	unshift(@returnTree, introStatsHTML($htmlName, $treeDBID, $pageTitle));
	return @returnTree;
}

sub archiveForWeb {
    my $success = 1;
    my $zipSuccess = 1;
    my $gmtime=gmtime;
    my @gmarray = split(' ', $gmtime);
    
    my $webPath = $settings{'webPath'};
    # Apr -> 04
    if ($gmarray[1] eq "Jan") { $gmarray[1] = "01"; }
    elsif ($gmarray[1] eq "Feb") { $gmarray[1] = "02"; }
    elsif ($gmarray[1] eq "Mar") { $gmarray[1] = "03"; }
    elsif ($gmarray[1] eq "Apr") { $gmarray[1] = "04"; }
    elsif ($gmarray[1] eq "May") { $gmarray[1] = "05"; }
    elsif ($gmarray[1] eq "Jun") { $gmarray[1] = "06"; }
    elsif ($gmarray[1] eq "Jul") { $gmarray[1] = "07"; }
    elsif ($gmarray[1] eq "Aug") { $gmarray[1] = "08"; }
    elsif ($gmarray[1] eq "Sep") { $gmarray[1] = "09"; }
    elsif ($gmarray[1] eq "Oct") { $gmarray[1] = "10"; }
    elsif ($gmarray[1] eq "Nov") { $gmarray[1] = "11"; }
    elsif ($gmarray[1] eq "Dec") { $gmarray[1] = "12"; }
    
    # 2007 -> 07
    $gmarray[4] =~ s/^\d\d(\d\d)$/$1/;
    
    # 3 -> 03
    if ($gmarray[2] < 10) { $gmarray[2] = "0$gmarray[2]"; }
    my $dirname = "$gmarray[4]-$gmarray[1]-$gmarray[2]";
    print $dirname;
    my $archivePath = getUniqueFile("$webPath/archive/$dirname");
    
    # make unique, "time-stamp-named" dir
    mkdir "$archivePath" or $success = 0; 
    chmod(0755,"$archivePath") or $success = 0;
    
    if (!$success) { 
	showError("Unable to create archive path [$archivePath]. Aborting archival process."); 
	return; 
    }
    
    # Get the ML tree object
    my $mlTree = $maindb->getTreeObj("trees", "PARTree");
    # Save the ML tree object
    my $mlOut = Bio::TreeIO->new(-file => ">$archivePath\/mltree.tree", -format => 'newick');
    $mlOut->write_tree($mlTree);
    # Get the NJ tree object
    my $njTree = $maindb->getTreeObj("trees", "NJTree");
    # Save the NJ tree object
    my $njOut = Bio::TreeIO->new(-file => ">$archivePath\/njtree.tree", -format => 'newick');
    $njOut->write_tree($njTree);
    # Get the ML tree object
    my $parTree = $maindb->getTreeObj("trees", "PARTree");
    # Save the ML tree object
    my $parOut = Bio::TreeIO->new(-file => ">$archivePath\/partree.tree", -format => 'newick');
    $parOut->write_tree($parTree);
    # Get the Easyrider aligned align object
    my $alnEasy = $maindb->getAlignObject("mor", "align_seq", "accepted = 1");
    # Save that object
    my $easyOut = Bio::AlignIO->new(-file => ">$archivePath\/easyriderAligned.fasta", -format => 'fasta');
    $easyOut->write_aln($alnEasy);
    # Get the Easyrider unaligned align object
    my $easyAll = $maindb->getAlignObject("mor", "original_seq", "accepted != 2");
    # Save that object
    my $allOut = Bio::AlignIO->new(-file => ">$archivePath\/easyriderAll.fasta", -format => 'fasta');
    $allOut->write_aln($easyAll);
print "copying html\n";    
    # Copy the current html tree files
    my $currentParsimony = $settings{'HTMLOutputDir'} . "/" . $settings{'ParsimonyHTML'};
    my $currentJackknife = $settings{'HTMLOutputDir'} . "/" . $settings{'JackknifeHTML'};
    my $currentMaximum   = $settings{'HTMLOutputDir'} . "/" . $settings{'MaximumLikelihoodHTML'};
    my @archiveFiles = ($currentParsimony, $currentJackknife, $currentMaximum);
    
    foreach my $file (@archiveFiles) {
	system("cp $file $archivePath/")
    }
    # Dump the database to the archive folder?
    
    # Zip the files
    my $zipFilename = getUniqueFile("MorFiles");
    my $zipProg = 'zip';
    
    # go into new dir
    chdir "$archivePath" or $success = 0;

    if ($success) {
	open (ZIP, "|$zipProg -9 $zipFilename *.fasta *.newick &> /dev/null") or $zipSuccess = 0;
	close (ZIP);
	if ($zipSuccess) {
	    # make the file readable
	    chmod(0755, "$zipFilename" . ".zip");
	    $log->print("Files were successfully zipped.", 1);
	} else {
	    # ok, zip failed. strange but no big deal.
	    $log->print("Failed to zip files.", 1);
	}
    } else {
	$log->print("Aborted archiving process.", 1); 
    }
    
    return $success;
}

sub cladeHTMLformat{
	my $file = shift;
	my $table = shift;
	my $linksRef = shift;
	my @links = @{$linksRef};

	open(LOG, "<$file");
	my @lines = <LOG>;
	my @returnTree;
	my $trim;
	until ($trim =~ m/\s*\/\-\s*/){
	        $trim = shift(@lines);
		if(!@lines){
			die "Not a valid paup log fed into HTMLformat";
		}
	}
	unshift(@lines, $trim);
	until ($trim =~ m/\s*\\\-\s*/){
		$trim = pop(@lines);
		if(!@lines){
			die "The top was OK the bottom was not, What happened?";
		}
	}
	push(@lines, $trim);
	my $decoration;
	foreach my $line (@lines){
	        chomp($line);
	        if ($line =~ s/(\w+)//){
	                my $name = $maindb->getValue("SELECT species FROM ".$table." WHERE accno = \"$1\"");
			$decoration = "$line $name ($1)";
			my $accno = $1;
			foreach my $url (@links){
			    my $text = $url;
			    $text =~ s/\.html.*$//;
			    $text =~ s/^.*\///;
			    $decoration .= "<a href=\"$url$accno\">$text</a> ";
			}
	                push(@returnTree, $decoration);
	        }else{
	                push(@returnTree, $line);
	        }
	}
	return @returnTree;
}


1;


__END__

=head1 NAME

mor::Utilities - various helper routines

=head1 SYNOPSIS

	use mor::Utilities;
	use mor::Log;
	my $log = mor::Log->new("log/file/name");
	$log->print("any message you want to record");
	
=head1 DESCRIPTION

Provides functions that do not exactly belong in one of the other modules, but
the functions are mostly involved in file creation

=head2 Methods

=over 6

=item B<getDateToday>

 Title   : getDateToday
 Usage   : $date = $mor::Utilities::getDateToday($baseName);
 Function: Returns a date string in the format YYYY/MM/DD (MySQL/Genbank compatible)
 Returns : A string representation of the date.
 Args    : None.


=item B<getUniqueFile>

 Title   : getUniqueFile
 Usage   : $uniqueFilename = $mor::Utilities::getUniqueFile($baseName);
 Function: Creates a unique file name by using the date and incrementing a
	   counter (if needed).  It does not simply append the date to the end
	   of the base name
 Returns : A string of a filename
 Args    : The base from which the function will create the unique filename.

=item B<lockProcessOrDie>

 Title   : lockProcessOrDie
 Usage   : $mor::Utilities::lockProcessOrDie($lockName);
 Function: Creates a concurrency file for a process that prevents other instances
	   of the program from running.  If it detects another process that is running,
	   it will die.
 Returns : Nothing
 Args    : The name of the lock in a string

=item B<clearProcessLocks>

 Title   : clearProcessLocks
 Usage   : $mor::Utilities::clearProcessLocks($lockName);
 Function: Clears a concurrency file for a process that prevents other instances
	   of the program from running.
 Returns : Nothing
 Args    : The name of the lock in a string


=item B<HTMLformat>

sub HTMLformat{
	my $htmlName = shift;
	my $treeDBID = shift;
	my $pageTitle = shift;
	my $file = shift;
	my @links = shift;
	my $table = shift;

 Title   : HTMLformat
 Usage   : my $log = mor::Log->new("log/file/name")
 Function: Creates a new mor::Log object
 Returns : mor::Log object that uses "log/file/name" as the data 
	   storage location
 Args    : A string in the form of a filepath

=item B<_introStatsHTML>

 Title   : _introStatsHTML
 Usage   : Internal function, should not be called.
 Function: Outputs a string containing HTML code to display introductory
	   statistics, e.g., the process id, number of species, etc.
 Returns : String of html for the beginning of html tree for the trees
 Args    : What type of tree it is in string form and the ID of it in the database (MLTree, NJTree, PARTree and backbone)

=item B<archiveForWeb>
                                                                                                                 
 Title   : archiveForWeb 
 Usage   : mor::Archive::archiveForWeb
 Function: Saves information from the database to a folder accessible from the website, it also
	   saves older html trees.
 Returns : Nothing
 Args    : Nothing

=back

=head1 AUTHOR

The I<mor> Team - L<http://mor.clarku.edu/morTeam.php>

=cut
