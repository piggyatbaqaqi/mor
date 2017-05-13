package mor::Archive;

use strict;

use Bio::TreeIO;
use Bio::AlignIO;
use mor::Settings;
use mor::Utilities;
use mor::DatabaseIO;

# The global log object and settings object are created dynamically at the start
# of the program and assumed to be valid by the time this module is used. I.e.,
# this module is not valid outside of the context of the full program.
use vars qw( $log *settings $maindb);
*log = \$mor::Log::log;
*settings = \%mor::Settings::settings;
*maindb = \%mor::DatabaseIO::maindb;


sub archiveForWeb {
	my $success = 1;
	my $zipSuccess = 1;
	my $gmtime=gmtime;
	my @gmarray = split(' ', $gmtime);

	my $webPath = $settings{'webPath'};
	print "$gmarray[2]\n$gmarray[1]\n$gmarray[0]\n";
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

	my $archivePath = getNewFilenameVersion("$webPath/archive/$dirname");

	# make unique, "time-stamp-named" dir
	mkdir "$archivePath" or $success = 0; 
	chmod(0755,"$archivePath") or $success = 0;

	if (!$success) { 
		showError("Unable to create archive path [$archivePath]. Aborting archival process."); 
		return; 
	}
	print "ARCHIVE PATH IS: $archivePath";
	# Get the ML tree object
	my $mlTree = $maindb->getTreeObj("trees", "MLTree");
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
	$parOut->write_tree($njTree);
	
	# Get the Easyrider aligned align object
	my $alnEasy = $maindb->getAlignObject("mor", "aligned_seq", "accepted = 1 AND HRI != 0");
	# Save that object
	my $easyOut = Bio::AlighIO->new(-file => ">$archivePath\/easyriderAligned.fasta", -format => 'fasta');
	$easyOut->write_tree($alnEasy);
	# Get the Easyrider unaligned align object
	my $easyAll = $maindb->getAlignObject("mor", "original_seq", "accepted != 2");
	# Save that object
	my $allOut = Bio::AlighIO->new(-file => ">$archivePath\/easyriderAll.fasta", -format => 'fasta');
	$allOut->write_tree($easyAll);
	

	# Copy the current html tree files
	my $currentParsimony = $settings{'HTMLOutputDir'} . "/" . $settings{'ParsimonyHTML'};
	my $currentJackknife = $settings{'HTMLOutputDir'} . "/" . $settings{'JackknifeHTML'};
	my $currentMaximum   = $settings{'HTMLOutputDir'} . "/" . $settings{'MaximumLikelihoodHTML'};
	my @archiveFiles = ($currentParsimony, $currentJackknife, $currentMaximum);

	foreach my $file (@archiveFiles) {
		copy($file, "$archivePath/");
	}

	# Dump the database to the archive folder?


	# Zip the files
	my $zipFilename = mor::Utilites::getUniqueFile("MorFiles");
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

1;


__END__

=head1 NAME
                                                                                                                    
mor::Archive - Archives the useful information from the database to the website
                                                                                                                    
=head1 SYNOPSIS
                                                                                                                    
	use Bio::TreeIO;
	use Bio::AlignIO;
	use mor::Settings;
	use mor::Utilities;
	use mor::DatabaseIO;
	use mor::Archive;
                                                                                                                    
=head1 DESCRIPTION
                                                                                                                    
This module saves information from the database to the website
                                                                                                                    
=head2 Methods
                                                                                                                    
=over 6
                                                                                                                    
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
