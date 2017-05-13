package mor::Analysis;

use strict;
use mor::Log;
use mor::Settings;
use mor::DatabaseIO;
use mor::Utilities;
use Bio::TreeIO;
use vars qw( $log *settings $maindb );

# The global log object and settings object are created dynamically at the start
# of the program and assumed to be valid by the time this module is used. I.e.,
# this module is not valid outside of the context of the full program.
*log = \$mor::Log::log;
*settings = \%mor::Settings::settings;
*maindb = \$mor::DatabaseIO::maindb;


sub raxml {

    my $alnObj = shift;

    my $newPHYLIP = mor::Utilities::getUniqueFile($settings{'raxmlPHY'},1);
    
    $log->print("Writing the alignment to a phylip file for raxml.");
    my $outPHY = Bio::AlignIO->new(-file => ">$newPHYLIP", -format => 'phylip');
    $outPHY->write_aln($alnObj);

    
    my $raxmlProg = $settings{'raxml'};
    my $raxmlParams = $settings{'raxmlParams'};
    my $raxmlFile = $settings{'raxmlFile'};
    my $raxmlLog = mor::Utilities::getUniqueFile("$settings{'raxmlLog'}");
    my $raxmlResult = mor::Utilities::getUniqueFile("$settings{'raxmlTree'}");
    my $raxmlInfo = mor::Utilities::getUniqueFile("$settings{'raxmlInfo'}");
    my $newickConst = $settings{'newickConst'};

    # Running RAxML
    $log->print("Starting RAxML process... sit back and relax.", 1);

    system ( "$raxmlProg -s $newPHYLIP -g $newickConst -n $raxmlFile $raxmlParams");
    
    # Moving a couple outfiles from raxml into an appropriate location
    system ( "mv RAxML_log.$raxmlFile $raxmlLog");
    system ( "mv RAxML_result.$raxmlFile $raxmlResult");
    system ( "mv RAxML_info.$raxmlFile $raxmlInfo");
    
    $log->print("RAxML analysis finished.", 1);	

    my @raxmlFiles = ($raxmlResult, $raxmlInfo, $raxmlLog);
    
    return @raxmlFiles;

}


sub createPAUPfile ($$$\@\@) {
    my $alignment = shift;
    my $outfile = shift;
    my $paupTemplate = shift;
    my (@keywords) = @{(shift)};
    my (@replacements) = @{(shift)};

#read in the pauptemplate into an array
    open(FILE, "$paupTemplate");
    my @paup = <FILE>;
    close(FILE);

#replace the keywords with the replacements
    for(my $i = 0; $i < scalar(@keywords); $i++){
	for( my $j = 0; $j < scalar(@paup); $j++) {
	    $paup[$j] =~ s/$keywords[$i]/$replacements[$i]/g;
	}
    }

#append the updated information onto the alignment
    open(ALIGNMENT, "$alignment");
    my @alignment = <ALIGNMENT>;
    close(ALIGNMENT);

    open(PAUP, ">$outfile");
    foreach my $line (@alignment) {
	print PAUP $line;
    }
    foreach my $line (@paup) {
	print PAUP $line;
    }
    close(ALIGNMENT);
    
    return $outfile;

}

sub runPAUP {
    # This is rather simplistic but paup might need to be wrapped up better later
    my $file = shift;
    my $cmd = "paup -n $file ";
    system($cmd);
}

# a mafft function would fit in well here now

1;


__END__

=head1 NAME
                                                                                                                    
mor::Analysis - Phylogenetic analysis functions
                                                                                                                    
=head1 SYNOPSIS

        use strict;
        use mor::Log;
        use mor::Settings;
        use mor::DatabaseIO;
        use mor::Utilities;
        use Bio::TreeIO;
        use vars qw( $log *settings $maindb );

=head1 DESCRIPTION

This module does phylogenetic analyses using the programs RAxML and PAUP*

=head2 Methods

=over 6

=item B<raxml>

 Title   : raxml
 Usage   : my @raxmlFiles = mor::Analysis::raxml($alnObj)
 Function: Runs a ML analysis on a passed in Bio::Align object
 Returns : Returns an array of the file names that are created
 Args    : a Bio::Align object.


=item B<createPAUPfile>

 Title   : createPAUPfile
 Usage   : my $runFile = mor::Analysis::createPAUPfile($nexus, $batchFile, $paupBlock, @keywords,@replacements)
 Function: Creates a PAUP* batch file, by appending a paup block from a nexus file to a character
	   matrix in a nexus file.  
 Returns : The file name of the batch file (same as $batchFile)
 Args    : The file name of an alignment in nexus format 
	   The file name of the batch file that will be run with mor::runPAUP
	   The file name of the paup command block
	   An array of words in the paup block that are to be replaced
	   An array of related words in that are filled in.

=item B<runPAUP>

 Title   : runPAUP
 Usage   : mor::Analysis::runPAUP($runFile)
 Function: Runs PAUP* on a supplied batch file
 Returns : Nothing
 Args    : The file name of the batch file
	
=back

=head1 AUTHOR

The I<mor> Team - L<http://mor.clarku.edu/morTeam.php>

=cut
