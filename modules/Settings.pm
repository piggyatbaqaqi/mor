package mor::Settings;

use strict;
use mor::Log;
use mor::Utilities;

use vars qw($log %settings);

# The global log object and settings object are created dynamically at the start
# of the program and assumed to be valid by the time this module is used. I.e.,
# this module is not valid outside of the context of the full program.
*log = \$mor::Log::log;

sub LoadSettings {
	my $settingsFile = shift;

	if ($settingsFile eq "") {
		$settingsFile = "settings/settings.dat";
	}

	if (! (-e $settingsFile)) {
		$log->print("Could not read settings file ($settingsFile)! Aborting program.", 1);
		exit;
	}


	open(SETTINGS, $settingsFile) || die "Can't open settings file.";
	my @settings_file = <SETTINGS>;
	close(SETTINGS);

	my $line;
	foreach $line (@settings_file) {
		$line =~ s/(^\s+)|(\s+$)//g;

		if ($line =~ /^#/ || $line =~ /^\s*$/) {
			next;
		}

		(my $varName, my $varVal) = split(/=/, $line, 2);
		$settings{"$varName"} = $varVal;
	}

}

sub SetSettings{
	my $setting = shift;
	my $value = shift;
	$setting =~ s/\s//g;
	$value =~ s/\s//g;
	my $settingsFile = "settings/settings.dat.new";
	open(SET, "<$settingsFile") or die "problem with SetSettings in Settings.pm: $!";
	
	my @settings = <SET>;
	close SET;

	foreach my $option (@settings){
		if($option =~ m/^$setting\=/){
			$option = "$setting\=$value\n";
		}
		else{next;}
	}
	open(SET, ">$settingsFile") or die "another error in SetSettings: $!";
	foreach my $line (@settings){
		print SET $line;
	}
}

1;



__END__


=head1 NAME
                                                                                                                    
mor::Settings - provides access to globally-accessible settings
                                                                                                                    
=head1 SYNOPSIS
                                                                                                                    
	use strict;
	use DBI;
	use mor::Log;
	use mor::Utilities;
	use vars qw($log %settings $dbmain);                                                                                                                    

=head1 DESCRIPTION
                                                                                                                    
This module provides a basis for reading in and storing settings in a globally-accessible hash.
                                                                                                                    
=head2 Methods
                                                                                                                    
=over 6
                                                                                                                    
=item B<LoadSettings>
                                                                                                                 
 Title   : LoadSettings
 Usage   : mor::Settings::LoadSettings($filename)
 Function: Loads all the settings for mor to run
 Returns : Nothing
 Args    : the location and name of the settings file

=item B<SetSettings>
                                                                                                                 
 Title   : SetSettings
 Usage   : mor::Settings::LoadSettings($key, $value)
 Function: Creates a new key with value for mor to use
 Returns : Nothing
 Args    : the new key for mor
	   the new value for that key

=back

=head1 AUTHOR

The I<mor> Team - L<http://mor.clarku.edu/morTeam.php>

=cut

