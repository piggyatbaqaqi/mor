package mor::Log;

use strict;

use mor::Utilities;
use vars qw($log);

sub new {
	my $class = shift;
	my $logName = shift;
	my $uniqueLog = mor::Utilities::getUniqueFile($logName);
	my $self = {filename=>"$uniqueLog"};
	open(LOG, ">$uniqueLog");
	
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $time = "$hour:$minute:$second";
	$month++;
	my $year = $yearOffset+1900;
	print LOG "Starting mor at $time on $month\/$dayOfMonth\/$year\n";
	close LOG;
	bless($self, $class);

	return $self;
}


sub print {
	my $self = shift;
	my $message = shift;
	my $priority = shift;

	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	$yearOffset += 1900;
	$month++;
	my $time = "$month-$dayOfMonth-$yearOffset $hour:$minute:$second";
	print "$time \=\> $message\n";
	open(LOG_HANDLE, ">>" . $self->{filename});
	print LOG_HANDLE "$time \=\>  $message\n";
	close(LOG_HANDLE);

}

1;

__END__

=head1 NAME

mor::Log - Log handler

=head1 SYNOPSIS

	use mor::Utilities;
	use mor::Log;
	my $log = mor::Log->new("log/file/name");
	$log->print("any message you want to record");
	
=head1 DESCRIPTION

A very simple object that takes a file name, in the form of a string, 
as its only parameter. It uses the getUniqueFileName found in 
mor::Utilites to always create a new log.

=head2 Methods

=over 6

=item B<new>

 Title   : new
 Usage   : my $log = mor::Log->new("log/file/name")
 Function: Creates a new mor::Log object
 Returns : mor::Log object that uses "log/file/name" as the data 
	   storage location
 Args    : A string in the form of a filepath

=item B<print>

 Title   : print
 Usage   : $log->print("message")
 Function: Prints "message" to both the screen and the file 
	   associated with the 
 Returns : nothing
 Args    : Any string

=back

=head1 AUTHOR

The I<mor> Team - L<http://mor.clarku.edu/morTeam.php>

=cut
