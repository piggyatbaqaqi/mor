package mor::DatabaseIO;

use DBI;

use mor::Settings;
use mor::Utilities;
use mor::Log;
use Bio::TreeIO;
use Bio::AlignIO;
use Bio::Seq::RichSeq;
use Bio::SeqIO;

# The global log object and settings object are created dynamically at the start
# of the program and assumed to be valid by the time this module is used. I.e.,
# this module is not valid outside of the context of the full program.
use vars qw($log %settings $maindb);

*log = \$mor::Log::log;
*settings = \%mor::Settings::settings;

sub new {
	my $class = shift;
	my $self = {};
	bless($self, $class);

	$self->{data_source} = shift;
	$self->{username} = shift;
	$self->{password} = shift;
	return $self;
}

sub LoadMainDB() {
	my $username = $settings{'username'};
	my $password = $settings{'password'};
	my $ds = $settings{'data_source'};

	$maindb = new mor::DatabaseIO ($ds, $username, $password);
	$maindb->connect();
}                    

sub connect {
	my $self = shift;
	if (! ($self->{dbh} = DBI->connect($self->{data_source}, $self->{username}, $self->{password})) ) {
		mor::Utilities::showError("Could not connect to MySQL database: [$!]", 1);
	}
	$self->{dbh}->{auto_reconnect} = 1;
}

sub query {
	my $self = shift;
	my $query = shift;

	$self->{dbh}->ping();

	unless ($self->{dbh} and $self->{dbh}->ping()) { 
		$self->connect();
	}

	$queryObj = $self->{dbh}->prepare($query);
	my $result = $queryObj->execute();
	if (!$result) {
	$log->print("Got error trying to do: $query", 1);
	}
	return $queryObj;
}

sub do{
	my $self = shift;
	my $query = shift;

	$self->{dbh}->ping();
	
	unless ($self->{dbh} and $self->{dbh}->ping()) { 
		$self->connect();
	}
	
	my $success = $self->{dbh}->do($query);
	return $success;
}

sub isInDatabase {

    my $self = shift;
    my $seqName = shift;
    my $table = shift;

    $self->{dbh}->ping();
    
    unless ($self->{dbh} and $self->{dbh}->ping()) { 
	$self->connect();
    }
    my $query = "SELECT COUNT(accno) FROM $table WHERE accno = \"$seqName\"";
    
    my $queryObj = $self->{dbh}->prepare($query);
    my $result = $queryObj->execute();
    if (!$result) {
	return 0;
    }
    if($result){
	return 1;
    }
}

sub removeFromDatabase {
    my $self = shift;
    my $seqName = shift;
    my $table = shift;

    $self->{dbh}->ping();
    
    unless ($self->{dbh} and $self->{dbh}->ping()) { 
	$self->connect();
    }
    my $query = "DELETE FROM $table WHERE accno = \"$seqName\"";
    $queryObj = $self->{dbh}->prepare($query);
    $queryObj->execute();
}

sub getValue {
	my $self = shift;
	my $query = shift;

	my $val = $self->getArray($query);
	my @retAr = @{$val};
	my $ret = $retAr[0];
	return $ret;
}

sub getArray {
	my $self = shift;
	my $query = shift;

	if ($query eq "") { return; }

	my $dbArray = $self->query($query);
	my @val;
	while (my @info = $dbArray->fetchrow_array()){
	    push(@val, @info);
	}
	return \@val;
}

sub getTreeObj {
    my $self = shift;
    my $table = shift;
    my $name = shift;

    my $tree = $self->getValue("SELECT newick FROM $table WHERE id = \'$name\' ORDER BY DATE LIMIT 0,1");
    
    open(my $fake_fh, '<', \$tree);
    my $treeI = new Bio::TreeIO(-fh => $fake_fh, -format => 'newick');
    my $treeObj = $treeI->next_tree();
    close($fake_fh);
    
    return $treeObj;

}

sub saveTreeObj {
    my $self = shift;
    my $table = shift;
    my $name = shift;
    my $tree = shift;
    
    my $output = "tempfiles/" . mor::Utilities::getUniqueFile($name);
    my $out = Bio::TreeIO->new(-file => ">$output", -format => 'newick');
    $out->write_tree($tree);
    $self->saveTreeFile($table, $output, "newick", $name);
    
}

sub saveTreeFile {
    my $self = shift;
    my $table = shift;
    my $file = shift;
    my $format = shift;
    my $name = shift;

    if ($format ne "newick") {
	my $output = "tempfiles/". mor::Utilities::getUniqueFile($name);
	my $in  = Bio::TreeIO->new(-file => $file, -format => "$format");
	my $out = Bio::TreeIO->new(-file => ">$output", -format => 'newick');

	while ( my $obj = $in->next_tree() ) {
	    $out->write_tree($obj);
	}
	$file = $output;
	
    }
    
    open(NEWICK, $file);
    my @tree = <NEWICK>;
    close(NEWICK);
    
    my $query = "INSERT INTO $table (id,newick,date) VALUES (\"". $name ."\", \"". $tree[0] ."\", CURRENT_DATE)";
    my $results = $self->do($query);
    
}

sub getAlignObject {
    my $self = shift;
    my $table = shift;
    my $seqType = shift;
    my $conditions = shift;
    
    open(FASTA, ">tempfiles/temp.fasta");
    my $query = "SELECT accno, ".$seqType." FROM ". $table ." WHERE ". $conditions ." ORDER BY date ASC";
    my $sth = $self->query($query);
    my $accno = "";
    my $seq = "";
    while ( ($accno, $seq) = $sth->fetchrow_array()){
	print FASTA ">$accno\n$seq\n";
    }
    close(FASTA);
    
    my $alnIn = new Bio::AlignIO(-file => "<tempfiles/temp.fasta", -format => "fasta");
    my $alnObj = $alnIn->next_aln();
    return $alnObj;

}

sub getEasyRider {
    my $self = shift;
    my $table = shift;
    my $seqType = shift;
    my $conditions = shift;

    my $accnoRef = $self->getArray("SELECT accno FROM ". $table ." WHERE ". $conditions ." ORDER BY date ASC");
    my @accnos = @{$accnoRef};
    my $seqRef = $self->getArray("SELECT ". $seqType ." FROM ". $table ." WHERE ". $conditions ." ORDER BY date ASC");
    my @sequences = @{$seqRef};
    my $nameRef = $self->getArray("SELECT species FROM ". $table ." WHERE ". $conditions ." ORDER BY date ASC");
    my @names = @{$seqRef};

    my $temp = "tempfiles/easyrider.fasta";
    open(FASTA, ">$temp");
    for( my $i = 0; $i < scalar(@accnos); $i++) {
	my $name = $names[$i];
	$name = s/\s/\_/g;
	print FASTA ">$accnos[$i]_$name\n$sequences[$i]\n";
    }
    close(FASTA);
    
    my $aln = new Bio::AlignIO(-file => $temp, -format => "fasta");
    my $alnOut = new Bio::AlignIO(-file => ">testing.fasta", -format => "fasta");
    my $alnObj = $aln->next_aln();
    $alnOut->write_aln($alnObj);

    return $alnObj;

}


1;

__END__

=head1 NAME

mor::DatabaseIO - An object that can be used to access a database.

=head1 SYNOPSIS

	use mor::Settings;
	use mor::Utilities;
	use mor::Log;
	use Bio::TreeIO;
	use Bio::AlignIO;
	use Bio::SeqIO;

	my $DBIO = mor::DatabaseIO->new($databaseName, $userName, $password);
	$DBIO->connect();
	my $results = $DBIO->query($MySQLquery);

=head1 DESCRIPTION

=over 3

=item B<Precondition>

 This module requires the settings hash and the log file to exist and be
 initialized.

=back

This is pretty much a wrapper for DBI. It saves lots of code in other modules
and allows for more specialized code relating to I<mor>.

=head2 Methods

=over 6

=item B<new>

 Title   : new
 Usage   : mor::DatabaseIO->new($database,$username,$password)
 Function: Creates a new mor::DatabaseIO object.
 Returns : mor::DatabaseIO object
 Args    : Database name
	   User name
	   Password

=item B<LoadMainDB>

 Title   : LoadMainDB
 Usage   : Internal Method should not be used
 Function: Loads the $maindb variable as a global variable.
	   It is called by many other modules
 Returns : nothing
 Args    : none

=item B<connect>

 Title   : connect
 Usage   : $mor::DatabaseIO->connect()
 Function: Connects the mor::DatabaseIO object to the database
 Returns : nothing 
 Args    : nothing

=item B<query>

 Title   : query
 Usage   : $retVal = $mor::DatabaseIO->query($queryString)
	   @resutls = $retVal->fetchRowArray()
 Function: This executes a query that will have a return table 
	   (ie: SELECT)
 Returns : A query object which you can call fetchRowArray on
 Args    : A MySQL query string

=item B<do>

 Title   : do
 Usage   : $mor::DatabaseIO->do($query)
 Function: Executes a query with no return table (ie: UPDATE)
 Returns : 1 if successful and 0 if unsuccessful
 Args    : A MySQL query string

=item B<isInDatabase>

 Title   : isInDatabase
 Usage   : $mor::DatabaseIO->isInDatabase($table, $accessionNumber)
 Function: Tests whether or not a given accession number is in the
	   table provided.
 Returns : 1 if the accession number exists and 0 if it does not
 Args    : A database table
	   A genbank accession number

=item B<removeFromDatabase>

 Title   : removeFromDatabase
 Usage   : $mor::DatabaseIO->removeFromDatabase($table, $accessionNumber)
 Function: Deletes an entry for any given accession number from
	   the given table
 Returns : nothing
 Args    : A database table
	   A genbank accession number

=item B<getValue>

 Title   : getValue
 Usage   : $mor::DatabaseIO->getValue($query)
 Function: Executes a query with returning table then returns
	   the first result
 Returns : The first result from a MySQL query
 Args    : A MySQL query string

=item B<getArray>

 Title   : getArray
 Usage   : $mor::DatabaseIO->getArray($query)
 Function: Executes a query and stores each row in an array. It then 
  	   adds each array to an array
 Returns : An array of arrays where the internal arrays are rows from the
	   resulting table of a MySQL query
 Args    : A MySQL query string

=item B<getTreeObj>

 Title   : getTreeObj
 Usage   : $mor::DatabaseIO->getTreeObj($table, $treeName)
 Function: This grabs the latest newick format tree from the database and puts it
	   into a Bio::Tree object
 Returns : Bio::Tree object
 Args    : A database table name
	   A tree name in the table

=item B<saveTreeObj>

 Title   : saveTreeObj
 Usage   : $mor::DatabaseIO->saveTreeObj($table, $treeName, $Bio::TreeOBJ) 
 Function: Stores the tree object into the database
 Returns : nothing
 Args    : A database table name
	   A tree name
	   A Bio::Tree object

=item B<saveTreeFile>

 Title   : saveTreeFile
 Usage   : $mor::DatabaseIO->saveTreeObj($table, $treeFilename, $treeFormat, $treeName) 
 Function: Changes a tree into newick format to be stored in the database 
 Returns : nothing
 Args    : A database table name
	   A file name
	   A format type
	   A tree name 

=item B<getAlignObject>

 Title   : getAlignObject
 Usage   : $Bio::Align = $mor::DatabaseIO->getAlignObject($table, $sequenceType, $conditions)
 Function: Gets all the sequences in the database and puts them into an 
	   alignment object
 Returns : $Bio::Align object
 Args    : A database table name
	   The type of sequence from the database
	   Conditions (what follows WHERE in MySQL)

=item B<getEasyRider>

 Title   : getEasyRider
 Usage   : $Bio::Align = $mor::DatabaseIO->getEasyRider($table, $sequenceType, $conditions)
 Function: Gets all the sequences in the database and puts them into an 
	   alignment object
 Returns : $Bio::Align object
 Args    : A database table name
	   The type of sequence from the database
	   Conditions (what follows WHERE in MySQL)

=back

=head1 AUTHOR

The I<mor> Team - L<http://mor.clarku.edu/morTeam.php>

=cut
