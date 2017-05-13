package mor::Clades;

use strict;

use mor::Analysis;
use mor::Settings;
use mor::DatabaseIO;
use Bio::TreeIO;
use Bio::SeqIO;
use vars qw(*settings $maindb );

# The global log object and settings object are created dynamically at the start
# of the program and assumed to be valid by the time this module is used. I.e.,
# this module is not valid outside of the context of the full program.
*settings = \%mor::Settings::settings;
*maindb = \$mor::DatabaseIO::maindb;

sub _getCladeConstraint {
        my $cladeName = shift;
        my $table = shift;
	my $tree = shift;
	my $results = $maindb->query("SELECT const FROM $table WHERE cladename = \"$cladeName\"");
        my @row = $results->fetchrow_array();
	my @const = split(/:/,$row[0]);
	my @constraints;
	my @nodes;
	foreach my $a (@const){
		eval{
			 @nodes = $tree->find_node(-id=>$a);
		};
		if($@){
			print $@."\n";
			next;
		}
		@nodes = $tree->find_node(-id=>$a);
		if($nodes[0]){
			push (@constraints, $nodes[0]);
		}
	}
	return \@constraints;
}

sub _getSubtree {
        my $constraintsref = shift;
        my $tree = shift;
	my $clade = shift;
	my $singleClade = shift;
	my @constraints = @{$constraintsref};
        if(scalar(@constraints) >= 2){
                my $lca = $tree->get_lca(-nodes=>\@constraints);
		my $subtree = Bio::Tree::Tree->new(-root => $lca);
		my @newnodes = $subtree->get_nodes();
		if(!$singleClade){
			foreach my $anewnode (@newnodes){
					if($anewnode->is_Leaf()){
						$maindb->do("UPDATE mor SET clade = \"$clade\" WHERE accno = \"".$anewnode->id."\"");
					}
			}
		}
		return $subtree;
        }
        else{
                return 0;
        }
}

sub updateClades {
	my $table = shift;
        my $treeObj = shift;
	my $query = shift;
	my $singleClade = shift;
	if(!$query){
		$query = "SELECT cladename FROM $table WHERE public=1";
	}
	if(!$treeObj){
		$treeObj=_getTree("PARTree", "trees");
	}
        my $results = $maindb->query($query);
	if($singleClade){
		my @row=$results->fetchrow_array();
                my $constraintsref = _getCladeConstraint($row[0],$table, $treeObj);
                my $cladeObj = _getSubtree($constraintsref, $treeObj, $row[0],$singleClade);
		if($cladeObj){
        		my $out = Bio::TreeIO->new(-file   => '>tmp/clade.new' ,
						   -format => 'newick');
                        $out->write_tree($cladeObj);
                        open(CLADE, "<tmp/clade.new");
                        my @file = <CLADE>;
                        close CLADE;
			chop($file[0]);
                	$maindb->do("UPDATE $table SET newick = \"$file[0]\" WHERE cladename = \"$row[0]\"");	
                }
	} else {
        	while (my @row = $results->fetchrow_array()){
			my $constraintsref = _getCladeConstraint($row[0],$table,$treeObj);
       		        my $cladeObj = _getSubtree($constraintsref, $treeObj, $row[0],$singleClade);
			if($cladeObj){
        			my $out = Bio::TreeIO->new(-file   => '>tmp/clade.new' ,
							   -format => 'newick');
                	        $out->write_tree($cladeObj);
                	        open(CLADE, "<tmp/clade.new");
                	        my @file = <CLADE>;
                	        close CLADE;
				chop($file[0]);
                	        $maindb->do("UPDATE $table SET newick = \"$file[0]\" WHERE cladename = \"$row[0]\"");	
                	}
        	}
	}
}

sub _getTree {
	my $treeName = shift;
	my $table = shift;
	my $res = $maindb->query("SELECT newick FROM $table WHERE id=\"$treeName\"");
	my @row = $res->fetchrow_array();
	open(CLADE,">tmp/clade.new");
	print CLADE $row[0];
	close CLADE;
	my $in = Bio::TreeIO->new(-file => 'tmp/clade.new',
				-format => 'newick');
	my $tree = $in->next_tree();
	return $tree;
}

sub updateSingleClade {
	my $cladeName = shift;
	my $tree = _getTree("PARTree", "trees");
	updateClades("clades",$tree,"SELECT cladename FROM clades WHERE cladename=\"$cladeName\"",1);
	my $nexusfile = getNexusFile($cladeName);
	my $log = drawTree($tree,$nexusfile,$cladeName,'/home/www/treegraphs/logTree.nxs');
	return $log;
}

sub getNexusFile{
	my $cladename = shift;
	my $res = $maindb->query("SELECT newick FROM clades WHERE cladename=\"$cladename\"");
	my @row = $res->fetchrow_array();
	$row[0] =~ s/[(|)|;]//g;
	my @accnos = split(/,/,$row[0]);
	my $query = "SELECT accno,original_seq FROM mor WHERE accno=\"1\"";
	foreach my $a (@accnos){
		$query.=" OR accno=\"$a\"";
	}
	$res = $maindb->query("$query");
	open(FAST, ">/home/www/treegraphs/$cladename.fasta");
	while (my @ans = $res->fetchrow_array()){
		print FAST ">$ans[0]\n";
		print FAST "$ans[1]\n\n";
	}
	close FAST;
	my $in = Bio::AlignIO->new(-file => "/home/www/treegraphs/$cladename.fasta",
				-format => 'Fasta');
	open(J, "/home/www/treegraphs/$cladename-taxa.nexus");
    chmod 0777, "/home/www/treegraphs/$cladename-taxa.nexus";
    my $out = Bio::AlignIO->new(-file=>">/home/www/treegraphs/$cladename-taxa.nexus",
				-format=>'nexus');
	while(my $obj = $in->next_aln()){
		$out->write_aln($obj);
	}
	return "/home/www/treegraphs/$cladename-taxa.nexus";
}

sub drawTree {
	my $tree = shift;
	my $taxablock = shift;
	my $cladename = shift;
	my $paupTemplate = shift;
	my $res = $maindb->query("SELECT id FROM clades WHERE cladename=\"$cladename\"");
	my @row = $res->fetchrow_array();
	my $cladeTree=_getTree($row[0],"clades");
	my $savefile = "/home/www/treegraphs/$cladename.nexus";
	my $outfile = "/home/www/treegraphs/$cladename.paup";
	my $bio = Bio::TreeIO->new(-file => ">$savefile",
			 	 -format => 'nexus');
	$bio->write_tree($cladeTree);
	editNexusFile($savefile);
	my @keywords = ("REPLACELOG","REPLACETAXA","REPLACETREE");
	my @replace = ("/home/www/treegraphs/$cladename.log",$taxablock,"/home/www/treegraphs/$cladename.nexus");
	open(LOG, ">".$replace[0]);
	chmod 0777, $replace[0];
	close LOG;
	my $run = mor::Analysis::createPAUPfile($taxablock,$outfile,$paupTemplate,@keywords,@replace);
	mor::Analysis::runPAUP($run);
	return $replace[0];	
}

sub editNexusFile{
	my $file = shift;
	open(NEX, "<$file");
	my @lines = <NEX>;
	close NEX;
	for(my $i=0;$i<scalar(@lines);$i++){
		if($lines[$i]=~ /tree Bioperl/){
			$lines[$i]=~s/\(//;
			$lines[$i]=~s/\[\]\)//;
		}
	}
	open(F, ">$file");
	foreach my $e (@lines){
		print F $e;
	}
	close F;
	return;
}

1;

__END__

=head1 NAME

mor::Clades - A module for the I<mor> pacakge.

=head1 SYNOPSIS

	use mor::Settings;
	use mor::DatabaseIO;
	use Bio::TreeIO;
	use vars qw(*settings $maindb );
	
	my $cladeTree = $treeIO->next_tree();
	$MySQLtable = "table";
	mor::Clades::updateClades($MySQLtable, $cladeTree);

=head1 DESCRIPTION

=over 3

=item B<Precondition>

 This module requires the settings hash and maindb object to be creaeted
 before it is called.

=back

This package should be called with mor::Clades::updateClades. It will
take the constraints table and rip apart the tree passed in to create
as many different clades as are in the table.

=head2 Methods

=over 6

=item B<_getCladeConstraint>

 Title   : _getCladeConstraint
 Usage   : Internal method should not be used
 Function: Grabs a list of accession numbers that are constraints to the
	   given clade name. 
 Returns : Array reference to the list of accession numbers
 Args    : Clade name 
	   MySQL table

=item B<_getNodeIdHash>

 Title   : _getNodeIdHash
 Usage   : Internal method should not be used
 Function: Takes a Bio::Tree object and returns a hash of accession
	   numbers where the key is the accession number and the value
	   stored is 1
 Returns : Hash reference to the list of accession numbers
 Args    : Bio::Tree object

=item B<_getSubtree>

 Title   : _getSubtree
 Usage   : Internal method should not be used
 Function: Creates a Bio::Tree object that represents the clade 
 Returns : Either the Bio::Tree clade or 0 if the method failed
 Args    : Array reference of the constraints
	   Bio::Tree object that has the main tree
	   Hash reference to the exisiting accession numbers

=item B<updateClades>

 Title   : updateClades
 Usage   : mor::Clades::updateClades($MySQLtable, $Bio::TreeOBJ)
 Function: Updates the MySQL table to have newick trees that represent
	   the clade name.
 Returns : nothing
 Args    : the name of a MySQL table
	   Bio::Tree object with the main tree data

=item B<_getTree>

 Title   : _getTree
 Usage   : Internal method should not be used
 Function: Gets a tree from the database
 Returns : Bio::Tree object of the current tree stored in trees
 Args    : a string that is the name of the tree i.e. "PARTree"

=item B<updateSingleClade>

 Title   : updateSingleClade
 Usage   : mor::Clades::updateClades("cladeName")
 Function: This will update the database with the new clade. Or it will
	   create a clade if it has not been run before.
 Returns : nothing
 Args    : A string that is the name of the clade to be updated

=back

=head1 AUTHOR

The I<mor> Team - L<http://mor.clarku.edu/morTeam.php>

=cut
