#!/usr/bin/perl
#
# This script is responsible for analyzing sequences that have been downloaded
# and aligned.
use strict;

chdir '/misc/mor/';
use vars qw($log %settings $maindb);

use mor::Analysis;
use mor::Archive;
use mor::DatabaseIO;
use mor::Log;
use mor::Settings;
use mor::Utilities;
use mor::Clades;
use Bio::AlignIO;
use Bio::Align::AlignI;
use Bio::TreeIO;
use Bio::SeqIO;

*log = \$mor::Log::log;
$log = new mor::Log("logs/.log");

*maindb = \$mor::DatabaseIO::maindb;
*settings = \%mor::Settings::settings;

$log->print("Starting mor analysis",1);

# Load the settings from the settings file
mor::Settings::LoadSettings("settings/settings.dat");

mor::Utilities::lockProcessOrDie("analysis");

# Connect the database
mor::DatabaseIO::LoadMainDB();

my $lsu = $settings{'LSUtable'};

# Going to create the nexus file from the database to work with
$log->print("Getting the accepted sequences from ". $lsu,1);

my $alignment = mor::Utilities::getUniqueFile($settings{'NXSfile'});

my $alnObj = $maindb->getAlignObject($lsu, "align_seq", "accepted = 1");

my $out = Bio::AlignIO->new(-file => ">$alignment", -format => "nexus");
$out->write_aln($alnObj);

# Create the ML tree
$log->print("Creating a maximum likelihood using RAxML.",1);

my @raxmlFiles = mor::Analysis::raxml($alnObj);

$maindb->saveTreeFile("trees", $raxmlFiles[0], "newick", "MLTree");

$log->print("Creating the NJ tree using PAUP\*", 1);

# Create the paup batch file from a template
my @njKeywords = ("REPLACENJTREE", "REPLACELOG", "REPLACETREEOUT", "REPLACETREEFINAL");
my @njReplace;

push(@njReplace, mor::Utilities::getUniqueFile($settings{'NJNXS'}));
push(@njReplace, mor::Utilities::getUniqueFile($settings{'NJLog'}));
push(@njReplace, mor::Utilities::getUniqueFile($settings{'NJJKreplicates'}));
push(@njReplace, mor::Utilities::getUniqueFile($settings{'ConsensusNJTree'}));

my $njRunFile = mor::Utilities::getUniqueFile($settings{'NJNXS'});
my $runNJ = mor::Analysis::createPAUPfile($alignment, $njRunFile, $settings{'paupNJrun'},@njKeywords,@njReplace);

# Run PAUP on the created file
mor::Analysis::runPAUP($runNJ);

# Load in the tree after a decade's wait and write it
$maindb->saveTreeFile("trees", $njReplace[3], "nexus", "NJTree");

# Repeat the process one more time, only this time for the parsimony tree
$log->print("Creating the parsimony tree using PAUP\*", 1);

my $bbTree = $maindb->getTreeObj("trees", "backbone");
my $bbFile = mor::Utilities::getUniqueFile($settings{'constraints'});
my $treeOut = Bio::TreeIO->new(-file => ">$bbFile", -format => "nexus");
$treeOut->write_tree($bbTree);

my @parKeywords = ("REPLACEBACKBONE", "REPLACETIME", "REPLACEPARTREE", "REPLACESCORE", "REPLACETREEFINAL", "REPLACELOG");
my @parReplace;
push(@parReplace, "backbone.nex");
push(@parReplace, $settings{'parsimonyTime'});
push(@parReplace, mor::Utilities::getUniqueFile($settings{'ParsimonyNXS'}));
push(@parReplace, mor::Utilities::getUniqueFile($settings{'scorefile'}));
push(@parReplace, mor::Utilities::getUniqueFile($settings{'TreePars'}));
push(@parReplace, mor::Utilities::getUniqueFile($settings{'ParsimonyLog'}));

my $runPAR = mor::Analysis::createPAUPfile($alignment, mor::Utilities::getUniqueFile($settings{'ParsimonyNXS'}), $settings{'paupPARSrun'},@parKeywords,@parReplace);

# Run paup again, you can come back in a couple days now if you'd like
mor::Analysis::runPAUP($runPAR);

$maindb->saveTreeFile("trees", $parReplace[4], "nexus", "PARTree");

$log->print("Analysis finished -- saving time log.", 1);

# This next step is needed for the website, it creates a paup log with the tre
my $MLnxsFile = $settings{'mlNXSTree'};
my $mlObj = $maindb->getTreeObj("trees", "MLTree");
my $out = Bio::TreeIO->new(-file => ">$MLnxsFile", -format => 'nexus');
$out->write_tree($mlObj);

my @logKeywords = ("REPLACELOG","REPLACETREE");
my @MLreplacements = (mor::Utilities::getUniqueFile($settings{'paupMLLog'}), $MLnxsFile);

my $runML = mor::Analysis::createPAUPfile($alignment, mor::Utilities::getUniqueFile($settings{'paupMLLog'}), $settings{'mlPaupTree'}, @logKeywords, @MLreplacements);

$log->print("Running paup to print ML tree.",1);

mor::Analysis::runPAUP($runML);


# Now it's time to create all the html pages for the trees, maybe this step should be done in a separate script, update_website.pl maybe?
$log->print("Dressing the trees for the website", 1);

my @MLlinks = ("http://mor.clarku.edu/nj.html#", "http://mor.clarku.edu/mp.html#");
my @MLTree = mor::Utilities::HTMLformat("Maximum Likelihood HTML Tree", "MLTree", "mor - Maximum Likelihood Tree", $MLreplacements[1], \@MLlinks, $lsu);
open(TREE, ">trees/ml.html");
foreach my $line (@MLTree) {
    print TREE "$line\n";
}
close TREE;

my @NJlinks = ("http://mor.clarku.edu/mp.html#","http://mor.clarku.edu/ml.html#");
my @njTree = mor::Utilities::HTMLformat("Neighbor-Joining Jackknife HTML Tree", "NJTree", "mor - Jackknife NJ Tree", $njReplace[1], \@NJlinks, $lsu);
open(NJTREE, ">trees/nj.html");
foreach my $line (@njTree){
	print NJTREE "$line\n";
}
close NJTREE;

my @PARlinks = ("http://mor.clarku.edu/nj.html#","http://mor.clarku.edu/ml.html#");
my @parTree = mor::Utilities::HTMLformat("Maximum Parsimony HTML Tree", "PARTree", "mor - Consensus MP Tree", $parReplace[5], \@PARlinks, $lsu);
open(PARTREE, ">trees/mp.html");
foreach my $line (@parTree){
	print PARTREE "$line\n";
}
close PARTREE;

#The module Clades.pm no longer requires or even uses phylip in any way.
$log->print("Updating the clades for the website");
my $cladeTree = $maindb->getTreeObj("trees", "PARtree");
mor::Clades::updateClades("clades", $cladeTree);

$log->print("Archiving info");
mor::Utilities::archiveForWeb();

$log->print("mor analysis completed.");

mor::Utilities::clearProcessLocks("analysis");

