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

mor::Utilities::lockProcessOrDie("html");

# Connect the database
mor::DatabaseIO::LoadMainDB();

my $lsu = $settings{'LSUtable'};

my @mlReplace;
my @mllinks = ("http://mor.clarku.edu/nj.html#","http://mor.clarku.edu/ml.html#");
my @mlTree = mor::Utilities::HTMLformat("Maximum Likelihood HTML Tree", "MLTree", "mor - Jackknife NJ Tree", "paupfiles/ml_humanRun.log", \@mllinks, $lsu);
open(MLTREE, ">htmlfiles/ml.html") or die;
foreach my $line (@mlTree){
	print MLTREE "$line\n";
}
close MLTREE;

my @NJlinks = ("http://mor.clarku.edu/mp.html#","http://mor.clarku.edu/ml.html#");
my @njTree = mor::Utilities::HTMLformat("JackknifeHTML", "NJTree", "mor - Jackknife NJ Tree", "paupfiles/njtree-2008-10-01.log", \@NJlinks, $lsu);
open(NJTREE, ">htmlfiles/nj.html") or die;
foreach my $line (@njTree){
	print NJTREE "$line\n";
}
close NJTREE;

my @PARlinks = ("http://mor.clarku.edu/nj.html#","http://mor.clarku.edu/ml.html#");
my @parTree = mor::Utilities::HTMLformat("ParsimonyHTML", "PARTree", "mor - Consensus MP Tree", "paupfiles/parsimony-2008-10-05.log", \@PARlinks, $lsu);
open(PARTREE, ">htmlfiles/mp.html") or die;
foreach my $line (@parTree){
	print PARTREE "$line\n";
}
close PARTREE;

`cp htmlfiles/* /home/www/`;

mor::Utilities::clearProcessLocks("html");
