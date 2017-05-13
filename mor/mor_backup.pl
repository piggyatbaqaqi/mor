#!/usr/bin/perl
#
# This script copies and saves all the important bits and pieces of mor...
# in case of bad things happening.

use strict;

chdir '/misc/mor/';
use vars qw($log %settings $maindb);

use mor::Acquisition;
use mor::DatabaseIO;
use mor::Log;
use mor::Settings;

*log = \$mor::Log::log;
$log = new mor::Log("logs/.log");
*maindb = \$mor::DatabaseIO::maindb;
*settings = \%mor::Settings::settings;

# Load the settings from the settings file
mor::Settings::LoadSettings("settings/settings.dat");

mor::Utilities::lockProcessOrDie("backup");

#Connect the database
mor::DatabaseIO::LoadMainDB();

#Dump all the tables in the database
system("mysqldump --user=$settings{'username'} --password=$settings{'password'} -h 127.0.0.1 mormain > $settings{'backupLoc'}\/mormain.sql");

#Copy the mor modules
system("cp -r /usr/lib/perl5/site_perl/5.8.5/mor $settings{'backupLoc'}\/mm");

#Copy the mor folder
system("cp -r * $settings{'backupLoc'}\/mor");

#tarball it up to save space
my $datedName = mor::Utilities::getUniqueFile("morBackup.tar.gz");
system("tar -czvf $datedName $settings{'backupLoc'}\/*");
system("mv $datedName $settings{'backupLoc'}\/");
system("rm -rf $settings{'backupLoc'}\/mor");
system("rm -rf $settings{'backupLoc'}\/mm");
system("rm -rf $settings{'backupLoc'}\/mormain.sql");

mor::Utilities::clearProcessLocks("backup");
