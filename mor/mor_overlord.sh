#!/bin/bash
#simple script to call all of mor's scripts in order without wasting time
#in the crontab
perl mor_backup.pl
perl mor_download.pl
perl mor_analysis.pl
