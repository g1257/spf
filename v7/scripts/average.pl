#!/usr/bin/perl -w

use strict;
use lib 'scripts';

use Average;
my ($label)=@ARGV;

Average::average($label);

