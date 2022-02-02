#! /usr/bin/perl
#
SUBDIRS += CommonUtils/TauPolSoftware/TauDecaysInterface/
system(sprintf("echo \"sed -i 's;SUBDIRS += TauSpiner/ CommonUtils/TauPolSoftware/TauDecaysInterface/;SUBDIRS += CommonUtils/TauPolSoftware/TauDecaysInterface/;g' Makefile\""));
