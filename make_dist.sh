#!/bin/bash
export LC_ALL=C
TESTDATA=testdata

timestamp=`date +"%Y%m%d"`

tarfile=source_spec_${timestamp}.tgz

tar cvfz ../$tarfile *.py README.md ChangeLog.txt $TESTDATA
