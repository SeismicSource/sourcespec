#!/bin/bash
export LC_ALL=C

timestamp=`date +"%Y%m%d"`

tarfile=source_spec_${timestamp}.tgz

tar cvfz ../$tarfile *.py README.txt testdata
