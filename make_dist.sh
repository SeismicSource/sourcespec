#!/bin/bash
export LC_ALL=C

timestamp=`date +"%Y%m%d"`

dir=source_spec_${timestamp}
tarfile=${dir}.tgz

mkdir $dir
rsync -a\
    --exclude='.*'\
    --exclude='*.pyc'\
    --exclude='*.sqlite'\
    --exclude='sspec_out'\
    --exclude='make_dist.sh'\
    --exclude="$dir"\
    . $dir
tar cvfz ../$tarfile $dir
rm -rf $dir
