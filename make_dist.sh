#!/bin/bash
export LC_ALL=C

timestamp=`date +"%Y%m%d"`

branch=`git branch | grep \* | awk '{print $2}'`
if [ $branch == "master" ]
then
    dir=source_spec_${timestamp}
else
    dir=source_spec_${branch}_${timestamp}
fi
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
