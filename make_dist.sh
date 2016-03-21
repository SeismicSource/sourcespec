#!/bin/bash
export LC_ALL=C

git_tree=HEAD
branch=`git branch | grep \* | awk '{print $2}'`
version=`git tag | tail -n1`
timestamp=`date +"%Y%m%d"`
suffix=$timestamp
if [ $# -ge 1 ]
then
    if [ $1 == "release" ]
    then
        if [ $branch != "master" ]
        then
            echo "Error: you have to be in master to make a release!"
            exit 1
        fi
        git_tree=$version
        suffix=$version
    fi
fi

if [ $branch == "master" ]
then
    dir=source_spec_${suffix}
else
    dir=source_spec_${branch}_${suffix}
fi
tarfile=${dir}.tgz

git archive --prefix=$dir/ $git_tree | tar xv
python source_spec/version.py
cp source_spec/RELEASE-VERSION $dir/source_spec

cd $dir/doc &&
    make html-install &&
    make latexpdf-install &&
    make clean &&
cd - || exit
find $dir -iname "*.pyc" | xargs rm

tar cvfz ../$tarfile $dir
rm -rf $dir
