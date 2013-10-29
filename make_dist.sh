#!/bin/bash
export LC_ALL=C

git_tree=HEAD
branch=`git branch | grep \* | awk '{print $2}'`
version=`git tag | tail -n1`
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
    fi
fi

timestamp=`date +"%Y%m%d"`

if [ $branch == "master" ]
then
    dir=source_spec_${timestamp}
else
    dir=source_spec_${branch}_${timestamp}
fi
tarfile=${dir}.tgz

git archive --prefix=$dir/ $git_tree | tar xv

cd $dir/doc_src &&
    make html-install &&
    make latexpdf-install &&
    make clean &&
cd -
find $dir -iname "*.pyc" | xargs rm

tar cvfz ../$tarfile $dir
rm -rf $dir
