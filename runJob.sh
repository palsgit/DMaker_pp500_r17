#!/bin/bash
cd workDir
productionId=`date +%F_%H-%M`_D0_DTlusty_thesis_magiccut
analyzer="palsp"

mkdir $productionId
cd $productionId
#copylist
cp ../../picoLists/part1.list  ./
list="part1.list"

#copy needed folders
cp -r ../../.sl73_gcc485 ./
cp -Lr ../../StRoot ./
cp ../../picoLists/picoList_bad.list ./
mkdir starSubmit
cp ../../starSubmit/submitPicoHFMaker.xml ./starSubmit

mkdir -p production
mkdir -p report
mkdir -p csh
mkdir -p list
mkdir -p jobs
mkdir -p jobs/log
mkdir -p jobs/err

path=`pwd -P`
path=$( echo $path | sed 's|//|/|g' )

baseFolder=${path}
rootMacro=runPicoD0AnaMaker.C
inputList=part1.list

star-submit-template -template ./starSubmit/submitPicoHFMaker.xml -entities listOfFiles=${inputList},basePath=${baseFolder},prodId=${productionId},rootMacro=${rootMacro}
