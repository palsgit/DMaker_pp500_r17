#!/bin/bash
cd workDir
productionId=`date +%F_%H-%M`_following_lastD0woHFTwVr_loose_part2
analyzer="palsp"

mkdir $productionId
cd $productionId
#copylist
##cp ../../picoLists/test.list  ./
##list="test.list"

#copy needed folders
cp -r ../../.sl73_gcc485 ./
cp -Lr ../../StRoot ./
cp ../../picoLists/picoList_bad.list ./
mkdir starSubmit
cp ../../starSubmit/submitPicoHFMaker_part2.xml ./starSubmit

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
##inputList=test.list

##star-submit-template -template ./starSubmit/submitPicoHFMaker_part2.xml -entities listOfFiles=${inputList},basePath=${baseFolder},prodId=${productionId},rootMacro=${rootMacro}
star-submit-template -template ./starSubmit/submitPicoHFMaker_part2.xml -entities basePath=${baseFolder},prodId=${productionId},rootMacro=${rootMacro}