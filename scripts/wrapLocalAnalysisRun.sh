#!/bin/bash

#determine CMSSW config
SCRIPT=$(readlink -f $0)
SCRIPTPATH=`dirname $SCRIPT`
BASE=${SCRIPTPATH%DMHiggsAnalysis*}


#configure environment
cd $BASE
sh rcSetup.sh

#run with the arguments passed
echo $1 + $2
$1 $2
