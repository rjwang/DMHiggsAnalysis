#!/bin/bash

#determine CMSSW config
SCRIPT=$(readlink -f $0)
SCRIPTPATH=`dirname $SCRIPT`
BASE=${SCRIPTPATH%DMHiggsAnalysis*}


#configure environment
cd $BASE
mkdir -p tmp/`whoami`/output_$SCRIPTPATH
cd /tmp/`whoami`/output_$SCRIPTPATH
sh rcSetup.sh


#run with the arguments passed
echo $1 + $2
$1 $2
cd /tmp/`whoami`/
rm -rf /tmp/`whoami`/output_$SCRIPTPATH
cd $BASE
