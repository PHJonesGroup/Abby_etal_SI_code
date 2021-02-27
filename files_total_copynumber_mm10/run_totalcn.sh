#!/bin/bash

samplename=$1
bamfile=$2
normalname=$3
normalbamfile=$4
kit="final_kit.bed"
mkdir -p output/${samplename}
cmd="R --no-save --no-restore --vanilla -f call_copynumber_offtarget.R --args ${samplename} ${bamfile} ${normalname} ${normalbamfile} ${kit}"
echo ${cmd}
${cmd}
