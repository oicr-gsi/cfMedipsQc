#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#enter the workflow's final output directory ($1)
cd $1

#get md5sum for output file ignoring two lines
cat qc_metrics.json | grep -v "maxEstCor\"\|maxTruCor\"" | md5sum
