 
#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#enter the workflow's final output directory ($1)
cd $1

#find all files, return their md5sums to std out
find . \( -type f -size +0 -iname "qc_metrics.json" \) -printf "json file exists %f\n";
