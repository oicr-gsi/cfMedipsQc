# cfMedipsQc

Medips

## Overview

## Dependencies

* [samtools 1.9](https://github.com/samtools/samtools)
* [bowtie2 2.1.0](https://github.com/BenLangmead/bowtie)
* [picard 2.21.2](https://github.com/broadinstitute/picard/)
* [rstats 3.5](https://www.r-project.org/)
* [python 3.6](https://www.python.org/)
* [cfmedips 1.5](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [trimmomatic 0.39](https://github.com/timflutre/trimmomatic)
* [bedops 2.4.37](https://github.com/bedops/bedops)
* [bc 2.1.3](https://github.com/gavinhoward/bc/)
* [hg19-thaliana 1.0](https://gitlab.oicr.on.ca/ResearchIT/modulator)


## Usage

### Cromwell
```
java -jar cromwell.jar run cfMedipsQc.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`fastq1`|File|
`fastq2`|File|


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`window`|Int|300|
`referenceModule`|String|"hg19-thaliana/1.0"|


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`trimming.fastq1Basename`|String|basename("~{fastq1}","_R1_001.fastq.gz")|Basename for the first fastq
`trimming.fastq2Basename`|String|basename("~{fastq2}","_R2_001.fastq.gz")|Basename for the second fastq
`trimming.headCrop`|Int|5|How many bases to crop
`trimming.threads`|Int|6|Requested CPU threads
`trimming.jobMemory`|Int|16|Memory (GB) allocated for this job
`trimming.timeout`|Int|6|Number of hours before task timeout
`trimming.modules`|String|"trimmomatic/0.39"|Module needed to run trimmomatic extract
`alignment.basename`|String|basename("~{fastq1Paired}",".R1_paired.fastq.gz")|Name to make output sam file
`alignment.referenceGenome`|String|"$HG19_THALIANA_ROOT/hg19_thaliana_random"|Using either HG19 or HG38 both with added chromosomes
`alignment.threads`|Int|8|Requested CPU threads
`alignment.jobMemory`|Int|16|Memory (GB) allocated for this job
`alignment.timeout`|Int|6|Number of hours before task timeout
`alignment.modules`|String|"bowtie2/2.1.0 ~{referenceModule}"|Module needed to run bowtie2 alignment
`preprocessing.basename`|String|basename("~{samFile}",".sam")|Name to make output files
`preprocessing.threads`|Int|8|Requested CPU threads
`preprocessing.jobMemory`|Int|16|Memory (GB) allocated for this job
`preprocessing.timeout`|Int|6|Number of hours before task timeout
`preprocessing.modules`|String|"samtools/1.9 picard/2.21.2"|Module needed to run preprocessing
`alignmentMetrics.basename`|String|basename("~{dedupBam}",".sorted.dedup.bam")|Name to make output files
`alignmentMetrics.referenceGenome`|String|"$HG19_THALIANA_ROOT/hg19_thaliana_random"|Using either HG19 or HG38 both with added chromosomes
`alignmentMetrics.threads`|Int|8|Requested CPU threads
`alignmentMetrics.jobMemory`|Int|16|Memory (GB) allocated for this job
`alignmentMetrics.timeout`|Int|6|Number of hours before task timeout
`alignmentMetrics.modules`|String|"samtools/1.9 picard/2.21.2 ~{referenceModule} bc/2.1.3 rstats/3.5"|Modules needed to run alignment metrics
`extractMedipsCounts.basename`|String|basename("~{dedupBam}",".sorted.dedup.bam")|
`extractMedipsCounts.threads`|Int|8|Requested CPU threads
`extractMedipsCounts.jobMemory`|Int|16|Memory (GB) allocated for this job
`extractMedipsCounts.timeout`|Int|6|Number of hours before task timeout
`extractMedipsCounts.modules`|String|"rstats/3.5 cfmedips/1.5 bedops/2.4.37"|Modules needed to run alignment metrics
`finalMetrics.threads`|Int|8|Requested CPU threads
`finalMetrics.jobMemory`|Int|16|Memory (GB) allocated for this job
`finalMetrics.timeout`|Int|6|Number of hours before task timeout
`finalMetrics.modules`|String|"cfmedips/1.5"|Modules needed to run alignment metrics


### Outputs

Output | Type | Description
---|---|---
`outputAlignmentSummaryMetrics`|File|Metrics for alignments
`outputBaseDistributionMetrics`|File|Metrics for base distributions
`outputInsertSizeMetrics`|File|Metrics for insert size
`outputQualityByCycleMetrics`|File|Quality by cycle metrics
`outputQualityDistributionMetrics`|File|Quality distribution metrics
`outputGcBiasMetrics`|File|gc bias metrics
`outputSummaryGcBiasMetrics`|File|Summary of gc bias metrics
`outputThaliaSummary`|File|Summary of thalia chromosomes
`outputDedupBam`|File|De-duplicated bam file
`outputqcMetrics`|File|Final output


## Niassa + Cromwell

This WDL workflow is wrapped in a Niassa workflow (https://github.com/oicr-gsi/pipedev/tree/master/pipedev-niassa-cromwell-workflow) so that it can used with the Niassa metadata tracking system (https://github.com/oicr-gsi/niassa).

* Building
```
mvn clean install
```

* Testing
```
mvn clean verify \
-Djava_opts="-Xmx1g -XX:+UseG1GC -XX:+UseStringDeduplication" \
-DrunTestThreads=2 \
-DskipITs=false \
-DskipRunITs=false \
-DworkingDirectory=/path/to/tmp/ \
-DschedulingHost=niassa_oozie_host \
-DwebserviceUrl=http://niassa-url:8080 \
-DwebserviceUser=niassa_user \
-DwebservicePassword=niassa_user_password \
-Dcromwell-host=http://cromwell-url:8000
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with wdl_doc_gen (https://github.com/oicr-gsi/wdl_doc_gen/)_
