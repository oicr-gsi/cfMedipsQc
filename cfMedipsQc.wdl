version 1.0
workflow cfMedipsQc {
  input {
      File fastq1
      File fastq2
      Int window = 300
      String referenceGenome
      String referenceModule
      String fastqFormat
  }
  
  call trimming {
    input: fastq1 = fastq1,
           fastq2 = fastq2,
           fastqFormat = fastqFormat
  }  

  call alignment {
    input: fastq1Paired = trimming.outputFastq1Paired,
           fastq2Paired = trimming.outputFastq2Paired,
           referenceGenome = referenceGenome,
           referenceModule = referenceModule
  }

  call preprocessing {
    input: samFile = alignment.samFile
  }

  call alignmentMetrics {
    input: dedupBam = preprocessing.dedupBam,
           referenceGenome = referenceGenome,
           referenceModule = referenceModule
  }

  call extractMedipsCounts {
    input: dedupBam = preprocessing.dedupBam,
           metricsDedup = preprocessing.metricsDedup,
           summaryGcBiasMetrics = alignmentMetrics.summaryGcBiasMetrics,
           alignmentSummaryMetrics = alignmentMetrics.alignmentSummaryMetrics,
           thaliaSummary = alignmentMetrics.thaliaSummary,
           window = window
  }
  call finalMetrics {
    input: dedupBam = preprocessing.dedupBam,
           metricsDedup = preprocessing.metricsDedup,
           summaryGcBiasMetrics = alignmentMetrics.summaryGcBiasMetrics,
           alignmentSummaryMetrics = alignmentMetrics.alignmentSummaryMetrics,
           thaliaSummary = alignmentMetrics.thaliaSummary,
           coverageWindows = extractMedipsCounts.coverageWindows,
           enrichmentData = extractMedipsCounts.enrichmentData,
           coverageCounts = extractMedipsCounts.coverageCounts, 
           saturationMetrics = extractMedipsCounts.saturationMetrics,
           nameFile = extractMedipsCounts.nameFile
  }

  parameter_meta {
    fastq1: "Read 1 input fastq file"
    fastq2: "Read 2 input fastq file"
    fastqFormat: "Quality encoding, default is phred33, but can be set to phred64"
    window:  "window length, over which to assess"
    referenceGenome: "reference genome to use"
    referenceModule: "module to load the reference genome"
  }

  meta {
    author: "Rishi Shah"
    email: "rshah@oicr.on.ca"
    description: "Medips"
    dependencies: 
    [
      {
      name: "samtools/1.9",
      url: "https://github.com/samtools/samtools"
      },
      {
      name: "bowtie2/2.1.0",
      url: "https://github.com/BenLangmead/bowtie"
      },
      {
      name: "picard/2.21.2",
      url: "https://github.com/broadinstitute/picard/"
      },
      {
      name: "rstats/3.5",
      url: "https://www.r-project.org/"
      },
      {
      name: "python/3.6",
      url: "https://www.python.org/"
      },
      {
      name: "cfmedips/1.5",
      url: "https://github.com/oicr-gsi/medips-tools.git"
      },
      {
      name: "trimmomatic/0.39",
      url: "https://github.com/timflutre/trimmomatic"
      },
      {
      name: "bedops/2.4.37",
      url: "https://github.com/bedops/bedops"
      },
      {
      name: "bc/2.1.3",
      url: "https://github.com/gavinhoward/bc/"
      }
    ]   
    output_meta: {
       outputAlignmentSummaryMetrics: "Metric for alignments",
       outputBaseDistributionMetrics: "Metrics for base distributions",
       outputInsertSizeMetrics: "Metrics for insert size",
       outputQualityByCycleMetrics: "Quality by cycle metrics",
       outputQualityDistributionMetrics: "Quality distribution metrics",
       outputGcBiasMetrics: "gc bias metrics",
       outputSummaryGcBiasMetric: "Summary of gc bias metrics",
       outputThaliaSummary: "Summary of thalia chromosomes",
       outputDedupBam: "De-duplicated bam file",
       outputqcMetrics: "Final output"
    }
  }
  output {
    File outputAlignmentSummaryMetrics = alignmentMetrics.alignmentSummaryMetrics 
    File outputBaseDistributionMetrics = alignmentMetrics.baseDistributionMetrics
    File outputInsertSizeMetrics = alignmentMetrics.insertSizeMetrics
    File outputQualityByCycleMetrics = alignmentMetrics.qualityByCycleMetrics
    File outputQualityDistributionMetrics = alignmentMetrics.qualityDistributionMetrics
    File outputGcBiasMetrics = alignmentMetrics.gcBiasMetrics
    File outputSummaryGcBiasMetrics = alignmentMetrics.summaryGcBiasMetrics
    File outputThaliaSummary = alignmentMetrics.thaliaSummary
    File outputDedupBam = preprocessing.dedupBam
    File outputqcMetrics = finalMetrics.qcMetrics
  }
}
task trimming {
  input {
    File fastq1
    File fastq2
    String fastq1Basename =basename("~{fastq1}", "_R1_001.fastq.gz")
    String fastq2Basename =basename("~{fastq2}", "_R2_001.fastq.gz")
    String fastqFormat
    Int headCrop = 5
    Int threads = 6
    Int jobMemory = 16
    Int timeout = 6  
    String modules = "trimmomatic/0.39"
  }

  parameter_meta {
    fastq1: "First fastq input file containing reads"
    fastq2: "Second fastq input file containing reads"
    fastq1Basename: "Basename for the first fastq"
    fastq2Basename: "Basename for the second fastq"
    fastqFormat: "Quality encoding, default is phred33, but can be set to phred64"
    headCrop: "How many bases to crop"
    modules: "Module needed to run trimmomatic extract"
    jobMemory: "Memory (GB) allocated for this job"
    threads: "Requested CPU threads"
    timeout: "Number of hours before task timeout"
  }
  command <<<
    set -euo pipefail
    trimmomatic PE \
                ~{fastq1} ~{fastq2} \
                "-~{fastqFormat}" \
                "~{fastq1Basename}.R1_paired.fastq.gz" "~{fastq1Basename}.R1_unpaired.fastq.gz" "~{fastq2Basename}.R2_paired.fastq.gz" "~{fastq2Basename}.R2_unpaired.fastq.gz" \
                HEADCROP:~{headCrop}
  >>>

  runtime {
    modules: "~{modules}"
    memory:  "~{jobMemory} GB"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  }

  output {
    File outputFastq1Paired = "~{fastq1Basename}.R1_paired.fastq.gz"
    File outputFastq1Unpaired = "~{fastq1Basename}.R1_unpaired.fastq.gz"
    File outputFastq2Paired = "~{fastq2Basename}.R2_paired.fastq.gz"
    File outputFastq2Unpaired = "~{fastq2Basename}.R2_unpaired.fastq.gz"
  }

  meta {
    output_meta: {
      outputFastq1Paired: "Read 1 Paired",
      outputFastq1Unpaired: "Read 1 Unpaired",
      outputFastq2Paired: "Read 2 Paired",
      outputFastq2Unpaired: "Read 2 Unpaired"
    }
  }
}

task alignment {
  input {
    File fastq1Paired
    File fastq2Paired
    String basename = basename("~{fastq1Paired}", ".R1_paired.fastq.gz")
    String referenceModule
    String referenceGenome
    Int threads = 8
    Int jobMemory = 16
    Int timeout = 6  
    String modules = "bowtie2/2.1.0 ~{referenceModule}"
  } 
  parameter_meta {
    fastq1Paired: "First fastq input file containing reads"
    fastq2Paired: "Second fastq input file containing reads"
    basename: "Name to make output sam file"
    referenceGenome: "Using either HG19 or HG38 both with added chromosomes"
    modules: "Module needed to run bowtie2 alignment"
    jobMemory: "Memory (GB) allocated for this job"
    threads: "Requested CPU threads"
    timeout: "Number of hours before task timeout"
  }
  command <<<
    set -euo pipefail
    bowtie2 -p 8 -x "~{sub(referenceGenome, "\.fa(sta)?$", "")}" \
            -1 ~{fastq1Paired} \
            -2 ~{fastq2Paired} \
            -S "~{basename}.sam"
  >>>

  runtime {
    modules: "~{modules}"
    memory:  "~{jobMemory} GB"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  }
  output {
    File samFile = "~{basename}.sam"
  }

  meta {
    output_meta: {
     samFile: "Stores the aligned genome"
    }
  }
}

task preprocessing {
  input {
    File samFile
    String basename =basename("~{samFile}", ".sam")
    Int threads = 8
    Int jobMemory = 16
    Int timeout = 6  
    String modules = "samtools/1.9 picard/2.21.2"
  } 

  parameter_meta {
    samFile: "Sam file"
    basename: "Name to make output files"
    modules: "Module needed to run preprocessing"
    jobMemory: "Memory (GB) allocated for this job"
    threads: "Requested CPU threads"
    timeout: "Number of hours before task timeout"
  }

  command <<<
    set -euo pipefail
    samtools view -bS ~{samFile} | samtools sort -o "~{basename}.sorted.bam"
    java -jar ${PICARD_ROOT}/picard.jar MarkDuplicates \
      I="~{basename}.sorted.bam" \
      O="~{basename}.sorted.dedup.bam" \
      M="~{basename}.sorted.dedup.metrics" \
      ASSUME_SORTED=true \
      VALIDATION_STRINGENCY=SILENT \
      REMOVE_DUPLICATES=false
  >>>
  runtime {
    modules: "~{modules}"
    memory:  "~{jobMemory} GB"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  } 

  output {
    File bamFile = "~{basename}.sorted.bam"
    File dedupBam = "~{basename}.sorted.dedup.bam"
    File metricsDedup = "~{basename}.sorted.dedup.metrics"
  }
  meta {
    output_meta: {
     bamFile: "Stores the aligned genome",
     dedupBam: "Stores the genome de-duplicated",
     metricsDedup: "Stores the dedup metrics"
    }
  }

}

task alignmentMetrics {
  input {
    File dedupBam
    String basename =basename("~{dedupBam}", ".sorted.dedup.bam")
    String referenceModule
    String referenceGenome
    Int threads = 8
    Int jobMemory = 32
    Int timeout = 6  
    String modules = "samtools/1.9 picard/2.21.2 ~{referenceModule} bc/2.1.3 rstats/3.5"
  } 
 parameter_meta {
    dedupBam: "De-Duplicated Bam file"
    basename: "Name to make output files"
    referenceGenome: "Using either HG19 or HG38 both with added chromosomes"
    modules: "Modules needed to run alignment metrics"
    jobMemory: "Memory (GB) allocated for this job"
    threads: "Requested CPU threads"
    timeout: "Number of hours before task timeout"
  }

  command <<< 
    set -euo pipefail
    java -jar ${PICARD_ROOT}/picard.jar CollectMultipleMetrics \
      R="~{referenceGenome}"\
      I=~{dedupBam} \
      O="~{basename}" \
      VALIDATION_STRINGENCY=SILENT
    java -jar ${PICARD_ROOT}/picard.jar CollectGcBiasMetrics \
      R="~{referenceGenome}" \
      I=~{dedupBam} \
      O="~{basename}.gc_bias_metrics.txt" \
      S="~{basename}.summary_gc_bias_metrics.txt" \
      CHART="~{basename}.gc_bias_metrics.pdf"
    samtools view ~{dedupBam} | cut -f 3 | sort | uniq -c | sort -nr | sed -e 's/^ *//;s/ /\t/' | awk 'OFS="\t" {print $2,$1}' | sort -n -k1,1 > thalia.counts
    total=$(samtools view ~{dedupBam} | wc -l)
    unmap=$(cat thalia.counts | grep "^\*" | cut -f2); if [[ -z $unmap ]]; then unmap="0"; fi
    if [[ $(cat thalia.counts | grep F19K16 | cut -f2) == "" ]]; then methyl=0; else methyl=$(cat thalia.counts | grep F19K16 | cut -f2); fi
    if [[ $(cat thalia.counts | grep F24B22 | cut -f2) == "" ]]; then unmeth=0; else unmeth=$(cat thalia.counts | grep F24B22 | cut -f2); fi
    if [[ $total == "" ]]; then pct_thalia=0; else pct_thalia=$(echo "scale=3; ($methyl + $unmeth)/$total * 100" | bc -l); fi
    if [[ -z $pct_thalia ]]; then pct_thalia="0"; fi
    if [[ $methyl == 0 && $unmeth == 0 ]]; then bet_thalia=0; else bet_thalia=$(echo "scale=3; $methyl/($methyl + $unmeth)" | bc -l); fi 
    if [[ -z $bet_thalia ]]; then bet_thalia="0"; fi
    echo -e "total\tunmap\tmethyl\tunmeth\tPCT_THALIANA\tTHALIANA_BETA" > thalia_summary.txt
    echo -e "$total\t$unmap\t$methyl\t$unmeth\t$pct_thalia\t$bet_thalia" >> thalia_summary.txt
  >>>
 
  runtime {
    modules: "~{modules}"
    memory:  "~{jobMemory} GB"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  } 

  output {
    File alignmentSummaryMetrics = "~{basename}.alignment_summary_metrics"
    File baseDistributionMetrics = "~{basename}.base_distribution_by_cycle_metrics"
    File insertSizeMetrics = "~{basename}.insert_size_metrics"
    File qualityByCycleMetrics = "~{basename}.quality_by_cycle_metrics"
    File qualityDistributionMetrics = "~{basename}.quality_distribution_metrics"
    File gcBiasMetrics = "~{basename}.gc_bias_metrics.txt"
    File summaryGcBiasMetrics = "~{basename}.summary_gc_bias_metrics.txt"
    File thaliaSummary = "thalia_summary.txt"

  }

  meta {
    output_meta: {
     alignmentSummaryMetrics: "Stores the alignment summary metrics",
     baseDistributionMetrics: "Stores the base distribution metrics",
     insertSizeMetrics: "Stores the insert size metrics",
     qualityByCycleMetrics: "Stores the quality by cycle metrics",
     qualityDistributionMetrics: "Stores the quality by distribution metrics",
     gcBiasMetrics: "Stores the gc bias metrics",
     summaryGcBiasMetrics: "Summary of gc bias metrics",
     thaliaSummary: "Summary of thalia metrics"
    }
  }

}

task extractMedipsCounts {
  input {
    File dedupBam
    File metricsDedup
    File summaryGcBiasMetrics
    File alignmentSummaryMetrics
    File thaliaSummary
    Int window
    String basename =basename("~{dedupBam}", ".sorted.dedup.bam")
    String convert2bed = "$BEDOPS_ROOT/convert2bed"
    Int threads = 8
    Int jobMemory = 32
    Int timeout = 6  
    String modules = "rstats/3.5 cfmedips/1.5 bedops/2.4.37"
  }
  parameter_meta {
    dedupBam: "Dedup bam file"
    metricsDedup: "Metrics of dedup bam file"
    summaryGcBiasMetrics: "GC metrics summary"
    alignmentSummaryMetrics: "Alignment summary metrics"
    thaliaSummary: "Summary of the thalia data"
    window: "value of window"
    basename: "basename for the sample"
    convert2bed: "path to conver2bed program"
    modules: "Modules needed to run alignment metrics"
    jobMemory: "Memory (GB) allocated for this job"
    threads: "Requested CPU threads"
    timeout: "Number of hours before task timeout"

  }
  command <<<
    set -euo pipefail
      medips.R \
        --basedir . \
        --bamfile ~{dedupBam} \
        --samplename ~{basename} \
        --ws  ~{window}\
        --outdir .
      NAME=""
      count0=$(awk '$1 == 0' genome_count.txt | wc -l)
      count1=$(awk '$1 >= 1' genome_count.txt | wc -l)
      count10=$(awk '$1 >= 10' genome_count.txt | wc -l)
      count50=$(awk '$1 >= 50' genome_count.txt | wc -l)
      count100=$(awk '$1 >= 100' genome_count.txt | wc -l)
      echo -e "sample\tcount0\tcount1\tcount10\tcount50\tcount100" > coverage_windows.txt
      echo -e "$NAME\t$count0\t$count1\t$count10\t$count50\t$count100" >> coverage_windows.txt
      echo -e "samples\n~{basename}" > name.txt
      ~{convert2bed} -d --input wig < medips.wig > medips.bed

  >>>
  runtime {
    modules: "~{modules}"
    memory:  "~{jobMemory} GB"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  } 
  output {
    File coverageCounts = "coverage_counts.txt"
    File enrichmentData = "enrichment_data.txt"
    File genomeCount = "genome_count.txt"
    File windowPerChr = "MEDIPS_window_per_chr.csv"
    File saturationMetrics = "saturation_metrics.txt"
    File medipsBed = "medips.bed"
    File nameFile = "name.txt"
    File coverageWindows = "coverage_windows.txt"
  }

  meta {
    output_meta: {
      coverageCounts: "Text files with the coverage counts",
      enrichmentData: "Text file with enrichment data",
      genomeCount: "Text file with the genome count",
      windowPerChr: "CSV file which shows window per chromosome",
      saturationMetrics: "Text file with saturation metrics",
      medipsBed: "Bed file with the medips data",
      nameFile: "File that contains the sample name",
      coverageWindows: "File that contains the coverage windows"
    }
  }
}  
task finalMetrics {
  input {
    File dedupBam
    File metricsDedup
    File summaryGcBiasMetrics
    File alignmentSummaryMetrics
    File thaliaSummary
    File coverageCounts
    File enrichmentData
    File coverageWindows
    File saturationMetrics
    File nameFile
    Int threads = 8
    Int jobMemory = 16
    Int timeout = 6
    String modules = "cfmedips/1.5"
  }
  parameter_meta {
    dedupBam: "Dedup bam file"
    metricsDedup: "Metrics of dedup bam file"
    summaryGcBiasMetrics: "GC metrics summary"
    alignmentSummaryMetrics: "Alignment summary metrics"
    thaliaSummary: "Summary of the thalia data"
    coverageCounts: "Text files with the coverage counts"
    enrichmentData: "Text file with enrichment data"
    coverageWindows: "File that contains the coverage windows"
    saturationMetrics: "Text file with saturation metrics"
    nameFile: "File that contains the sample name"
    modules: "Modules needed to run alignment metrics"
    jobMemory: "Memory (GB) allocated for this job"
    threads: "Requested CPU threads"
    timeout: "Number of hours before task timeout"
  }  
  command <<<
    set -euo pipefail
    txt-to-json.py -n ~{nameFile} \
                   -e ~{enrichmentData} \
                   -c ~{coverageCounts} \
                   -w ~{coverageWindows} \
                   -s ~{saturationMetrics} \
                   -d ~{metricsDedup} \
                   -u ~{summaryGcBiasMetrics} \
                   -a ~{alignmentSummaryMetrics} \
                   -t ~{thaliaSummary}
  >>>
  runtime {
  modules: "~{modules}"
  memory:  "~{jobMemory} GB"
  cpu:     "~{threads}"
  timeout: "~{timeout}"
  }
  output {
    File qcMetrics = "qc_metrics.json"
  } 
  meta {
    output_meta: {
      qcMetrics: "File that combines all the metrics to be used"
    }
  }  
}  
