# 1.2.0 - 2022-06-10
[GRD-635](https://jira.oicr.on.ca/browse/GRD-635) : moved assembly-specific choices into the workflow, updated README.md
# 1.1.3 - 2022-12-22
Patching aggregate function so it handle empty inputs
# 1.1.2 - 2022-12-23
[GDI-2555](https://jira.oicr.on.ca/browse/GDI-2555) : introduced a parameter to set a minimum threshold for reads in input bam files
# 1.1.1 - 2022-12-07
[GDI-2555](https://jira.oicr.on.ca/browse/GDI-2555) : scattering by chromosome also brought in a possibility to encounter empty inputs in shards of extractMedipsCounts task. This is a fix to that.
# 1.1.0 - 2022-11-20
[GDI-2615](https://jira.oicr.on.ca/browse/GDI-2615) : scattering the task extractMedipsCounts by chromosome, updated default module for medips script. Additional tasks for handling splitted outputs
# 1.0.7 - 2022-06-07
[GP-3386](https://jira.oicr.on.ca/browse/GP-3386) : restoring reference argument for extract medips counts task
# 1.0.6 - 2022-05-13
[GP-2942](https://jira.oicr.on.ca/browse/GP-2942) : making insert size metrics an optional output (to prevent fails when not enough data are available)
## 1.0.5 - 2022-04-19
[GC-9245](https://jira.oicr.on.ca/browse/GC-9245) : small bug in the final shell script, to deal with grep causing pipefail
[GP-2872](https://jira.oicr.on.ca/browse/GP-2872) Making RT scripts more robust
## 1.0.4 - 2021-06-01
- Migration to Vidarr
