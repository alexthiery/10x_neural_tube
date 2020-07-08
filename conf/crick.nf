#!/usr/bin/env nextflow

singularity {
  enabled = true
  autoMounts = true
}

process {
  beforeScript = 'module load Singularity/2.6.0-foss-2016b'
  executor = 'slurm'
}