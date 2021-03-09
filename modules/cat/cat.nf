#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process cat {

    label 'process_low'

    input:
        path(file1)
        path(file2)

    output:
        path('concatenated.gtf')

    script:

        """
        #!/bin/bash

        cat ${file1} ${file2} > concatenated.gtf
        """
}