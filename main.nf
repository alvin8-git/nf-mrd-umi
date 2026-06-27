#!/usr/bin/env nextflow
// nf-mrd-umi entry router. Pipeline B (mrd_monitor) is implemented; Pipeline A
// (panel_design) is on the roadmap (see TODO.md).
nextflow.enable.dsl = 2

include { MRD_MONITOR } from './workflows/mrd_monitor'

workflow {
    if( params.workflow == 'mrd_monitor' ) {
        MRD_MONITOR()
    } else if( params.workflow == 'panel_design' ) {
        error "panel_design (Pipeline A) is not implemented yet — see TODO.md"
    } else {
        error "Set --workflow mrd_monitor | panel_design (got: '${params.workflow}')"
    }
}
