#!/usr/bin/env nextflow
// nf-mrd-umi entry router. Pipeline B (mrd_monitor) is implemented; Pipeline A
// (panel_design) is on the roadmap (see TODO.md).
nextflow.enable.dsl = 2

include { MRD_MONITOR  } from './workflows/mrd_monitor'
include { PANEL_DESIGN } from './workflows/panel_design'

workflow {
    if( params.workflow == 'mrd_monitor' ) {
        MRD_MONITOR()
    } else if( params.workflow == 'panel_design' ) {
        PANEL_DESIGN()
    } else {
        error "Set --workflow mrd_monitor | panel_design (got: '${params.workflow}')"
    }
}
