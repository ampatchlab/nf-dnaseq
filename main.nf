#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
 *
 * ampatchlab/nf-dnaseq: DNA variant calling and annotation Nextflow pipeline
 *
 * Copyright (C) 2021 QIMR Berghofer Medical Research Institute
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


nextflow.enable.dsl=2

import nextflow.config.ConfigParser

nextflow_config = file( "${baseDir}/nextflow.config" ).text
parsed_config = new ConfigParser().setIgnoreIncludes( true ).parse( nextflow_config )
defaults = parsed_config.params

check_params()


/*
 * Imports
 */

// functions
include { parse_readgroup_csv } from './functions/input_csv_parsers.nf' params( params )
include { parse_somatic_csv } from './functions/input_csv_parsers.nf' params( params )
include { parse_germline_csv } from './functions/input_csv_parsers.nf' params( params )

// modules
include { bwa_index } from './modules/bwakit.nf' params( params )
include { gunzip as gunzip_fasta } from './modules/gzip.nf' params( params )
include { multiqc } from './modules/multiqc.nf' params( params )
include { samtools_faidx } from './modules/samtools.nf' params( params )

include { mosdepth } from './modules/mosdepth.nf' params( params )
include { qualimap } from './modules/qualimap.nf' params( params )

include { bcftools_subset_regions as subset_germline_variants } from './modules/bcftools.nf' params( params )
include { bcftools_subset_regions as subset_somatic_variants } from './modules/bcftools.nf' params( params )

include { extract as unpack_vep_cache } from './modules/tar.nf' params( params )

// workflows
include { dna_alignment } from './workflows/alignment.nf' params( params )
include { germline_variant_calling } from './workflows/variant_calling.nf' params( params )
include { somatic_variant_calling } from './workflows/variant_calling.nf' params( params )

include { ensembl_vep as germline_annotation } from './workflows/annotation.nf' params( params )
include { ensembl_vep as somatic_annotation } from './workflows/annotation.nf' params( params )


/*
 * Params
 */

// Cutadapt adapter files
params.r1_adapter_file = params.adapters in params.adapter_files
    ? params.adapter_files[ params.adapters ].r1_adapters
    : "${baseDir}/resource-adapters/null-1.fa"
params.r2_adapter_file = params.adapters in params.adapter_files
    ? params.adapter_files[ params.adapters ].r2_adapters
    : "${baseDir}/resource-adapters/null-2.fa"

// Reference genome files
params.ref_fasta = params.genome in params.genomes
    ? params.genomes[ params.genome ].ref_fasta
    : null
params.vep_cache = params.genome in params.genomes
    ? params.genomes[ params.genome ].vep_cache
    : null

// Mosdepth BED file
params.mosdepth_bed_file = "${baseDir}/assets/null"

// QualiMap feature files
params.qualimap_feature_file = "${baseDir}/assets/null"

// Manta and Strelka call regions
params.call_regions = "${baseDir}/assets/null"

params.manta_call_regions = params.call_regions
params.strelka_call_regions = params.call_regions

// VEP chromosome synonyms
params.vep_synonyms = "${baseDir}/assets/null"


/*
 * Workflow
 */

workflow {

    readgroup_inputs = parse_readgroup_csv( params.readgroup_csv )

    cutadapt_adapter_files = [ params.r1_adapter_file, params.r2_adapter_file ]

    ref_fasta = gunzip_fasta( params.ref_fasta ).ifEmpty( params.ref_fasta )


    // STEP 1 - Index the reference FASTA file
    samtools_faidx( ref_fasta )
    bwa_index( ref_fasta )

    indexed_ref_fasta = ref_fasta \
        | concat( samtools_faidx.out ) \
        | collect


    // STEP 2 - Perform adapter trimming and alignment to the reference
    dna_alignment(
        readgroup_inputs,
        bwa_index.out,
        cutadapt_adapter_files,
    )


    // STEP 3 - Run Mosdepth
    mosdepth_inputs = ! params.skip_mosdepth
        ? dna_alignment.out.alignments.map { sample, indexed_bam -> indexed_bam }
        : Channel.empty()

    mosdepth( mosdepth_inputs, params.mosdepth_bed_file )


    // STEP 4 - Run QualiMap
    qualimap_inputs = ! params.skip_qualimap
        ? dna_alignment.out.alignments.map { sample, indexed_bam -> indexed_bam.first() }
        : Channel.empty()

    qualimap( qualimap_inputs, params.qualimap_feature_file )


    // STEP 5 - Extract the indexed Ensembl VEP cache files
    vep_cache = ( params.germline_csv || params.somatic_csv )
        ? unpack_vep_cache( params.vep_cache )
        : Channel.empty()

    vep_species = ( params.germline_csv || params.somatic_csv )
        ? vep_cache.map { it.baseName }
        : Channel.empty()


    // STEP 6 - Call and annotate germline variants
    if( params.germline_csv ) {

        parse_germline_csv( params.germline_csv ) \
            | map { analysis, samples ->
                tuple( groupKey(analysis, samples.size()), samples )
            } \
            | transpose() \
            | map { analysis, sample -> tuple( sample, analysis ) } \
            | combine( dna_alignment.out.alignments, by: 0 ) \
            | map { sample, analysis, indexed_bam -> tuple( analysis, indexed_bam ) } \
            | groupTuple() \
            | map { analysis, indexed_bam_files ->
                tuple( analysis.toString(), indexed_bam_files.flatten() )
            } \
            | set { germline_inputs }

        germline_variant_calling(
            germline_inputs,
            indexed_ref_fasta,
            params.strelka_call_regions,
        )

        germline_annotation(
            germline_variant_calling.out,
            indexed_ref_fasta,
            vep_cache,
            params.vep_synonyms,
            vep_species,
        )
    }


    // STEP 7 - Call and annotate somatic variants
    if( params.somatic_csv ) {

        test_control_inputs = parse_somatic_csv( params.somatic_csv )

        test_control_inputs \
            | map { analysis, test, control ->
                tuple( test, tuple( analysis, test, control ) )
            } \
            | combine( dna_alignment.out.alignments, by: 0 ) \
            | map { sample, key_tuple, indexed_bam -> tuple( *key_tuple, indexed_bam ) } \
            | set { test_inputs }

        test_control_inputs \
            | map { analysis, test, control ->
                tuple( control, tuple( analysis, test, control ) )
            } \
            | combine( dna_alignment.out.alignments, by: 0 ) \
            | map { sample, key_tuple, indexed_bam -> tuple( *key_tuple, indexed_bam ) } \
            | set { control_inputs }

        test_control_inputs \
            | join( test_inputs, by: [0,1,2] ) \
            | join( control_inputs, by: [0,1,2] ) \
            | map { analysis, test, control, indexed_test_bam, indexed_control_bam ->
                tuple( analysis, indexed_test_bam, indexed_control_bam )
            } \
            | set { somatic_inputs }

        somatic_variant_calling(
            somatic_inputs,
            indexed_ref_fasta,
            params.manta_call_regions,
            params.strelka_call_regions,
        )

        somatic_annotation(
            somatic_variant_calling.out,
            indexed_ref_fasta,
            vep_cache,
            params.vep_synonyms,
            vep_species,
        )
    }


    // STEP 8 - Create a MultiQC report
    logs = Channel.empty() \
        | mix( dna_alignment.out.logs ) \
        | mix( mosdepth.out.dists ) \
        | mix( mosdepth.out.summary ) \
        | mix( qualimap.out ) \
        | collect

    multiqc( logs, params.multiqc_config )
}


workflow.onComplete {

    log.info "Workflow completed at: ${workflow.complete}"
    log.info "Time taken: ${workflow.duration}"
    log.info "Execution status: ${workflow.success ? 'success' : 'failed'}"
    log.info "Output directory: ${params.publish_dir}"
}


workflow.onError {

    log.info "Execution halted: ${workflow.errorMessage}"
}


/*
 * Functions
 */

def check_params() {

    if( params.help ) {
        usage()
        exit 0
    }

    if( params.version ) {
        log.info( workflow.manifest.version )
        exit 0
    }

    if( !params.readgroup_csv ) {
        log.info "Readgroup CSV not specified. Please using the `--readgroup_csv` parameter"
        exit 1
    }
}

def usage() {

    log.info"""
    Usage:
        nextflow run -profile <profile> -revision <revision> ampatchlab/nf-dnaseq [options]


    Nextflow execution options:

        -profile STR
            Nextflow configuration profile to use. Available profiles include:
            'awsbatch', 'conda', 'docker' and 'singularity'

        -revision STR
            Git branch/tag (version) of this workflow to use

        -work-dir DIR
            Directory where intermediate result files are stored

        -help
            Show additional execution options and exit


    Workflow input params:

        --readgroup_csv FILE
            Comma-separated list of sample and readgroup inputs
            Please see: https://github.com/ampatchlab/nf-dnaseq#readgroup-inputs

        --germline_csv FILE
            Comma-separated list of analysis identifiers and sample inputs
            Please see: https://github.com/ampatchlab/nf-dnaseq#germline-inputs

        --somatic_csv FILE
            Comma-separated list of analysis identifiers and test and control sample inputs
            Please see: https://github.com/ampatchlab/nf-dnaseq#somatic-inputs


    Reference genome params:

        --genome STR
            Genome name [Either: ${defaults.genomes.keySet().join(", ")}; Default: ${defaults.genome}]

        --ref_fasta FILE
            Override the reference FASTA file with FILE [Default: ${defaults.ref_fasta ?: null}]

        --vep_cache FILE
            Override the VEP cache with FILE [Default: ${defaults.vep_cache ?: null}]


    Adapter trimming params:

        --adapters STR
            The adapters to trim [Either: ${defaults.adapter_files.keySet().join(", ")}; Default: ${defaults.adapters}]

        --r1_adapter_file FILE
            Override the R1 adapter file with FILE [Default: ${defaults.r1_adapter_file ?: null}]

        --r2_adapter_file FILE
            Override the R2 adapter file with FILE [Default: ${defaults.r2_adapter_file ?: null}]


    Mosdepth params:

        --mosdepth_bed_file FILE
            Override the Mosdepth BED file of regions [Default: ${defaults.mosdepth_bed_file ?: null}]

        --skip_mosdepth
            Skips Mosdepth process execution


    Qualimap params:

        --qualimap_feature_file FILE
            Override the QualiMap feature file of regions in GFF/GTF or BED format [Default: ${defaults.qualimap_feature_file ?: null}]

        --skip_qualimap
            Skips QualiMap process execution


    Strelka and Manta params:

        --call_regions FILE
            Restrict variant calling to the regions in BED FILE [Default: ${defaults.call_regions ?: null}]

        --exome
            Provide appropriate settings for WES and other regional enrichment analyses


    Ensembl-VEP params:

        --vep_cache_type STR
            Apply when using a merged or refseq VEP cache [Either: 'merged', 'refseq'; Default: ${defaults.vep_cache_type ?: null}]

        --vep_synonyms FILE
            Load chromosome synonyms from FILE [Default: ${defaults.vep_synonyms ?: null}]


    Output params:

        --publish_dir DIR
            Path where the results will be published [Default: ${defaults.publish_dir}]

        --publish_mode STR
            The mode used to publish files to the target directory [Default: ${defaults.publish_mode}]


    Standard params:

        --help
            Show this message and exit

        --version
            Show the pipeline version and exit
    """.stripIndent()
}
