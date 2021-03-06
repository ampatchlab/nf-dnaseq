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


/*
 * Params
 */

params.publish_dir = './results'
params.publish_everything = false
params.publish_mode = 'copy'

params.publish_multiqc = false


/*
 * Processes
 */

process multiqc {

    label 'multiqc'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_multiqc,
        mode: params.publish_mode,
    )

    input:
    path 'data/*'
    path config

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    """
    multiqc \\
        --config "${config}" \\
        .
    """
}
