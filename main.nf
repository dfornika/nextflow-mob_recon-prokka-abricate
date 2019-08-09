#!/usr/bin/env nextflow

Channel
    .fromPath(params.input)
    .splitCsv(header:true, sep:'\t', quote:'"')
    .map{ row-> tuple(row.sample_id, file(row.assembly)) }
    .set { samples_mob_recon_ch; }


def summary = [:]
summary['Pipeline Name']  = 'mob_recon-prokka-abricate'
summary['Input']          = params.input
// summary['MOB-Suite DB']   = params.mobdb
summary['Output dir']     = params.outdir
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile


summary.each{ k, v -> println "${k}: ${v}" }

/*
process template {
    tag "$sample_id"
    cpus 1
    conda '/home/dfornika/miniconda3/envs/'
    publishDir "${params.outdir}", mode: 'copy', pattern: ""
    
    input:
    set sample_id, file('') from _ch

    output:
    set sample_id, file('') into _ch
    
    script:
    """
    
    """
}
*/

/*
 * Reconstruct plasmids with mob_recon
 */
process mob_recon {
    tag "$sample_id"
    cpus 8
    conda '/home/dfornika/miniconda3/envs/mob_suite-1.4.9.1'
    publishDir "${params.outdir}", mode: 'copy', pattern: "$sample_id"
 
    input:
    set sample_id, file(assembly) from samples_mob_recon_ch

    output:
    file(sample_id)
    set sample_id, file("${sample_id}/mobtyper_aggregate_report.txt") into add_sample_id_ch
    file("${sample_id}/${sample_id}_plasmid*.fasta") into plasmids_abricate_ch
    file("${sample_id}/${sample_id}_chromosome.fasta") into chromosome_abricate_ch

    script:
    """
    mob_recon \
    --num_threads 8 \
    --unicycler_contigs \
    --run_circlator \
    --run_typer \
    --infile $assembly \
    --outdir $sample_id
    rename plasmid ${sample_id}_plasmid ${sample_id}/plasmid*.fasta
    rename chromosome ${sample_id}_chromosome ${sample_id}/chromosome.fasta
    sed -i 's/^>.*|/>/g' ${sample_id}/*.fasta
    """
}

process add_sample_id_to_mobtyper_aggregate_report {
    tag "$sample_id"
    cpus 1

    input:
    set sample_id, file(mobtyper) from add_sample_id_ch

    output:
    file('mobtyper_aggregate_report_with_sample_id.txt') into concatenate_reports_ch
    
    script:
    """
    head -n 1 "mobtyper_aggregate_report.txt" | sed 's/^/sample_id\\t/' > header.tsv
    tail -n +2 "mobtyper_aggregate_report.txt" | sed 's/^/${sample_id}\\t/' > data.tsv
    cat header.tsv data.tsv > mobtyper_aggregate_report_with_sample_id.txt
    """
}

process concatenate_mobtyper_reports {
    cpus 1
    publishDir "${params.outdir}", mode: 'copy', pattern: "mobtyper_aggregate_report.txt"
    
    input:
    file('mobtyper_aggregate_report_with_sample_id_*.txt') from concatenate_reports_ch.collect()

    output:
    file('mobtyper_aggregate_report.txt')
    
    script:
    """
    head -n 1 mobtyper_aggregate_report_with_sample_id_1.txt > header.txt
    tail -qn +2 mobtyper_aggregate_report_with_sample_id_*.txt > data.txt
    cat header.txt data.txt > mobtyper_aggregate_report.txt
    """
}


/*
process plasmids_prokka {
    cpus 8
    conda '/home/dfornika/miniconda3/envs/prokka-1.14.0'
    publishDir "${params.outdir}/prokka_output", mode: 'copy', pattern: "prokka_output/*.gbk"

    input:
    file(plasmid_fasta) from plasmids_prokka_ch.flatten()

    output:
    file('prokka_output/*.gbk')

    script:
    """
    prokka \
    --centre BCCDC-PHL \
    --compliant \
    --proteins /home/dfornika/analyses/download_plasmid_genbank/plasmids.gbk \
    --outdir prokka_output \
    "${plasmid_fasta}"
    """
}
*/

process plasmids_abricate {
    cpus 12
    conda '/home/dfornika/miniconda3/envs/abricate-0.8.13'
    publishDir "${params.outdir}", mode: 'copy', pattern: "*_abricate_card.tsv"
    
    input:
    file(plasmid_assembly) from plasmids_abricate_ch.flatten()

    output:
    set file(plasmid_assembly), file('*_abricate_card.tsv') into identify_kpc_contig_ch
    
    script:
    """
    abricate \
    --threads 12 \
    --db card \
    --nopath \
    ${plasmid_assembly} \
    > ${plasmid_assembly}_abricate_card.tsv
    rename .fasta_abricate _abricate ${plasmid_assembly}_abricate_card.tsv
    """
}

/*
process chromosome_abricate {
    cpus 12
    conda '/home/dfornika/miniconda3/envs/abricate-0.8.13'
    publishDir "${params.outdir}", mode: 'copy', pattern: "*_abricate_card.tsv"
    
    input:
    file(chromosome_assembly) from chromosome_abricate_ch

    output:
    file('*_abricate_card.tsv')
    
    script:
    """
    abricate \
    --threads 12 \
    --db card \
    ${chromosome_assembly} \
    > ${chromosome_assembly}_abricate_card.tsv
    """
}
*/


process identify_kpc_contig {
    tag "$sample_id"
    cpus 1
    conda '/home/dfornika/miniconda3/envs/seqtk-1.3'
    publishDir "${params.outdir}", mode: 'copy', pattern: "*_kpc.fasta"
    
    input:
    set file(plasmid_assembly), file(abricate_output) from identify_kpc_contig_ch

    output:
    file('*_kpc.fasta') optional true
    
    script:
    """
    grep -h 'KPC' ${abricate_output} | cut -f 2 > name.lst
    seqtk subseq ${plasmid_assembly} name.lst > ${plasmid_assembly}_kpc.fasta
    rename .fasta_kpc.fasta _kpc.fasta *.fasta_kpc.fasta
    find . -size 0 -type f -name *_kpc.fasta -delete
    """
}
