process BOWTIE2_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.5.1--py39hfa7eb30_4' :
        'quay.io/biocontainers/bowtie2:2.5.1--py39hfa7eb30_4' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bam")                   , emit: bam
    tuple val(meta), path("*.log")                   , emit: summary
    tuple val(meta), path("*fastq.gz"), optional:true, emit: fastq
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def seq_center = params.seq_center ? "--rg-id ${prefix} --rg SM:${prefix} --rg CN:${params.seq_center.replaceAll('\\s','_')}" : "--rg-id ${prefix} --rg SM:${prefix}"
    def unaligned_opt = ''

    if (meta.single_end) {
        unaligned_opt = params.save_unaligned || params.contaminant_screening ? "--un-gz ${prefix}.unmapped.fastq.gz" : ''
        """
        INDEX=`find -L ./ -name "*.1.bt2" -o -name "*.1.bt2l" | head -n 1 | sed 's/\\.1.bt2l\$//' | sed 's/\\.1.bt2\$//'`
        bowtie2 \\
            -x \$INDEX \\
            -U ${reads} \\
            --threads ${task.cpus} \\
            $seq_center \\
            $unaligned_opt \\
            $args \\
            2> ${prefix}.bowtie2.summary.log \\
            | samtools view -bS -F 4 -F 256 - > ${prefix}.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bowtie2: \$(bowtie2 --version | head -n 1 | sed 's/^bowtie2 version //')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        unaligned_opt = params.save_unaligned || params.contaminant_screening ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ''
        """
        INDEX=`find -L ./ -name "*.1.bt2" -o -name "*.1.bt2l" | head -n 1 | sed 's/\\.1.bt2l\$//' | sed 's/\\.1.bt2\$//'`
        bowtie2 \\
            -x \$INDEX \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            --threads ${task.cpus} \\
            $seq_center \\
            $unaligned_opt \\
            --no-mixed \\
            --no-discordant \\
            $args \\
            2> ${prefix}.bowtie2.summary.log \\
            | samtools view -bS -F 4 -F 8 -F 256 - > ${prefix}.bam

        if [ -f ${prefix}.unmapped.fastq.1.gz ]; then
            mv ${prefix}.unmapped.fastq.1.gz ${prefix}.unmapped_1.fastq.gz
        fi
        if [ -f ${prefix}.unmapped.fastq.2.gz ]; then
            mv ${prefix}.unmapped.fastq.2.gz ${prefix}.unmapped_2.fastq.gz
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bowtie2: \$(bowtie2 --version | head -n 1 | sed 's/^bowtie2 version //')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def unaligned = params.save_unaligned || params.contaminant_screening ? "echo '' | gzip > ${prefix}.unmapped_1.fastq.gz\necho '' | gzip > ${prefix}.unmapped_2.fastq.gz" : ''
    """
    ${unaligned}

    touch ${prefix}.bowtie2.summary.log
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | head -n 1 | sed 's/^bowtie2 version //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
