process BOWTIE2_BUILD {
    tag "${fasta}"
    label 'process_high'
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.5.1--py39hfa7eb30_4' :
        'quay.io/biocontainers/bowtie2:2.5.1--py39hfa7eb30_4' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bowtie2") , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bowtie2
    bowtie2-build \\
        -f \\
        ${args} \\
        ${fasta} \\
        bowtie2/${fasta.baseName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | head -n 1 | sed 's/^bowtie2 version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p bowtie2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | head -n 1 | sed 's/^bowtie2 version //')
    END_VERSIONS
    """
}
