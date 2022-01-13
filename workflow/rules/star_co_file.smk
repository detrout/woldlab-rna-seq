import os


def make_comment_file(lib_id, genome_dir, output_dir):
    pathname = os.path.join(output_dir, 'COfile.txt')
    with open(pathname, 'wt') as outstream:
        if lib_id is not None:
            outstream.write('@CO\tLIBID:{}\n'.format(lib_id))
        for comment_file in glob(os.path.join(genome_dir, '*_bamCommentLines.txt')):
            with open(comment_file, 'rt') as instream:
                for line in instream:
                    outstream.write(line)


rule generate_co_file:
    input:
        genome_index = config['genome_dir'],
    output:
        temp("COfile.txt"),
    params:
        library_accession = config['library_accession'],
    singularity:
        config['star_container']
    resources:
        mem_mb = 100
    threads: 1
    shell:
        """
        echo -e "@CO\tLIBID:{params.library_accession}" > {output}
        cat {input.genome_index}/*_bamCommentLines.txt >> {output}
        """
