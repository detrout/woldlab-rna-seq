import shutil

# still need to make CO file
# 
#+PreCmd="pre_star"
#+PreArguments="--genome-dir $(GENOME_DIR) --library-id $(library_id) --output-dir $(CURDIR)"

from woldrnaseq.snakeutils import (
    get_genome_dir,
    get_library_by_analysis_dir,
)


def generate_star_read_argument(libraries, read_name, analysis_dir):
    assert read_name in ["read_1", "read_2"], "Invalid read name {}".format(read_name)
    library = libraries.set_index("analysis_dir").loc[analysis_dir]
    fastqs = []
    if read_name in libraries.columns:
        return library[read_name]
    else:
        return []


def get_paired_args(libraries, analysis_dir):
    library = libraries.set_index("analysis_dir").loc[analysis_dir]
    if "read_2" in libraries.columns:
        return ""
    else:
        return "--outSAMstrandField intronMotif"


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
        genome_index = lambda wildcards: get_genome_dir(
            config, libraries, wildcards.analysis_dir),
    output:
        temp("{analysis_dir}/COfile.txt"),
    params:
        library_accession = lambda wildcards: get_library_by_analysis_dir(
            libraries, wildcards.analysis_dir).name,
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
    
rule encode_bulk_star:
    input:
        read1 = lambda wildcards: generate_star_read_argument(libraries, "read_1", wildcards.analysis_dir),
        read2 = lambda wildcards: generate_star_read_argument(libraries, "read_2", wildcards.analysis_dir),
        genome_index = lambda wildcards: get_genome_dir(
            config, libraries, wildcards.analysis_dir),
        cofile = "{analysis_dir}/COfile.txt",
    output:
        genome_bam = temp("{analysis_dir}/Aligned.sortedByCoord.out.bam"),
        transcriptome_bam = temp("{analysis_dir}/Aligned.toTranscriptome.out.bam"),
        log_out = "{analysis_dir}/Log.out",
        log_final = "{analysis_dir}/Log.final.out",
        log_progress = temp("{analysis_dir}/Log.progress.out"),
        splice_junctions = "{analysis_dir}/SJ.out.tab",
    params:
        read1_formatted = lambda wildcards: ",".join(generate_star_read_argument(libraries, "read_1", wildcards.analysis_dir)),
        read2_formatted = lambda wildcards: ",".join(generate_star_read_argument(libraries, "read_2", wildcards.analysis_dir)),
        stranded = libraries["stranded"],
        star_tmp = temp(directory("{analysis_dir}/_STARtmp")),
        paired_args = lambda wildcards: get_paired_args(libraries, wildcards.analysis_dir)
    wildcard_constraints:
        analysis_name = r"[^/]+"
    resources:
#        mem_mb = config['mem_mb'],
#        mem_bytes = config['mem_mb'] * (2 ** 20),
#        disk_mb = config['disk_mb'],
        mem_mb = 60000,
        mem_bytes = 60000 * (2 ** 20),
        disk_mb = 40000,
    threads: 16
    log: "{analysis_dir}/align-star.out"
    singularity:
        config['star_container']
    shell:
        """
        if [ -e {params.star_tmp} ]; then rm -rf {params.star_tmp} ; fi; \
        STAR --genomeDir {input.genome_index} \
           --readFilesIn {params.read1_formatted} {params.read2_formatted} \
           --readFilesCommand zcat \
	   --runThreadN {threads} \
	   --genomeLoad NoSharedMemory \
	   --outFilterMultimapNmax 20 \
	   --alignSJoverhangMin 8 \
	   --alignSJDBoverhangMin 1 \
	   --outFilterMismatchNmax 999 \
	   --outFilterMismatchNoverReadLmax 0.04 \
	   --alignIntronMin 20 \
	   --alignIntronMax 1000000 \
	   --alignMatesGapMax 1000000 \
	   --outSAMheaderCommentFile {input.cofile} \
	   --outSAMheaderHD @HD VN:1.4 SO:coordinate \
	   --outSAMunmapped Within \
	   --outFilterType BySJout \
	   --outSAMattributes NH HI AS NM MD \
	   --outSAMtype BAM SortedByCoordinate \
	   --quantMode TranscriptomeSAM GeneCounts \
	   --sjdbScore 1 \
           --outTmpDir {params.star_tmp}/ \
	   --outFileNamePrefix {wildcards.analysis_dir}/ \
	   --limitBAMsortRAM {resources.mem_bytes} \
           {params.paired_args};
           ls {wildcards.analysis_dir}/*.bam;
"""


rule rename_genome_bam:
    input:
        "{analysis_dir}/Aligned.sortedByCoord.out.bam",
    output:
        "{analysis_dir}/{analysis_name}-{genome_name}_genome.bam"
    threads: 1
    resources:
        mem_mb = 100
    run:
        source = Path(input[0])
        destination = Path(output[0])
        shutil.copy(source, destination)
        #if source.exists():
        #    # Do I need to delete destination if it exists?
        #    os.rename(source, destination)
        #elif destination.exists():
        #    # this step was already done
        #    print("{} already exists.".format(destination))
        #else:
        #    raise FileNotFoundError("Neither {} or {} exists".format(source, destination))


rule sort_transcriptome_bam:
    input:
        "{analysis_dir}/Aligned.toTranscriptome.out.bam"
    output:
        "{analysis_dir}/{analysis_name}-{genome_name}_anno.bam"
    singularity:
        config["rsem_container"]
    shell:
        """
        rsem-sam-validator {input} | grep "is valid\!"
        if [ $? -eq 0 ]; then 
           mv {input} {output}
        else
           convert-sam-for-rsem {input} {output}
        fi
        """
       
