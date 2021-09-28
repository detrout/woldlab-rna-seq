import netrc
import os
from pathlib import Path
import pprint
import requests

CONFIG_JSON = Path("config.json")
CONFIG_YAML = Path("config.yaml")

if CONFIG_JSON.exists():
    import json
    with open("config.json", "rt") as instream:
        config = json.load(instream)
elif CONFIG_YAML.exists():
    import yaml
    with open("config.yaml", "rt") as instream:
        config = yaml.load(instream)
else:
    raise RuntimeError("No recognized configuration file")


for key in ['genome_dir', 'allow_list', 'star_command', 'outdir']:
    config[key] = Path(config[key]).expanduser()

authdb = netrc.netrc()
username, _, password = authdb.hosts['www.encodeproject.org']
auth = (username, password)

def generate_read_argument(config, read):
    argument = []
    lane = 1
    read_num = int(read[-1])
    for accession in config[read]:
        argument.append(
            "{accession}_S1_L00{lane}_R{read_num}_001.fastq.gz".format(
                accession=config["accession"],
                lane=lane,
                read_num=read_num,
            ))
        lane += 1
    return argument


rule ALL:
    input:
        #config['outdir'] / "Solo.out" / "Gene" / "raw" / "matrix.mtx"
        config['outdir'] / "Log.final.out"


rule get_fastq:
    output:
        temp("{params.accession}_S1_L00{lane}_R{read}_001.fastq.gz")
    threads: 1
    params:
        accession = config['accession']
    run:
        shell("curl -L -o {params.accession}_S1_L00{wildcards.lane}_R{wildcards.read}.fastq.gz https://{params.username}:{params.password}@www.encodeproject.org/files/{wildcards.accession}/@@download/{wildcards.accession}.fastq.gz")

rule star_solo:
    input:
        sequence_reads = generate_read_argument(config, "read2"),
        barcode_reads = generate_read_argument(config, "read1"),
        genome_index = config['genome_dir'],
        allow_list = config['allow_list'],
    params:
        star_command = config['star_command'],
        stranded = config['stranded'],
        gene_model = "GeneFull" if config['include_intron'] else "Gene",
        sequence_reads = ",".join(generate_read_argument(config, "read2")),
        barcode_reads = ",".join(generate_read_argument(config, "read1")),
        umi_length = int(config["umi_length"]),
    resources:
        mem_mb = config['mem_mb'],
        mem_bytes = config['mem_mb'] * (2 ** 20),
        disk_mb = config['disk_mb'],
    threads: 16
    output:
        aligned_bam = config['outdir'] / "Aligned.sortedByCoord.out.bam",
        log_final = config['outdir'] / "Log.final.out",
#        barcode_stats = config['outdir'] / "Solo.out" / "Barcodes.stats",
#        features_stats = config['outdir'] / "Solo.out" / "Gene" / "Features.stats",
        # config['outdir'] / "Solo.out" / "Gene" / "UMIperCellSorted",
        star_tmp = temp(directory(config['outdir'] / "_STARtmp")),
#        filtered_barcodes = config['outdir'] / "Solo.out" / "Gene" / "filtered" / "barcodes.tsv",
#        filtered_features = config['outdir'] / "Solo.out" / "Gene" / "filtered" / "features.tsv",
#        filtered_matrix = config['outdir'] / "Solo.out" / "Gene" / "filtered" / "matrix.mtx",
#        raw_barcodes = config['outdir'] / "Solo.out" / "Gene" / "raw" / "barcodes.tsv",
#        raw_features = config['outdir'] / "Solo.out" / "Gene" / "raw" / "features.tsv",
#        raw_matrix = config['outdir'] / "Solo.out" / "Gene" / "raw" / "matrix.mtx",
    shell:
#           --quantMode TranscriptomeSAM  \ # Turns on outputing transcript bam file
        "{params.star_command} --genomeDir {input.genome_index} \
           --readFilesIn {params.sequence_reads} {params.barcode_reads} \
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
           --outSAMheaderCommentFile COfile.txt \
           --outSAMheaderHD @HD VN:1.4 SO:coordinate \
           --outSAMunmapped Within \
           --outFilterType BySJout \
           --outSAMattributes NH HI AS NM MD \
           --outSAMstrandField intronMotif \
           --outSAMtype BAM SortedByCoordinate \
           --sjdbScore 1 \
           --soloType Droplet \
           --soloUMIlen {params.umi_length} \
           --soloCBwhitelist {input.allow_list} \
           --soloStrand {params.stranded} \
           --soloFeatures Gene GeneFull SJ  \
           --limitBAMsortRAM {resources.mem_bytes} \
           --outTmpDir {output.star_tmp} \
q           --outFileNamePrefix ./ \
"
# my best guess for what we want but needs star 2.7.8
#           --outSAMtype BAM SortedByCoordinate \
#           --outSAMattributes NH HI AS NM MD CB UB \

#           --soloType CB_UMI_Simple \
#           --soloCellFilter EmptyDrops_CR \
#           --soloCBwhitelist {input.allow_list} \
#           --soloStrand {params.stranded} \
#           --soloFeatures SJ {params.gene_model} \
#           --soloMultiMappers Unique PropUniq \
# --soloFeatures SJ {params.gene_model}


#rule to_h5ad:
#    input:
#        filtered_barcodes = config['outdir'] / "Solo.out" / "Gene" / "filtered" / "barcodes.tsv",
#        raw_barcodes = config['outdir'] / "Solo.out" / "Gene" / "raw" / "barcodes.tsv",
#        raw_features = config['outdir'] / "Solo.out" / "Gene" / "raw" / "features.tsv",
#        raw_matrix = config['outdir'] / "Solo.out" / "Gene" / "raw" / "matrix.mtx",
#    output:
#        config['outdir'] / (config['analysis_name'] + '.h5ad')
