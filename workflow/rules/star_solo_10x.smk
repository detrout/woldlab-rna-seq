import gzip
import requests
import shutil
from urllib.parse import urlparse
from encoded_client.encoded import ENCODED

DEFAULT_MEM_MB = 1000


def compute_inclusion_list_name(url):
    """The resulting allow list file need be uncompressed
    """
    parts = urlparse(url)
    filename = Path(parts.path).name
    return filename.replace(".gz", "")


rule download_inclusion_list:
    output:
        allow_file = temp(compute_inclusion_list_name(config['inclusion_list_url']))
    params:
        inclusion_list_url = config['inclusion_list_url']
    resources:
        mem_mb = DEFAULT_MEM_MB,
    run:
        server = ENCODED(get_submit_host(config))
        with server.get_response(params.inclusion_list_url, stream=True) as response:
            with open(output.allow_file, "wb") as outstream:
                if (params.inclusion_list_url.endswith('.gz') or
                   response.headers.get("Content-Encoding") == "gzip"):
                    instream = gzip.GzipFile(fileobj=response.raw)
                else:
                    instream = response.raw
                shutil.copyfileobj(instream, outstream)


def generate_star_read_argument(config, read):
    """Convert accessions into list of {accession}_R{read_num}.fastq.gz

    this makes formatting for the star paired readFilesIn argument a bit
    simpler.
    """
    argument = []
    read_num = int(read[-1])
    for accession in config[read]:
        argument.append(
            "{accession}_R{read_num}.fastq.gz".format(
                accession=accession,
                read_num=read_num,
            ))
    return argument


rule star_solo_10x:
    input:
        sequence_reads = generate_star_read_argument(config, "read2"),
        barcode_reads = generate_star_read_argument(config, "read1"),
        genome_index = config['genome_dir'],
        inclusion_list = compute_inclusion_list_name(config['inclusion_list_url']),
        cofile = "COfile.txt",
    params:
        stranded = config['stranded'],
        gene_model = get_gene_model(),
        sequence_reads = ",".join(generate_star_read_argument(config, "read2")),
        barcode_reads = ",".join(generate_star_read_argument(config, "read1")),
        umi_length = int(config["umi_length"]),
        cb_length = int(config.get("cb_length", DEFAULT_10X_CB_LENGTH)),
        star_tmp = temp(directory("_STARtmp")),
    resources:
        mem_mb = config['mem_mb'],
        mem_bytes = config['mem_mb'] * (2 ** 20),
        disk_mb = config['disk_mb'],
    threads: 16
    log: "star_solo_10x.out"
    output:
        aligned_bam = "Aligned.sortedByCoord.out.bam",
        log_final = "Log.final.out",
        log_progress = "Log.progress.out",
        log_out = "Log.out",
        splice_junctions = "SJ.out.tab",
        barcode_stats = SOLO_ROOT / "Barcodes.stats",
        features_stats = SOLO_ROOT / get_gene_model() / "Features.stats",
        umis = SOLO_ROOT / get_gene_model() / "UMIperCellSorted.txt",
        gene_summary = SOLO_ROOT / get_gene_model() / "Summary.csv",
        filtered_barcodes = temp(SOLO_ROOT / get_gene_model() / "filtered" / "barcodes.tsv"),
        filtered_features = temp(SOLO_ROOT / get_gene_model() / "filtered" / "features.tsv"),
        filtered_uniq_matrix = temp(SOLO_ROOT / get_gene_model() / "filtered" / "matrix.mtx"),
        raw_barcodes = temp(SOLO_ROOT / get_gene_model() / "raw" / "barcodes.tsv"),
        raw_features = temp(SOLO_ROOT / get_gene_model() / "raw" / "features.tsv"),
        raw_unique_matrix = temp(SOLO_ROOT / get_gene_model() / "raw" / "matrix.mtx"),
        raw_em_matrix = temp(SOLO_ROOT / get_gene_model() / "raw" / "UniqueAndMult-EM.mtx"),
        sj_feature_stats = temp(SOLO_ROOT / "SJ" / "Features.stats"),
        sj_summary = SOLO_ROOT / "SJ" / "Summary.csv",
        sj_barcodes = temp(SOLO_ROOT / "SJ" / "raw" / "barcodes.tsv"),
        sj_features = temp(SOLO_ROOT / "SJ" / "raw" / "features.tsv"),
        sj_matrix = temp(SOLO_ROOT / "SJ" / "raw" / "matrix.mtx"),
    singularity:
        config['star_container']
    shell:
        "STAR --genomeDir {input.genome_index} \
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
           --outSAMattributes NH HI AS NM MD CB CR CY UB UR UY gx gn \
           --outSAMstrandField intronMotif \
           --outSAMtype BAM SortedByCoordinate \
           --sjdbScore 1 \
           --clipAdapterType CellRanger4 \
           --outFilterScoreMin 30 \
           --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
           --soloUMIdedup 1MM_CR \
           --soloUMIfiltering MultiGeneUMI_CR \
           --soloType CB_UMI_Simple \
           --soloCellFilter EmptyDrops_CR \
           --soloUMIlen {params.umi_length} \
           --soloCBlen {params.cb_length} \
           --soloBarcodeReadLength 0 \
           --soloCBwhitelist '{input.inclusion_list}' \
           --soloStrand {params.stranded} \
           --soloFeatures {params.gene_model} SJ \
           --soloMultiMappers Unique EM \
           --limitBAMsortRAM {resources.mem_bytes} \
           --outTmpDir {params.star_tmp} \
           --outFileNamePrefix ./ 2>&1 >> {log} ; \
           rm -rf {params.star_tmp} \
"
