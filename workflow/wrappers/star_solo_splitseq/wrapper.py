from pathlib import Path
import shutil

from snakemake.shell import shell

from woldrnaseq.splitseq_merger import write_merged_splitseq_matrix


def make_comment_file(library_id, genome_dir, cofile):
    genome_dir = Path(genome_dir)
    with open(cofile, "wt") as outstream:
        if library_id is not None:
            outstream.write("@CO\tLIBID:{}\n".format(library_id))
        for comment_file in genome_dir.glob("*_bamCommentLines.txt"):
            with open(comment_file, "rt") as instream:
                for line in instream:
                    outstream.write(line)


genome_index = snakemake.input["genome_index"]
sequence_reads = snakemake.input["sequence_reads"]
barcode_reads = snakemake.input["barcode_reads"]

if isinstance(sequence_reads, list):
    sequence_reads = ",".join(sequence_reads)
if isinstance(barcode_reads, list):
    barcode_reads = ",".join(barcode_reads)

if sequence_reads.endswith(".gz"):
    read_files_command = "--readFilesCommand zcat"
elif sequence_reads.endswith(".bz2"):
    read_files_command = "--readFilesCommand bzcat"
elif sequence_reads.endswith(".xz"):
    read_files_command = "--readFilesCommand xzcat"
elif sequence_reads.endswith(".zstd"):
    read_files_command = "--readFilesCommand zstdcat"
else:
    read_files_command = ""

stranded = snakemake.params.stranded
gene_model = snakemake.params.gene_model
library_id = snakemake.params.library_id

aligned_bam = Path(snakemake.output.get("aligned_bam"))
output_dir = aligned_bam.parent
cofile = output_dir / "COfile.txt"
star_tmp = output_dir / "_STARtmp"
solo_dir = output_dir / "Solo.out"
filtered_dir = solo_dir / gene_model / "filtered"
raw_dir = solo_dir / gene_model / "raw"
raw_2bc_dir = solo_dir / gene_model / "raw_2bc"
sj_raw_dir = solo_dir / "SJ" / "raw"
sj_2bc_raw_dir = solo_dir / "SJ" / "raw_2bc"

barcode_23 = Path(snakemake.input.inclusion_list) / "CB23.txt"
barcode_1 = Path(snakemake.input.inclusion_list) / "CB1.txt"

if not cofile.exists():
    make_comment_file(library_id, genome_index, cofile)

if "mem_bytes" in snakemake.resources:
    mem_bytes = snakemake.resources.mem_bytes
    bam_sort_limit = "--limitBAMsortRAM {}".format(mem_bytes)
else:
    bam_sort_limit = ""

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "STAR --genomeDir {genome_index} \
           --readFilesIn {sequence_reads} {barcode_reads} {read_files_command} \
           --runThreadN {snakemake.threads} \
           --genomeLoad NoSharedMemory \
           --outFilterMultimapNmax 20 \
           --alignSJoverhangMin 8 \
           --alignSJDBoverhangMin 1 \
           --outFilterMismatchNmax 999 \
           --outFilterMismatchNoverReadLmax 0.04 \
           --alignIntronMin 20 \
           --alignIntronMax 1000000 \
           --alignMatesGapMax 1000000 \
           --outSAMheaderCommentFile {cofile} \
           --outSAMheaderHD @HD VN:1.4 SO:coordinate \
           --outSAMunmapped Within \
           --outFilterType BySJout \
           --outSAMattributes NH HI AS NM MD CB CR CY UB UR UY gx gn \
           --outSAMstrandField intronMotif \
           --outSAMtype BAM SortedByCoordinate \
           --sjdbScore 1 \
           --clipAdapterType CellRanger4 \
           --clip5pAdapterSeq AAGCAGTGGTATCAACGCAGAGTGAATGGG \
           --outFilterScoreMin 30 \
           --soloUMIdedup 1MM_CR \
           --soloUMIfiltering MultiGeneUMI_CR \
           --soloType CB_UMI_Complex \
           --soloCellFilter EmptyDrops_CR \
           --soloCBmatchWLtype EditDist_2 \
           --soloCBposition 0_10_0_17 0_48_0_55 0_78_0_85 \
           --soloUMIposition 0_0_0_9 \
           --soloCBwhitelist {barcode_23} {barcode_23} {barcode_1} \
           --soloStrand {stranded} \
           --soloFeatures {gene_model} SJ \
           --soloMultiMappers Unique EM \
           {bam_sort_limit} \
           --outTmpDir {star_tmp} \
           --outFileNamePrefix {output_dir}/"
)

if star_tmp.exists():
    shutil.rmtree(star_tmp)

if filtered_dir.exists():
    shutil.rmtree(filtered_dir)

if raw_2bc_dir.exists():
    print(
        "This shouldn't happen. Intermediate raw barcode directory exists",
        list(solo_dir.glob("*")),
    )
    shutil.rmtree(raw_2bc_dir)

if sj_2bc_raw_dir.exists():
    print(
        "This shouldn't happen. intermediate SJ 2 barcode directory exists",
        list(solo_dir.glob("*")),
    )
    shutil.rmtree(sj_2bc_raw_dir)

# Merge gene count barcodes
raw_dir.rename(raw_2bc_dir)

raw_2bc_unique_matrix = raw_2bc_dir / "matrix.mtx"
raw_2bc_em_matrix = raw_2bc_dir / "UniqueAndMult-EM.mtx"

raw_unique_matrix = raw_dir / "matrix.mtx"
raw_em_matrix = raw_dir / "UniqueAndMult-EM.mtx"

write_merged_splitseq_matrix(raw_2bc_unique_matrix, raw_dir, library_id=library_id)
write_merged_splitseq_matrix(raw_2bc_em_matrix, raw_dir, library_id=library_id)

# Merge the SJ barcodes
sj_raw_dir.rename(sj_2bc_raw_dir)
sj_2bc_matrix = sj_2bc_raw_dir / "matrix.mtx"
sj_raw_matrix = sj_raw_dir / "matrix.mtx"
write_merged_splitseq_matrix(sj_2bc_matrix, sj_raw_dir, library_id=library_id)

assert raw_unique_matrix.exists(), "{} does not exist".format(raw_unique_matrix)
assert raw_em_matrix.exists(), "{} does not exist".format(raw_em_matrix)
assert sj_raw_matrix.exists(), "{} does not exist".format(sj_raw_matrix)

shutil.rmtree(raw_2bc_dir)
shutil.rmtree(sj_2bc_raw_dir)
