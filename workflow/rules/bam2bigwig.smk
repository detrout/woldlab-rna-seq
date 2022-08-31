import sys
import pysam
from woldrnaseq.snakeutils import (
    get_library_stranded,
    get_library_by_analysis_dir
)

def get_star_stranded(config, libraries, analysis_dir):
    stranded = get_library_stranded(config, libraries, analysis_dir)
    if stranded.lower() == "unstranded":
        return "Unstranded"
    else:
        return "Stranded"

def get_reference_prefix(config, libraries, analysis_dir):
    if "reference_prefix" not in libraries.columns:
        return "chr"
    else:
        return get_library_by_analysis_dir(libraries, analysis_dir).reference_prefix


rule alignments_from_bam_stranded:
    input:
        genome_bam = "{analysis_dir}/{analysis_name}-{genome_name}_genome.bam",
    output:
        allminus = temp("{analysis_dir}/{analysis_name}-{genome_name}_minusAll.out.bg"),
        uniqminus = temp("{analysis_dir}/{analysis_name}-{genome_name}_minusUniq.out.bg"),
        allplus = temp("{analysis_dir}/{analysis_name}-{genome_name}_plusAll.out.bg"),
        uniqplus = temp("{analysis_dir}/{analysis_name}-{genome_name}_plusUniq.out.bg"),
    params:
        stranded = lambda wildcards: get_star_stranded(
            config, libraries, wildcards.analysis_dir),
        reference_prefix = lambda wildcards: get_reference_prefix(
            config, libraries, wildcards.analysis_dir),
    singularity:
        config["star_container"]
    shell:
        """STAR --runMode inputAlignmentsFromBAM \
           --inputBAMfile {input.genome_bam} \
           --outWigType bedGraph \
           --outWigReferencesPrefix {params.reference_prefix} \
           --outWigStrand {params.stranded} \
	   --outFileNamePrefix {wildcards.analysis_dir}/ ;
        mv -v {wildcards.analysis_dir}/Signal.Unique.str1.out.bg {output.uniqminus} ;
        mv -v {wildcards.analysis_dir}/Signal.UniqueMultiple.str1.out.bg {output.allminus}
        mv -v {wildcards.analysis_dir}/Signal.Unique.str2.out.bg {output.uniqplus} ;
        mv -v {wildcards.analysis_dir}/Signal.UniqueMultiple.str2.out.bg {output.allplus}
        """

rule alignments_from_bam_unstranded:
    input:
        genome_bam = "{analysis_dir}/{analysis_name}-{genome_name}_genome.bam",
        genome_bai = "{analysis_dir}/{analysis_name}-{genome_name}_genome.bam.bai",
    output:
        allreads = temp("{analysis_dir}/{analysis_name}-{genome_name}_all.out.bg"),
        uniqreads = temp("{analysis_dir}/{analysis_name}-{genome_name}_uniq.out.bg"),
    params:
        stranded = lambda wildcards: get_star_stranded(
            config, libraries, wildcards.analysis_dir),
    singularity:
        config["star_container"]
    shell:
        """STAR --runMode inputAlignmentsFromBAM \
           --inputBAMfile {input.genome_bam} \
           --outWigType bedGraph \
           --outWigReferencesPrefix {params.reference_prefix} \
           --outWigStrand {params.stranded} ;
        mv {analysis_dir}/Signal.Unique.str1.out.bg {output.uniqreads} ;
        mv {analysis_dir}/Signal.UniqueMultiple.str1.out.bg {output.allreads}
        """
        
rule sort_bedgraph:
    input:
        "{analysis_dir}/{analysis_name}-{genome_name}_{strand}.out.bg"
    output:
        temp("{analysis_dir}/{analysis_name}-{genome_name}_{strand}.sorted.bg")
    singularity:
        config["ucsc_tools_container"]
    shell:
        "bedSort {input} {output}"


def bigwig_names(config, libraries, analysis_dir):
    stranded = get_library_stranded(config, libraries, analysis_dir)
    if stranded.lower() != "unstranded":
        return [
            "{analysis_dir}/{analysis_name}-{genome_name}_all.bw",
            "{analysis_dir}/{analysis_name}-{genome_name}_uniq.bw"
        ]
    else:
        return [
            "{analysis_dir}/{analysis_name}-{genome_name}_minusAll.bw",
            "{analysis_dir}/{analysis_name}-{genome_name}_minusUniq.bw",
            "{analysis_dir}/{analysis_name}-{genome_name}_plusAll.bw"
            "{analysis_dir}/{analysis_name}-{genome_name}_plusUniq.bw"
        ]


rule chrominfo:
    input:
        genome_bam = "{analysis_dir}/{analysis_name}-{genome_name}_genome.bam",
        genome_bai = "{analysis_dir}/{analysis_name}-{genome_name}_genome.bam.bai",
    output:
        temp("{analysis_dir}/{analysis_name}-{genome_name}_genome.chrominfo"),
    run:
        with pysam.AlignmentFile(input.genome_bam, 'rb') as alignment, open(output[0], 'wt') as chrom_info:
            for row in alignment.header['SQ']:
                name = row['SN']
                length = row['LN']
                chrom_info.write(name)
                chrom_info.write('\t')
                chrom_info.write(str(length))
                chrom_info.write(os.linesep)

        return output[0]
        
        
rule bedgraph2bigwig:
    input:
        chrominfo = rules.chrominfo.output,
        bedgraph = "{analysis_dir}/{analysis_name}-{genome_name}_{strand}.sorted.bg",
    output:
        bigwig = "{analysis_dir}/{analysis_name}-{genome_name}_{strand}.bw"
    singularity:
        config["ucsc_tools_container"]
    shell:
        """
        bedGraphToBigWig {input.bedgraph} {input.chrominfo} {output.bigwig}
        """
    
