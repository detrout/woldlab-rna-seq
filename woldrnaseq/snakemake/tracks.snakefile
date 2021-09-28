from woldrnaseq import models
from pathlib import Path


config = models.config_from_condor_initialdir()


def generate_bigwig_targets(config):
    targets = []
    bedgraph_suffixes = [".Unique.str1.out.bg", ".UniqueMulti.str1.out.bg"]
    # the order of bigwig suffixes, and of bedgraphs files need to match
    # Unstranded
    #   Unique.str1 -> _uniq.bw
    #   UniqueMulti.str1 -> all.bw
    # Stranded
    #   Unique.str1 -> _minusUniq.bw
    #   UniqueMulti.str1 -> _minusAll.bw
    #   Unique.str2 -> _plusUniq.bw
    #   UniqueMulti.str2 -> _plusAll.bw

    if config["stranded"].lower() == "unstranded":
        suffixes = ["_uniq.bw", "_all.bw"]
    else:
        suffixes = ["_minusUniq.bw", "_minusAll.bw", "_plusUniq.bw", "_plusAll.bw"]
        bedgraph_suffixes.extend([".Unique.str2.out.bg", ".UniqueMulti.str2.out.bg"])
    config["bigwig_suffixes"] = suffixes
    config["bedgraph_suffixes"] = bedgraph_suffixes
    config["chrom_info"] = Path(config["genome_dir"]) / config["genome_triple"] / "chrNameLength.txt"

    for suffix in suffixes:
        targets.append("{analysis_name}-{genome_triple}{suffix}".format(
            analysis_name=config["analysis_name"],
            genome_triple=config["genome_triple"],
            suffix=suffix))
    return targets


rule ALL:
    input:
        generate_bigwig_targets(config)


rule bam_to_bedgraph:
    input:
        expand("{analysis_name}-{genome_triple}_genome.bam",
               analysis_name=config["analysis_name"],
               genome_triple=config["genome_triple"])
    output:
        temp(expand(["Signal{signal}"], signal=config["bedgraph_suffixes"]))
    resources:
        mem_mb = 2000
    run:
        shell(expand("{star_dir}/STAR --runMode inputAlignmentsFromBAM \
           --inputBAMfile {{input}} \
           --outWigType bedGraph \
           --outWigStrand {stranded} \
           --outWigReferencesPrefix {reference_prefix}",
                     star_dir=config["star_dir"],
                     stranded=config["stranded"],
                     reference_prefix=config["reference_prefix"]))

rule sorted_bedgraph:
    input:
        "Signal.{name}.bg"
    output:
        temp("Sorted.{name}.bg")
    shell:
        expand("{ucsc_tools_dir}/bedSort {{input}} {{output}}",
               ucsc_tools_dir=config["ucsc_tools_dir"])

rule bedgraph_to_bigwig:
    input:
#rule bedgraph_to_bigwig:
#    input:
#        expand(["Sorted{signal}"], signal=config["bedgraph_suffixes"])
#    output:
#        mutliext("{}-{}".format(config["analysis_name"],config["genome_triple"]),
#                 *config["bigwig_suffixes"])
#    resources:
#        mem_mb = 4000
#    version: "0.1"
#    run:
#        for i in range(len(output)):
#            shell("{ucsc_tools_dir}/bedGraphToBigWig {input} {chrom_info} {output}".format(
#                ucsc_tools_dir=config["ucsc_tools_dir"],
#                input=input[i],
#                chrom_info=config["chrom_info"],
#                output=output[i]))
#

