

rule sort_genome_annotation:
    input:
        "{annotation}"
    output:
        "{annotation}.sorted.gz"
    run:
        '(grep ^"#" {input}; grep -v ^"#" {input} | sort -k1,1 -k4,4n) | bgzip > {output}'

rule index_genome_annotation:
    input:
        "{annotation}.sorted.gz"
    output:
        "{annotation}.sorted.gz.tbi"
    run:
        "tabix -p gff {input}"
        
