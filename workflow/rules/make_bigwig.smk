import sys

rule make_bigwig:
    input:
        pass
    output:
        pass
    singularity:
        config["uscs_tools"]
    rule:
        # this doesn't work because I'd need the make_bigwig script to be in the container.
        """{sys.executable} -m """
