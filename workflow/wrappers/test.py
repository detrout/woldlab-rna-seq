
from pathlib import Path
import os
import subprocess
import tempfile
import shutil

# I copied this from the snakemake-wrappers test.py module and then hacked on it.
def run(wrapper, cmd, check_log=None):
    wrapper = Path(wrapper)
    origdir = os.getcwd()
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        dst = d / "main"
        dst.mkdir(parents=True, exist_ok=True)
        print("Destdir", dst)
        copy = lambda pth, src: shutil.copy(
            os.path.join(pth, src), os.path.join(dst, pth)
        )

        used_wrappers = []
        wrapper_file = "used_wrappers.yaml"
        if (wrapper / wrapper_file).exists():
            # is meta wrapper
            with open(wrapper / wrapper_file, "r") as wf:
                wf = yaml.load(wf, Loader=yaml.BaseLoader)
                used_wrappers = wf["wrappers"]
        else:
            used_wrappers.append(wrapper)

        for w in used_wrappers:
            success = False
            for ext in ("py", "R", "Rmd"):
                script = "wrapper." + ext
                if os.path.exists(os.path.join(w, script)):
                    os.makedirs(os.path.join(dst, w), exist_ok=True)
                    copy(w, script)
                    success = True
                    break
            assert success, "No wrapper script found for {}".format(w)
            environment = w / "environment.yaml"
            if environment.exists():
                copy(w, "environment.yaml")

        print("Used wrappers", used_wrappers)
        #if (DIFF_MASTER or DIFF_LAST_COMMIT) and not any(
        #    any(f.startswith(w) for f in DIFF_FILES)
        #    for w in chain(used_wrappers, [wrapper])
        #):
        #    raise Skipped("wrappers not modified")

        testdir = d / "test"
        # pkgdir = os.path.join(d, "pkgs")
        shutil.copytree(wrapper / "test", testdir)
        # prepare conda package dir
        # os.makedirs(pkgdir)
        # switch to test directory
        os.chdir(testdir)
        if os.path.exists(".snakemake"):
            shutil.rmtree(".snakemake")
        cmd = cmd + [
            "--wrapper-prefix",
            "file://{}/".format(d),
            "--conda-cleanup-pkgs",
            "--printshellcmds",
        ]

        #if CONTAINERIZED:
        #    # run snakemake in container
        #    cmd = [
        #        "sudo",
        #        "docker",
        #        "run",
        #        "-it",
        #        "-v",
        #        "{}:{}".format(os.getcwd(), "/workdir"),
        #        "snakemake/snakemake",
        #        " ".join(cmd),
        #    ]

        # env = dict(os.environ)
        # env["CONDA_PKGS_DIRS"] = pkgdir
        try:
            subprocess.check_call(cmd)
        except Exception as e:
            # go back to original directory
            os.chdir(origdir)
            logfiles = [
                os.path.join(d, f)
                for d, _, files in os.walk(os.path.join(testdir, "logs"))
                for f in files
            ]
            for path in logfiles:
                with open(path) as f:
                    msg = "###### Logfile: " + path + " ######"
                    print(msg, "\n")
                    print(f.read())
                    print("#" * len(msg))
            if check_log is not None:
                for f in logfiles:
                    check_log(open(f).read())
            else:
                raise e
        finally:
            # cleanup environments to save disk space
            subprocess.check_call(
                "for env in `conda env list | grep -P '\.snakemake/conda' | "
                "cut -f1 | tr -d ' '`; do conda env remove --prefix $env; done",
                shell=True,
            )
            # go back to original directory
            os.chdir(origdir)


def test_get_fastq():
    run(
        "getfastqs",
        ["snakemake", "--cores", "1", "--use-conda"]
    )


def test_star_solo_splitseq():
    run(
        "star_solo_splitseq",
        ["snakemake", "--cores", "1", "--use-conda"]
    )
