from setuptools import setup, find_packages
from woldrnaseq.version import get_git_version

setup(
    name='long-rna-seq-condor',
    version=get_git_version(),
    packages=find_packages(),
    package_data={
        'woldrnaseq': ['RELEASE-VERSION'],
        'woldrnaseq.templates': ['*.dagman', '*.html']
    },
    entry_points={
        'console_scripts': [
            'madqc = woldrnaseq.madqc:main',
            'make_dag = woldrnaseq.make_dag:main',
            'make_rsem_csv = woldrnaseq.makersemcsv:main',
            'make_star_csv = woldrnaseq.makestarcsv:main',
            'qc_report = woldrnaseq.report:main',
            'make_trackhub = woldrnaseq.make_trackhub:main',
            'plot_genes_detected = woldrnaseq.plot_genes_detected:main',
            'merge_encode_annotations = woldrnaseq.merge_encode_annotations:main',
        ],
    },
    install_requires=[
        'bokeh>=0.12.10,<1.0',
        'jinja2>=2.8',
        'matplotlib>=3.0',
        'numpy>=1.16',
        'pandas>=0.23',
        'scipy>=1.1',
        'tables>=3.4',
        'xopen',
    ],
    author='Diane Trout',
    author_email='diane@caltech.edu',
    description='Implementation of "ENCODE long rna-seq pipeline" using condor',
    zip_safe=False,
    test_suite='woldrnaseq.tests',
)
