from setuptools import setup, find_packages
from woldrnaseq.version import get_git_version

setup(
    name='long-rna-seq-condor',
    version=get_git_version(),
    packages = find_packages(),
    package_data={
        'woldrnaseq': ['RELEASE-VERSION'],
        'woldrnaseq.templates': ['*.dagman', '*.html']
    },
    entry_points={
        'console_scripts': [
            'madqc = woldrnaseq.madqc:main',
            'make_dag = woldrnaseq.make_dag:main',
            'makersemcsv = woldrnaseq.makersemcsv:main',
            'qcreport = woldrnaseq.report:main',
        ],
    },
    install_requires=[
        'bokeh>=0.9.3',
        'numpy>=1.10',
        'pandas>=0.17',
        'matplotlib>=1.4',
        'jinja2>=2.8',
        'scipy>=0.17',
    ],
    author='Diane Trout',
    author_email='diane@caltech.edu',
    description='Implementation of "ENCODE long rna-seq pipeline" using condor',
    zip_safe=False,
    test_suite='woldrnaseq.tests',
)
