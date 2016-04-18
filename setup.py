from setuptools import setup, find_packages

setup(
    name='long-rna-seq-condor',
    version='0.1',
    packages = find_packages(),
    entry_points={
        'console_scripts': [
            'madqc = woldrnaseq.madqc.py:main',
            'make_dag = woldrnaseq.make_dag.py:main',
            'makersemcsv = woldrnaseq.makersemcsv.py:main'],
    },
    install_requires=[
        'bokeh>=0.9.3',
        'numpy>=1.10',
        'pandas>=0.17,<0.18',
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
