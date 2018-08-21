from collections.abc import Mapping
import logging
import os

import pandas
from .models import (
    genome_name_from_library,
    load_gtf_cache,
)


logger = logging.getLogger(__name__)


class GTFCache(Mapping):
    """Map library IDs to the correct GTF annotation table
    """
    def __init__(self, libraries, genome_dir):
        """Initializat class
        :Parameters:
          - libraries: (DataFrame) Library Table
          - genome_dir: (directory) root path where indexes are stored
        """
        assert isinstance(libraries, pandas.DataFrame)
        self._libraries = libraries
        self._genome_dir = genome_dir
        self._gtf_cache = {}

    def __getitem__(self, key):
        """Return gtf cache for a library_id
        :Parameters:
          - key: library ID
        """
        row = self._libraries.loc[key]
        genome_name = genome_name_from_library(row)
        return self.get_gtf_cache(genome_name)

    def __iter__(self):
        """iterator of keys
        """
        for k in self._libraries.index:
            yield k

    def __len__(self):
        """Return number of records
        """
        return len(self._libraries.index)

    def get_gtf_cache(self, genome_name):
        """Look through list of libraries and attempt to load GTF caches
        """
        if genome_name not in self._gtf_cache:
            cache_pathname = self._get_gtf_cache_filename(genome_name)
            logging.debug('Searching for %s', cache_pathname)

            try:
                self._gtf_cache[genome_name] = load_gtf_cache(cache_pathname)
            except FileNotFoundError as e:
                logging.error('Unable to load gene cache %s', cache_pathname)

        return self._gtf_cache.get(genome_name)

    def _get_gtf_cache_filename(self, genome_name):
        """Return the expected h5 cache file name
        :Paramters:
          - genome_name: (string) like mm10-M4-male
        :Returns: Filename
        """
        if self._genome_dir is None:
            logger.error("genome_dir is not specified. Please configure")
            raise ValueError("genome_dir is not set")
        return os.path.join(self._genome_dir, genome_name, genome_name + '.h5')
        
def protein_coding_gene_ids(annotation):
    """Filter GTF just protein coding genes
    """
    entry_type = (annotation['type'] == 'gene')
    gene_type = (annotation['gene_type'] == 'protein_coding')
    return annotation[entry_type & gene_type]['gene_id']
