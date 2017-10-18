import os
from unittest import TestCase, main
from pkg_resources import resource_filename

from woldrnaseq.models import load_library_tables
from woldrnaseq import make_tracks

class TestMakeTracks(TestCase):
    def setUp(self):
        mm10tsv = resource_filename(__name__, 'library-mm10-se.tsv')
        hg38tsv = resource_filename(__name__, 'library-hg38-se.tsv')
        self.root = os.path.abspath(os.path.join(os.path.dirname(mm10tsv)))
        self.mm10 = load_library_tables([mm10tsv])
        self.hg38 = load_library_tables([hg38tsv])


    def test_make_bigwig_custom_track_subdir(self):
        row = self.mm10.loc['12304']
        web_root = 'https://woldlab.caltech.edu/~user/project/'
        track = make_tracks.make_bigwig_custom_tracks(row, web_root, self.root)
        print(track)
        
    def test_make_bigwig_track_name_subdir(self):
        """Can we find bigwigs in analysis subdirectories
        """
        row = self.mm10.loc['12304']
        track = make_tracks.make_bigwig_track_name(row, 'uniq', self.root)
        self.assertEqual(track, 'lib_12304/lib_12304-mm10-M4-female_uniq.bw')
        track = make_tracks.make_bigwig_track_name(row, 'all', self.root)
        self.assertEqual(track, 'lib_12304/lib_12304-mm10-M4-female_all.bw')

    def test_make_bigwig_track_name_flat(self):
        """Can we find bigwigs when they've all been moved into a single directory
        """
        row = self.hg38.loc['12310']
        track = make_tracks.make_bigwig_track_name(row, 'uniq', self.root)
        self.assertEqual(track, 'lib_12310-hg38-V38-male_uniq.bw')
        track = make_tracks.make_bigwig_track_name(row, 'all', self.root)
        self.assertEqual(track, 'lib_12310-hg38-V38-male_all.bw')

    def test_make_bam_custom_track_subdir(self):
        row = self.mm10.loc['12304']
        web_root = 'https://woldlab.caltech.edu/~user/project/'
        track = make_tracks.make_bam_custom_track(row, web_root, self.root)
        print(track)
        
    def test_make_bam_track_name_subdir(self):
        """Can we find bams when still in subdirectories

        This is made harder by a later version of long-rna-seq-condor
        renaming the default STAR filename from Aligned.sortedByCoord.out.bam
        to <analysis_name>-<genome_triplet>_genome.bam
        """
        row = self.mm10.loc['12304']
        track = make_tracks.make_bam_track_name(row, self.root)
        self.assertEqual(track, 'lib_12304/lib_12304-mm10-M4-female_genome.bam')

        row = self.mm10.loc['12307']
        track = make_tracks.make_bam_track_name(row, self.root)
        self.assertEqual(track, 'lib_12307/Aligned.sortedByCoord.out.bam')

    def test_make_bam_track_name_flat(self):
        """Can we find bams when still in a single directory

        This has to use the newer naming convention as
        Aligned.sortedByCoord.out.bam would collide.
        """
        row = self.hg38.loc['12310']
        track = make_tracks.make_bam_track_name(row, self.root)
        self.assertEqual(track, 'lib_12310-hg38-V38-male_genome.bam')

    def test_return_subpath(self):
        """Properly strip out the analysis_root from our returned path
        """
        path = '/home/user/analysis/library/uniq.bw'
        result = 'library/uniq.bw'

        for root in ['/home/user/analysis/',
                     '/home/user/analysis']:
            self.assertEqual(make_tracks.return_subpath(path, root),
                             'library/uniq.bw')
if __name__ == '__main__':
    main()
