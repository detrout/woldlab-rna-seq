import pandas
from unittest import TestCase
from unittest.mock import patch
import pytest

pytest.importorskip("htsworkflow")

from woldrnaseq.downloader import (
    FastqFragment,
    get_library_flowcells,
    guess_runfolder_type,
    parse_apache_dirindex,
    search_for_runfolders,
)

import lxml.html

runfolder_index = lxml.html.fromstring(
    """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">
<html>
 <head>
  <title>Index of /runfolders/host02</title>
 </head>
 <body>
<h1>Index of /runfolders/host02</h1>
<table><tr><th><img src="/icons/blank.gif" alt="[ICO]"></th><th><a href="?C=N;O=D">Name</a></th><th><a href="?C=M;O=A">Last modified</a></th><th><a href="?C=S;O=A">Size</a></th><th><a href="?C=D;O=A">Description</a></th></tr><tr><th colspan="5"><hr></th></tr>
<tr><td valign="top"><img src="/icons/back.gif" alt="[DIR]"></td><td><a href="/runfolders/">Parent Directory</a></td><td>&nbsp;</td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/text.gif" alt="[TXT]"></td><td><a href="2018_06_11_samplesheet_fixed.csv">2018_06_11_samplesheet_fixed.csv</a></td><td align="right">11-Jun-2018 12:11  </td><td align="right">5.6K</td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="190108_SN787_0947_AHBBB222CX2/">190108_SN787_0947_AHBBB222CX2/</a></td><td align="right">11-Jan-2019 12:20  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="190114_SN787_0952_BHAAA111CX2/">190114_SN787_0952_BHAAA111CX2/</a></td><td align="right">15-Jan-2019 08:15  </td><td align="right">  - </td><td>&nbsp;</td></tr>
</table>
<address>Apache/2.2.0 Server at host.example.edu Port 443</address>
</body></html>
"""
)

hiseq_runfolder = lxml.html.fromstring(
    """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">
<html>
 <head>
  <title>Index of /runfolders/host02/191106_SN787_1115_AHCCDDDCX3</title>
 </head>
 <body>
<h1>Index of /runfolders/host02/191106_SN787_1115_AHCCDDDCX3</h1>
<table><tr><th><img src="/icons/blank.gif" alt="[ICO]"></th><th><a href="?C=N;O=D">Name</a></th><th><a href="?C=M;O=A">Last modified</a></th><th><a href="?C=S;O=A">Size</a></th><th><a href="?C=D;O=A">Description</a></th></tr><tr><th colspan="5"><hr></th></tr>
<tr><td valign="top"><img src="/icons/back.gif" alt="[DIR]"></td><td><a href="/runfolders/volvox02/">Parent Directory</a></td><td>&nbsp;</td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="BarcodeImages/">BarcodeImages/</a></td><td align="right">06-Nov-2019 21:22  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/text.gif" alt="[TXT]"></td><td><a href="Basecalling_Netcopy_complete.txt">Basecalling_Netcopy_complete.txt</a></td><td align="right">06-Nov-2019 21:49  </td><td align="right"> 45 </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/text.gif" alt="[TXT]"></td><td><a href="Basecalling_Netcopy_complete_Read1.txt">Basecalling_Netcopy_complete_Read1.txt</a></td><td align="right">06-Nov-2019 21:22  </td><td align="right"> 45 </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/text.gif" alt="[TXT]"></td><td><a href="Basecalling_Netcopy_complete_Read2.txt">Basecalling_Netcopy_complete_Read2.txt</a></td><td align="right">06-Nov-2019 21:49  </td><td align="right"> 45 </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="Config/">Config/</a></td><td align="right">06-Nov-2019 21:22  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="Data/">Data/</a></td><td align="right">06-Nov-2019 21:22  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/text.gif" alt="[TXT]"></td><td><a href="First_Base_Report.htm">First_Base_Report.htm</a></td><td align="right">06-Nov-2019 14:46  </td><td align="right">1.6K</td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/text.gif" alt="[TXT]"></td><td><a href="ImageAnalysis_Netcopy_complete.txt">ImageAnalysis_Netcopy_complete.txt</a></td><td align="right">06-Nov-2019 21:49  </td><td align="right"> 45 </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/text.gif" alt="[TXT]"></td><td><a href="ImageAnalysis_Netcopy_complete_Read1.txt">ImageAnalysis_Netcopy_complete_Read1.txt</a></td><td align="right">06-Nov-2019 21:21  </td><td align="right"> 45 </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/text.gif" alt="[TXT]"></td><td><a href="ImageAnalysis_Netcopy_complete_Read2.txt">ImageAnalysis_Netcopy_complete_Read2.txt</a></td><td align="right">06-Nov-2019 21:49  </td><td align="right"> 45 </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="InterOp/">InterOp/</a></td><td align="right">06-Nov-2019 21:48  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="Logs/">Logs/</a></td><td align="right">06-Nov-2019 21:49  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="PeriodicSaveRates/">PeriodicSaveRates/</a></td><td align="right">06-Nov-2019 21:22  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/text.gif" alt="[TXT]"></td><td><a href="RTAComplete.txt">RTAComplete.txt</a></td><td align="right">06-Nov-2019 21:49  </td><td align="right"> 45 </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="Recipe/">Recipe/</a></td><td align="right">06-Nov-2019 21:22  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/unknown.gif" alt="[   ]"></td><td><a href="RunInfo.xml">RunInfo.xml</a></td><td align="right">06-Nov-2019 14:51  </td><td align="right">571 </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/text.gif" alt="[TXT]"></td><td><a href="SampleSheet110619.csv">SampleSheet110619.csv</a></td><td align="right">06-Nov-2019 13:53  </td><td align="right">3.4K</td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="Thumbnail_Images/">Thumbnail_Images/</a></td><td align="right">06-Nov-2019 14:43  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="Unaligned/">Unaligned/</a></td><td align="right">07-Nov-2019 07:51  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/unknown.gif" alt="[   ]"></td><td><a href="runParameters.xml">runParameters.xml</a></td><td align="right">06-Nov-2019 14:51  </td><td align="right">5.5K</td><td>&nbsp;</td></tr>
<tr><th colspan="5"><hr></th></tr>
</table>
<address>Apache/2.2.0 Server at host.example.edu Port 443</address>
</body></html>
"""
)

nextseq_runfolder = lxml.html.fromstring(
    """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">
<html>
 <head>
  <title>Index of /runfolders/volvox02/nextseq/211011_VH00472_2_AAAWXYZHV</title>
 </head>
 <body>
<h1>Index of /runfolders/volvox02/nextseq/211011_VH00472_2_AAAWXYZHV</h1>
<table><tr><th><img src="/icons/blank.gif" alt="[ICO]"></th><th><a href="?C=N;O=D">Name</a></th><th><a href="?C=M;O=A">Last modified</a></th><th><a href="?C=S;O=A">Size</a></th><th><a href="?C=D;O=A">Description</a></th></tr><tr><th colspan="5"><hr></th></tr>
<tr><td valign="top"><img src="/icons/back.gif" alt="[DIR]"></td><td><a href="/runfolders/volvox02/nextseq/">Parent Directory</a></td><td>&nbsp;</td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="Analysis/">Analysis/</a></td><td align="right">14-Oct-2021 11:21  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="Autocenter/">Autocenter/</a></td><td align="right">11-Oct-2021 16:01  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="Autofocus/">Autofocus/</a></td><td align="right">13-Oct-2021 14:41  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="Config/">Config/</a></td><td align="right">11-Oct-2021 11:59  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/text.gif" alt="[TXT]"></td><td><a href="CopyComplete.txt">CopyComplete.txt</a></td><td align="right">13-Oct-2021 14:49  </td><td align="right">  0 </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="Data/">Data/</a></td><td align="right">11-Oct-2021 11:59  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="FocusModelGeneration/">FocusModelGeneration/</a></td><td align="right">13-Oct-2021 09:07  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="InstrumentAnalyticsLogs/">InstrumentAnalyticsLogs/</a></td><td align="right">13-Oct-2021 14:42  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="InterOp/">InterOp/</a></td><td align="right">13-Oct-2021 13:51  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="Logs/">Logs/</a></td><td align="right">13-Oct-2021 14:42  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="PrimaryAnalysisMetrics/">PrimaryAnalysisMetrics/</a></td><td align="right">13-Oct-2021 14:42  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/unknown.gif" alt="[   ]"></td><td><a href="RTA3.cfg">RTA3.cfg</a></td><td align="right">11-Oct-2021 11:59  </td><td align="right">4.8K</td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/text.gif" alt="[TXT]"></td><td><a href="RTAComplete.txt">RTAComplete.txt</a></td><td align="right">13-Oct-2021 13:51  </td><td align="right">  1 </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="Recipe/">Recipe/</a></td><td align="right">11-Oct-2021 11:59  </td><td align="right">  - </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/unknown.gif" alt="[   ]"></td><td><a href="RunCompletionStatus.xml">RunCompletionStatus.xml</a></td><td align="right">13-Oct-2021 14:42  </td><td align="right">344 </td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/unknown.gif" alt="[   ]"></td><td><a href="RunInfo.xml">RunInfo.xml</a></td><td align="right">11-Oct-2021 11:59  </td><td align="right">9.1K</td><td>&nbsp;</td></tr>
<tr><td valign="top"><img src="/icons/unknown.gif" alt="[   ]"></td><td><a href="RunParameters.xml">RunParameters.xml</a></td><td align="right">13-Oct-2021 14:41  </td><td align="right">2.8K</td><td>&nbsp;</td></tr>
<tr><th colspan="5"><hr></th></tr>
</table>
<address>Apache/2.2.0 (Debian) Server at host.example.org Port 443</address>
</body></html>
"""
)


class FakeHtswApi:
    def __init__(self, host, auth):
        self._host = host
        self._auth = auth

    def get_library(self, library_id):
        fake_data = [
            {
                "library_id": "23072",
                "library_name": "RUSH cDNA409",
                "lane_set": [
                    {
                        "status": "Unknown",
                        "paired_end": False,
                        "read_length": 100,
                        "flowcell": "HKHC7BCX3",
                    },
                ],
            },
            {
                "library_id": "23073",
                "library_name": "RUSH cDNA410",
                "lane_set": [
                    {
                        "status": "Unknown",
                        "paired_end": False,
                        "read_length": 100,
                        "flowcell": "HKHC7BCX3",
                    },
                ],
            },
            {
                "library_id": "23074",
                "library_name": "RUSH cDNA411",
                "lane_set": [
                    {
                        "status": "Unknown",
                        "paired_end": False,
                        "read_length": 100,
                        "flowcell": "HKHC7BCX3",
                    },
                ],
            },
            {
                "library_id": "23075",
                "library_name": "RUSH cDNA412",
                "lane_set": [
                    {
                        "status": "Failed",
                        "paired_end": False,
                        "read_length": 100,
                        "flowcell": "HKHCFBCX3",
                    },
                ],
            },
            {
                "library_id": "23196",
                "library_name": "gastrocnemius cDNA491",
                "lane_set": [
                    {
                        "status": "Unknown",
                        "paired_end": False,
                        "read_length": 100,
                        "flowcell": "HKGWTBCX3",
                    },
                    {
                        "status": "Good",
                        "paired_end": True,
                        "read_length": 150,
                        "flowcell": "AAAJ5NWHV",
                    },
                ],
            },
        ]

        for row in fake_data:
            if library_id == row["library_id"]:
                return row


def make_mini_library_row(library_id, paired=False):
    row = {
        "library_id": library_id,
        "analysis_name": library_id,
        "analysis_dir": library_id,
        "genome": "mm10",
        "annotation": "M21",
        "sex": "male",
        "stranded": "reverse",
        "read_1": ["{}/*_R1.fastq.gz".format(library_id)],
    }
    if paired:
        row["read_2"] = ["{}/*_R2.fastq.gz".format(library_id)]
    return row


class TestDownloader(TestCase):
    def test_make_short_fastq_name(self):
        ff = FastqFragment("22160_GAACGAAG_L002_R1_005.fastq.gz")
        self.assertEqual(ff.short_name, "22160_GAACGAAG_L002_R1.fastq.gz")

        ff = FastqFragment("22160_GAACGAAG_L002_R2_005.fastq.gz", merge_lanes=True)
        self.assertEqual(ff.short_name, "22160_GAACGAAG_R2.fastq.gz")

        ff = FastqFragment("22160_GAACGAAG_L002_R2_005.fastq.gz")
        self.assertEqual(ff.short_name, "22160_GAACGAAG_L002_R2.fastq.gz")

        ff = FastqFragment("23602_indexi7_257-i5_257_S1_R1_001.fastq.gz")
        self.assertEqual(ff.short_name, "23602_indexi7.257-i5.257_R1.fastq.gz")

        self.assertRaises(ValueError, FastqFragment, "22160_GAACGAAG_R1.fastq.gz")

        self.assertRaises(ValueError, FastqFragment, "Undetermined_S0_R1_001.fastq.gz")

    @patch("woldrnaseq.downloader.parse", return_value=runfolder_index)
    def test_parse_apache_dirindex_on_index_folder(self, parse):
        url = "http://example.org/runfolder/"
        runfolders = list(parse_apache_dirindex(url))
        assert len(runfolders) == 5
        assert runfolders[0].link_type == "header"
        assert runfolders[1].link_type == "parent"
        assert runfolders[2].link_type == "file"
        name = "190108_SN787_0947_AHBBB222CX2/"
        assert runfolders[3].link_type == "subdirectory"
        assert runfolders[3].name == name
        assert runfolders[3].url == url + name
        name = "190114_SN787_0952_BHAAA111CX2/"
        assert runfolders[4].link_type == "subdirectory"
        assert runfolders[4].name == name
        assert runfolders[4].url == url + name

    @patch("woldrnaseq.downloader.parse", return_value=runfolder_index)
    def test_search_for_runfolders(self, parse):
        url = "http://example.org/runfolders"
        flowcells = ["HBBB222CX2", "HAAA111CX2", "AAABCDEHV"]
        flowcell_urls = search_for_runfolders([url], flowcells)

        self.assertEqual(len(flowcell_urls), 2)
        self.assertEqual(
            flowcell_urls["HBBB222CX2"], url + "/190108_SN787_0947_AHBBB222CX2/"
        )
        self.assertEqual(
            flowcell_urls["HAAA111CX2"], url + "/190114_SN787_0952_BHAAA111CX2/"
        )

    @patch("woldrnaseq.downloader.parse", return_value=hiseq_runfolder)
    def test_guess_runfolder_type_hiseq(self, parse):
        url = "http://example/runfolder/191106_SN787_1115_AHCCDDDCX3"
        direntries = list(parse_apache_dirindex(url))
        folder_type = guess_runfolder_type(direntries)
        assert folder_type == "HiSeq"

    @patch("woldrnaseq.downloader.parse", return_value=nextseq_runfolder)
    def test_guess_runfolder_type_nextseq(self, parse):
        url = "http://example/runfolders/host02/191106_SN787_1115_AHCCDDDCX3"
        direntries = list(parse_apache_dirindex(url))
        folder_type = guess_runfolder_type(direntries)
        assert folder_type == "nextseq"

    def test_get_flowcells_one_flowcell(self):
        fake_api = FakeHtswApi("https://htsw.example.edu", (None, None))
        library_ids = ["23072", "23073"]
        libraries = pandas.DataFrame(
            [make_mini_library_row(x) for x in library_ids]
        ).set_index("library_id")

        flowcells = get_library_flowcells(fake_api, libraries)
        self.assertEqual(list(flowcells), ["HKHC7BCX3"])

    def test_get_flowcells_bad_flowcells(self):
        fake_api = FakeHtswApi("https://htsw.example.edu", (None, None))
        library_ids = ["23075"]
        libraries = pandas.DataFrame(
            [make_mini_library_row(x) for x in library_ids]
        ).set_index("library_id")

        flowcells = get_library_flowcells(fake_api, libraries)
        self.assertEqual(flowcells, set())

    def test_get_flowcells_two_flowcells(self):
        fake_api = FakeHtswApi("https://htsw.example.edu", (None, None))
        library_ids = ["23196"]
        libraries = pandas.DataFrame(
            [make_mini_library_row(x) for x in library_ids]
        ).set_index("library_id")

        flowcells = get_library_flowcells(fake_api, libraries)
        self.assertEqual(list(flowcells), ["HKGWTBCX3", "AAAJ5NWHV"])

    def test_get_flowcells_filtered_flowcells(self):
        fake_api = FakeHtswApi("https://htsw.example.edu", (None, None))
        library_ids = [
            "23196",
        ]
        libraries = pandas.DataFrame(
            [make_mini_library_row(x) for x in library_ids]
        ).set_index("library_id")

        flowcells = get_library_flowcells(
            fake_api, libraries, ["HKGWTBCX3", "BBBBBBBHV"]
        )
        self.assertEqual(list(flowcells), ["HKGWTBCX3"])

    def test_fastq_fragment_hiseq(self):
        url = "https://example.edu/flowcells/210416_SN787_1327_BHKH2HBCX3/23227_GGACATCA_L001_R1_001.fastq.gz"
        hiseq_frag = FastqFragment(url)
        self.assertEqual(hiseq_frag.url, url)
        self.assertEqual(hiseq_frag.merge_lanes, False)
        self.assertEqual(hiseq_frag.name, "23227_GGACATCA_L001_R1_001.fastq.gz")
        self.assertEqual(hiseq_frag.library_id, "23227")
        self.assertEqual(hiseq_frag.index, "GGACATCA")
        self.assertEqual(hiseq_frag.lane, "L001")
        self.assertEqual(hiseq_frag.read, "R1")
        self.assertEqual(hiseq_frag.chunk, "001")
        self.assertEqual(hiseq_frag.key, ("23227", "GGACATCA", "L001", "R1"))
        self.assertEqual(hiseq_frag.short_name, "23227_GGACATCA_L001_R1.fastq.gz")

        hiseq_frag = FastqFragment(url, merge_lanes=True)
        self.assertEqual(hiseq_frag.key, ("23227", "GGACATCA", "R1"))
        self.assertEqual(hiseq_frag.short_name, "23227_GGACATCA_R1.fastq.gz")

    def test_fastq_fragment_nextseq(self):
        url = "https://example.edu/flowcells/211011_VH00472_2_AAAJ5NWHV/23196_indexP27C3_2_S4_R2_001.fastq.gz"
        nextseq_frag = FastqFragment(url)
        nextseq_frag = FastqFragment(url)
        self.assertEqual(nextseq_frag.url, url)
        self.assertEqual(nextseq_frag.merge_lanes, False)
        self.assertEqual(nextseq_frag.name, "23196_indexP27C3_2_S4_R2_001.fastq.gz")
        self.assertEqual(nextseq_frag.library_id, "23196")
        self.assertEqual(nextseq_frag.index, "indexP27C3.2")
        self.assertEqual(nextseq_frag.lane, None)
        self.assertEqual(nextseq_frag.read, "R2")
        self.assertEqual(nextseq_frag.chunk, "001")
        self.assertEqual(nextseq_frag.key, ("23196", "indexP27C3.2", "R2"))
        self.assertEqual(nextseq_frag.short_name, "23196_indexP27C3.2_R2.fastq.gz")

        nextseq_frag = FastqFragment(url, merge_lanes=True)
        self.assertEqual(nextseq_frag.key, ("23196", "indexP27C3.2", "R2"))
        self.assertEqual(nextseq_frag.short_name, "23196_indexP27C3.2_R2.fastq.gz")

    def test_fastq_fragment_nextseq_undetermined(self):
        url = "https://example.edu/flowcells/211011_VH00472_2_AAAJ5NWHV/Undetermined_S0_R1_001.fastq.gz"
        self.assertRaises(ValueError, FastqFragment, url)

    def test_fastq_fragment_nextseq_too_short(self):
        url = "https://example.edu/flowcells/210416_SN787_1327_BHKH2HBCX3/22160_GAACGAAG_R1.fastq.gz"
        self.assertRaises(ValueError, FastqFragment, url)
