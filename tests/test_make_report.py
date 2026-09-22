# TBC ...
import subprocess as sp
from trana_vy.make_report import calculate_align_stats, get_align_stats
import pysam
from pytest import fail
import tempfile

def test_get_align_stats():
    alignment_file_path = create_alignment_file(ref_len=100, query_len=60)
    aln_file = pysam.AlignmentFile(alignment_file_path)

    want_perc_id = 1.0
    want_perc_cov = 0.6

    for aln in aln_file:
        perc_id, perc_cov = get_align_stats(aln, aln_file)

        if perc_id != want_perc_id:
            fail(f"Expected percent identity {want_perc_id} but got {perc_id}")
        if perc_cov != want_perc_cov:
            fail(f"Expected percent identity {want_perc_cov} but got {perc_cov}")

        pass

def test_calculate_align_stats():
    for t in [
        {
            "qlen": 10,
            "qalnlen": 10,
            "reflen": 20,
            "refalnlen": 10,
            "nins": 0,
            "ndel": 0,
            "editdist": 0,
            "softclips": 0,
            "want_identity": 1,
            "want_coverage": 0.5,
        },
    ]:
        id, cov = calculate_align_stats(t["qlen"], t["qalnlen"], t["reflen"],
                                        t["refalnlen"], t["nins"], t["ndel"],
                                        t["editdist"], t["softclips"])

        if id != t["want_identity"]:
            fail(f"Wanted identity {t['want_identity']} but got {id}")

        if cov != t["want_coverage"]:
            fail(f"Wanted coverage {t['want_coverage']} but got {cov}")


def create_alignment_file(ref_len, query_len):
    with tempfile.TemporaryDirectory(prefix="tranavy_tests_", delete=False) as tmpdir:
        seq = "CGGCTTAGAGGCGGCTGCGCGTAGTGCTGCTGATTATATTTCGGCGCTATATCGCTGATGATCGTCGATACGTACTGCTAGCTGCATCTGATCGATCGTC"
        # ------------------------------------------------
        # Create a Fasta file
        # ------------------------------------------------
        db_lines = [
            ">dbseq1\n",
            f"{seq[:ref_len]}\n",
        ]
        db_file = f"{tmpdir}/db.fa"
        with open(db_file, "w") as dbfile:
            dbfile.writelines(db_lines)

        # ------------------------------------------------
        # Create a FastQ file
        # ------------------------------------------------
        fastq_lines = [
            "@qseq1\n",
            f"{seq[:query_len]}\n",
            "+\n",
            "5"*query_len+"\n",
        ]
        reads_file = f"{tmpdir}/reads.fq"
        with open(reads_file, "w") as rfile:
            rfile.writelines(fastq_lines)

        sam_align_file = f"{tmpdir}/alignments.sam"

        sp.check_output((
            "minimap2 "
            "-a "         # Output in .sam format
            "-x map-ont " # map-ont
            "-t 1 "       # Threads
            "-N 50 "
            "-p .9 "
            "-u f "
            "-K 500000000 "
            f"{db_file} "
            f"{reads_file} "
            f"-o {sam_align_file}"
        ), shell=True)

        return sam_align_file
