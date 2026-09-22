# TBC ...
import subprocess as sp

from make_report import calculate_align_stats, get_align_stats
import pysam
from pytest import fail
import tempfile

# An example sequence to be used in multiple tests
seq = "CGGCTTAGAGGCGGCTGCGCGTAGTGCTGCTGATTATATTTCGGCGCTATATCGCTGATGATCGTCGATACGTACTGCTAGCTGCATCTGATCGATCGTC"


def test_get_align_stats():
    for i, tc in enumerate([
        {
            "qry": seq[0:60],
            "ref": seq,
            "want_id": 1.0,
            "want_cov": 0.6,
        },
        {
            # Diffs:|
            #       v
            "qry": "-GGCTTAGAGGCGGCTGCGCGTAGTGCTGCTGATTATATTTCGGCGCTATATCGCTGATGATCGTCGATACGTACTGCTAGCTGCATCTGATCGATCGTC".replace("-", ""),
            "ref": "CGGCTTAGAGGCGGCTGCGCGTAGTGCTGCTGATTATATTTCGGCGCTATATCGCTGATGATCGTCGATACGTACTGCTAGCTGCATCTGATCGATCGTC".replace("-", ""),
            "want_id": 1.0,
            "want_cov": 0.99,
        },
        {
            # Diffs:    |
            #           v
            "qry": "CGGC-TAGAGGCGGCTGCGCGTAGTGCTGCTGATTATATTTCGGCGCTATATCGCTGATGATCGTCGATACGTACTGCTAGCTGCATCTGATCGATCGTC".replace("-", ""),
            "ref": "CGGCTTAGAGGCGGCTGCGCGTAGTGCTGCTGATTATATTTCGGCGCTATATCGCTGATGATCGTCGATACGTACTGCTAGCTGCATCTGATCGATCGTC".replace("-", ""),
            "want_id": 0.99,
            "want_cov": 1.0,
        },
        {
            # Diffs:       |                            |                               |                      |
            #              v                            v                               v                      v
            "qry": "CGGCTTA-AGGCGGCTGCGCGTAGTGCTGCTGATTATATTTCGGCGCTATATCGCTGATGATCGTCGA-ACGTACTGCTAGCTGCATCTGATCGATCGTC".replace("-", ""),
            "ref": "CGGCTTAGAGGCGGCTGCGCGTAGTGCTGCTGATTA-ATTTCGGCGCTATATCGCTGATGATCGTCGATACGTACTGCTAGCTGCATCTGA-CGATCGTC".replace("-", ""),
            "want_id": 0.96,
            # The -2 are becase we have to subtract any bases missing in the
            # reference, for the coverage calculation
            "want_cov": (100-2)/(100-2),
        },
        {
            # Diffs:                                           ||||||||||
            #                                                  vvvvvvvvvv
            "qry": "CGGCTTAGAGGCGGCTGCGCGTAGTGCTGCTGATTATATTTCGGCGCTATATCGCTGATGATCGTCGATACGTACTGCTAGCTGCATCTGATCGATCGTC".replace("-", ""),
            "ref": "CGGCTTAGAGGCGGCTGCGCGTAGTGCTGCTGATTATATTTCG----------GCTGATGATCGTCGATACGTACTGCTAGCTGCATCTGATCGATCGTC".replace("-", ""),
            "want_id": 0.9,
            # The -10 are becase we have to subtract any bases missing in the
            # reference, for the coverage calculation
            "want_cov": (100-10)/(100-10),
        },
        {
            # Diffs:                                           ||||||||||
            #                                                  vvvvvvvvvv
            "qry": "CGGCTTAGAGGCGGCTGCGCGTAGTGCTGCTGATTATATTTCG----------GCTGATGATCGTCGATACGTACTGCTAGCTGCATCTGATCGATCGTC".replace("-", ""),
            "ref": "CGGCTTAGAGGCGGCTGCGCGTAGTGCTGCTGATTATATTTCGGCGCTATATCGCTGATGATCGTCGATACGTACTGCTAGCTGCATCTGATCGATCGTC".replace("-", ""),
            "want_id": 0.9,
            "want_cov": 1.0,
        },
        {
            # Diffs:            |          |                                          |
            #                   v          v                                          v
            "qry": "CGGCTTAGAGGC-GCTGCGCGTAGTGCTGCTGATTATATTTCGGCGCTATATCGCTGATGATCGTC-ATACGTACTGCTAGCTAGCTGCA----------".replace("-", ""),
            "ref": "CGGCTTAGAGGCGGCTGCGCGTA-TGCTGCTGATTATATTTCGGCGCTATATCGCTGATGATCGTCGATACGTACTGCTAGCTAGCTGCATCTGATCGTC".replace("-", ""),
            "want_id": (90-3)/90,
            # The -1 are becase we have to subtract any bases missing in the
            # reference, for the coverage calculation
            "want_cov": (90-1)/(100-1),
        },
        {
            # Diffs:                                      |
            #                                             v
            "qry": "--------------------GTAGTGCTGCTGATTATACTTCGGCGCTATATCGCTGATGATCGTCGATACGTACTGCTAGCTAGCTGCA----------".replace("-", ""),
            "ref": "CGGCTTAGAGGCGGCTGCGCGTAGTGCTGCTGATTATATTTCGGCGCTATATCGCTGATGATCGTCGATACGTACTGCTAGCTAGCTGCATCTGATCGTC".replace("-", ""),
            "want_id": (70-1)/70,
            "want_cov": 0.7,
        },
        {
            # Diffs:                       |                      |                   |
            #                              v                      v                   v
            "qry": "-----TAGAGGCGGCTGCGCGTAGTGCTGCTGATTATATTTCGGCGGTATATCGCTGATGATCGTC-ATACGTACTGCTAGCTAGCTGCA----------".replace("-", ""),
            "ref": "CGGCTTAGAGGCGGCTGCGCGTA-TGCTGCTGATTATATTTCGGCGCTATATCGCTGATGATCGTCGATACGTACTGCTAGCTAGCTGCATCTGATCGTC".replace("-", ""),
            "want_id": (85-3)/85,
            # The -1 are becase we have to subtract any bases missing in the
            # reference, for the coverage calculation
            "want_cov": (85-1)/(100-1),
        },
    ], start=1):
        aln_path = create_alignment_file(tc["qry"], tc["ref"])
        aln_file = pysam.AlignmentFile(aln_path)

        diff_delta = 1e-9

        for j, aln in enumerate(aln_file, start=1):
            perc_id, perc_cov = get_align_stats(aln, aln_file)
            if abs(perc_id - tc["want_id"]) >= diff_delta:
                fail(f"For test case {i} alignment {j}: Expected percent identity {tc['want_id']} but got {perc_id}")
            if abs(perc_cov - tc["want_cov"]) >= diff_delta:
                fail(f"For test case {i} alignment {j}: Expected percent coverage {tc['want_cov']} but got {perc_cov}")


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
        id, cov = calculate_align_stats(
            t["qlen"],
            t["qalnlen"],
            t["reflen"],
            t["refalnlen"],
            t["nins"],
            t["ndel"],
            t["editdist"],
            t["softclips"],
        )

        if id != t["want_identity"]:
            fail(f"Wanted identity {t['want_identity']} but got {id}")

        if cov != t["want_coverage"]:
            fail(f"Wanted coverage {t['want_coverage']} but got {cov}")


def create_alignment_file_subset(ref_len, query_len):
    ref_seq = seq[0:ref_len]
    query_seq = seq[0:query_len]

    return create_alignment_file(query_seq, ref_seq)


def create_alignment_file(read_seq, db_seq):
    with tempfile.TemporaryDirectory(prefix="tranavy_tests_", delete=False) as tmpdir:
        # ------------------------------------------------
        # Create a Fasta file
        # ------------------------------------------------
        db_lines = [
            ">dbseq1\n",
            f"{db_seq}\n",
        ]
        db_file = f"{tmpdir}/db.fa"
        with open(db_file, "w") as dbfile:
            dbfile.writelines(db_lines)

        # ------------------------------------------------
        # Create a FastQ file
        # ------------------------------------------------
        fastq_lines = [
            "@qseq1\n",
            f"{read_seq}\n",
            "+\n",
            "5" * len(read_seq) + "\n",
        ]
        reads_file = f"{tmpdir}/reads.fq"
        with open(reads_file, "w") as rfile:
            rfile.writelines(fastq_lines)

        sam_align_path = f"{tmpdir}/alignments.sam"

        sp.check_output(
            (
                "minimap2 "
                "-a "  # Output in .sam format
                "-x map-ont "  # map-ont
                "-t 1 "  # Threads
                "-N 50 "
                "-p .9 "
                "-u f "
                "-K 500000000 "
                f"{db_file} "
                f"{reads_file} "
                f"-o {sam_align_path}"
            ),
            shell=True,
        )

        return sam_align_path
