# TBC ...
import subprocess as sp
from importlib.metadata import version

from emuse.make_report import main
from emuse.alignment import calculate_align_stats, get_align_stats
from emuse.find_files import find_emu_file
from emuse.negative_control import compare_to_negative_controls
import pandas as pd
import pysam
from pytest import fail, mark, raises
import tempfile
from pathlib import Path
from jinja2 import Environment, FileSystemLoader

TEMPLATE_DIR = Path(__file__).parents[1] / "emuse" / "data" / "templates"
template = Environment(loader=FileSystemLoader(TEMPLATE_DIR)).get_template("report.html.j2")

# An example sequence to be used in multiple tests
seq = "CGGCTTAGAGGCGGCTGCGCGTAGTGCTGCTGATTATATTTCGGCGCTATATCGCTGATGATCGTCGATACGTACTGCTAGCTGCATCTGATCGATCGTC"


@mark.parametrize("flag", ["-v", "--version"])
def test_cli_version(flag, monkeypatch, capsys):
    monkeypatch.setattr("sys.argv", ["emuse", flag])

    with raises(SystemExit) as exit_info:
        main()

    assert exit_info.value.code == 0
    output = capsys.readouterr()
    assert output.out == f"emuse {version('emuse')}\n"
    assert output.err == ""


def test_find_emu_file(tmp_path):
    sample_name = "sample1"
    suffix = "rel-abundance.tsv"
    expected = tmp_path / f"{sample_name}_downsampled.fastq_{suffix}"
    expected.touch()

    assert find_emu_file(tmp_path, sample_name, suffix) == expected


def test_find_emu_file_multiple_matches(tmp_path):
    sample_name = "sample1"
    suffix = "rel-abundance.tsv"
    (tmp_path / f"{sample_name}_downsampled.fastq_{suffix}").touch()
    (tmp_path / f"{sample_name}_filtered.fastq_{suffix}").touch()

    with raises(ValueError, match="Expected exactly one"):
        find_emu_file(tmp_path, sample_name, suffix)


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


def test_compare_to_negative_controls():
    species = ["Staphylococcus aureus", "Escherichia coli", "Imtechella halotolerans"]
    sample = pd.DataFrame({
        "species": species + ["Cutibacterium acnes"],
        "abundance": [0.5, 0.3, 0.1, 0.1],
    })
    controls = [
        pd.DataFrame({"species": species, "abundance": [0.001, 0.001, 0.5]}),
        pd.DataFrame({"species": species, "abundance": [0.002, 0.4, 0.5]}),
    ]

    flags = compare_to_negative_controls(sample, controls, normalising_spike_species="Imtechella halotolerans")

    assert flags["enriched"].tolist() == [True, False, False, False]
    assert flags["absent_in_controls"].tolist() == [False, False, False, True]

# This test tests the behavior of compare_to_negative_controls when no negative controls are provided. It checks that the flags for "absent_in_controls" and "enriched" are correctly set to False for all species in the sample.

def test_compare_to_negative_controls_without_controls():
    sample = pd.DataFrame({
        "species": ["Species A", "Species B"],
        "abundance": [0.5, 0.3],
    })

    flags = compare_to_negative_controls(sample, [])

    assert not flags["absent_in_controls"].any()
    assert not flags["enriched"].any()

# Testing template rendering

def render_report(sample_name, neg_controls, neg_control_tables):
    stats = {
        "reads_before_filtering": 1000,
        "reads_after_filtering": 900,
        "mean_read_length": 1000,
        "mean_read_quality": 20,
        "median_read_length": 1000,
        "median_read_quality": 20,
        "read_length_n50": 1000,
        "stdev_read_length": 100,
        "total_bases": 900000,
    }

    return template.render(
        css="",
        table="<table>sample_table</table>",
        neg_control_tables=neg_control_tables,
        legend="sample_legend",
        legend_neg="neg_legend",
        today="2026-11-23",
        pipeline_version="1.0",
        stats=stats,
        input_dir="test_run",
        sample_name=sample_name,
        neg_controls=neg_controls,
        analysis_type="16s",
        low_abundance_cutoff=0.005,
    )


def test_report_one_negative_control():
    html = render_report(
        "sample1",
        ["NEG1"],
        {"NEG1": "<table>NEG1_table</table>"},
    )

    assert "sample_table" in html
    assert "NEG1_table" in html
    assert "Negative control: NEG1" in html


def test_report_two_negative_controls():
    html = render_report(
        "sample1",
        ["NEG1", "NEG2"],
        {
            "NEG1": "<table>NEG1_table</table>",
            "NEG2": "<table>NEG2_table</table>",
        },
    )

    assert "sample_table" in html
    assert "NEG1_table" in html
    assert "NEG2_table" in html
    assert "Negative control: NEG1" in html
    assert "Negative control: NEG2" in html


def test_report_negative_control_sample():
    html = render_report(
        "NEG1",
        ["NEG1", "NEG2"],
        {
            "NEG1": "<table>NEG1_table</table>",
            "NEG2": "<table>NEG2_table</table>",
        },
    )

    assert "NEG1_table" in html
    assert "sample_table" not in html
    assert "NEG2_table" not in html


def test_report_without_negative_controls():
    html = render_report("sample1", [], {})

    assert "sample_table" in html
    assert "Negative control:" not in html
    assert "neg_legend" not in html
