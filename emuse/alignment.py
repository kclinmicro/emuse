import statistics

import pandas as pd
import pysam

from emuse.taxonomy import TaxTranslator


def get_alignment_metrics(sample_name, input_dir):
    abundance_path = f"{input_dir}/results/{sample_name}_downsampled.fastq_rel-abundance.tsv"
    assignment_path = f"{input_dir}/results/{sample_name}_downsampled.fastq_read-assignment-distributions.tsv"
    alignments_path = f"{input_dir}/results/{sample_name}_downsampled.fastq_emu_alignments.sam"

    df_abundance_unsorted = pd.read_csv(abundance_path, sep="\t")
    df_abundance = df_abundance_unsorted.sort_values("abundance", ascending=False)

    df_reads = load_read_file(assignment_path, df_abundance)
    colnames_taxids = df_reads.columns

    align_file = pysam.AlignmentFile(alignments_path)

    alns_all = {}
    for aln in align_file:
        readid = aln.query_name  # Ex: b9bb144e-eb53-4509-8931-5f4477444a48
        refname = aln.reference_name  # Ex: 562:emu_db:23853
        taxid = str(aln.reference_name).split(":")[0]  # Ex: 562
        refid = aln.reference_id  # Ex: 23853

        if aln.is_secondary or aln.is_supplementary:
            # We don't count these
            continue

        if taxid not in alns_all:
            alns_all[taxid] = []
        alns_all[taxid].append(aln)

    taxtr = TaxTranslator()
    aln_infos = []
    i = 1
    for taxid in colnames_taxids:
        if taxid in alns_all:
            alns = alns_all[taxid]
            identities, coverages = collect_distribution(alns, align_file)
            median_id = statistics.median(identities)
            median_cov = statistics.median(coverages)
            abundance = float(df_abundance[df_abundance["tax_id"] == taxid]["abundance"].values[0])
            taxon = taxtr.taxid_to_label(taxid)
            aln_infos.append(
                {
                    "tax id": taxid,
                    "median aligned identity": median_id,
                    "median aligned coverage": median_cov,
                }
            )

    df_aln_metrics = pd.DataFrame(aln_infos)
    return df_aln_metrics


def load_read_file(readassmt_path, df_abundance):
    df_reads = pd.read_csv(readassmt_path, sep="\t", header=0)
    colnames_sorted = [cn for cn in df_abundance["tax_id"] if cn in df_reads.columns]
    df_reads = df_reads[colnames_sorted]
    return df_reads


def collect_distribution(alns, alignment_file):
    identities = []
    coverages = []
    for aln in alns:
        identity, coverage = get_align_stats(aln, alignment_file)
        identities.append(identity)
        coverages.append(coverage)
    return identities, coverages


def get_align_stats(alignment, alignment_file):
    CIGAROP_MATCH = 0
    CIGAROP_INS = 1
    CIGAROP_DEL = 2
    CIGAROP_REF_SKIP = 3
    CIGAROP_SOFTCLIP = 4
    CIGAROP_HARDCLIP = 5
    CIGAROP_PAD = 6
    CIGAROP_EQUAL = 7
    CIGAROP_DIFF = 8
    CIGAROP_BACK = 9

    cigar_stats = alignment.get_cigar_stats()[0]
    edit_distance = alignment.get_tag("NM") if alignment.has_tag("NM") else None

    cigar_stats = alignment.get_cigar_stats()[0]

    insertions = cigar_stats[CIGAROP_INS]
    deletions = cigar_stats[CIGAROP_DEL]
    soft_clips = cigar_stats[CIGAROP_SOFTCLIP]

    query_len = alignment.query_length
    query_alignment_len = alignment.query_alignment_length

    reference_len = alignment_file.get_reference_length(alignment.reference_name)
    reference_alignment_len = alignment.reference_length

    identity, ref_coverage = calculate_align_stats(query_len,
                                                   query_alignment_len,
                                                   reference_len,
                                                   reference_alignment_len,
                                                   insertions,
                                                   deletions,
                                                   edit_distance,
                                                   soft_clips)

    return identity, ref_coverage



def calculate_align_stats(query_len,
                          query_alignment_len,
                          reference_len,
                          reference_alignment_len,
                          insertions,
                          deletions,
                          edit_distance,
                          soft_clips):
    """
    Return identity and coverage against the reference sequence

    Note that the identity in this calculation differs from the one in BLAST.
    While BLAST's identity measure does not count indels and instead divides
    the number of identical bases over the number of aligned bases, the
    identity metric below includes each aligned column corresponding to a base
    pair position in the reference or read, and counts the number of matches
    (identical bases) across the total aligned length, which is the length
    where gaps in both the reference and query are included.

    Args:
        query_len: The length in base pairs of the query sequence.
        query_alignment_len: The length in base pairs of the part of the query
            sequence covered by the alignment.
        reference_len: The length in base pairs of the reference sequence.
        reference_alignment_len: The length in base pairs of the part of the
            reference sequence covered by the alignment.
        insertions: Number of insertions measured in base pairs.
        deletions: Number of deletions measured in base pairs.
        edit_distance: Edit distance corresponding to the NM tag in SAM files.
        soft_clips: Bases in the ends of the alignment which do not align.

    Returns:
        identity: Identity as matching bases across the alignment length which
            contains both insertions and deletions.
        coverage: Percent of bases of the reference length covered by the
            alignment.
    """

    # ------------------------------------------------
    # Calculate identity
    # ------------------------------------------------
    matches = 0
    mismatches = 0
    if edit_distance is not None:
        mismatches = edit_distance - insertions - deletions
        matches = query_alignment_len - insertions - mismatches

    # We can not normalize over query or reference length, as these
    # won't contain either insertions or deletions
    divisor = matches + mismatches + insertions + deletions

    identity = 0
    if divisor > 0:
        identity = matches / divisor

    # ------------------------------------------------
    # Calculate coverage
    # ------------------------------------------------
    coverage = 0
    if reference_len > 0:
        coverage = reference_alignment_len / reference_len

    return identity, coverage
