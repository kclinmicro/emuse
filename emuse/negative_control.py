import pandas as pd


def is_spike(row, spike_species):
    return row["species"] in spike_species


def absent_in_negative_control(row, neg_control):
    return row["species"] not in neg_control["species"].values


def is_low_abundance(row, cutoff):
    return row["abundance"] < cutoff


def is_enriched(row, sample, neg_control, normalising_spike_species, fold):
    # No spike configured -> do nothing
    if not normalising_spike_species:
        return False

    # Get abundance of spike in sample
    sample_spike = sample.loc[
        sample["species"] == normalising_spike_species, "abundance"
    ]

    # Get abundance of spike in neg control
    control_spike = neg_control.loc[
        neg_control["species"] == normalising_spike_species, "abundance"
    ]

    # Get abundance of species in neg control
    control_match = neg_control.loc[
        neg_control["species"] == row["species"], "abundance"
    ]

    # If any of these are empty, we can't do the calculation, so we return no highlight
    if sample_spike.empty or control_spike.empty or control_match.empty:
        return False

    # Normalise (species / spike) in sample and control, then compare
    sample_ratio = row["abundance"] / sample_spike.iloc[0]
    control_ratio = control_match.iloc[0] / control_spike.iloc[0]

    return control_ratio > 0 and sample_ratio > fold * control_ratio


def compare_to_negative_controls(sample_rows, controls, spike_species=(), normalising_spike_species=None,
                                 fold=25, low_abundance_cutoff=0.005):
    columns = ["species", "abundance"]
    sample = pd.DataFrame(sample_rows, columns=columns)
    control_tables = [pd.DataFrame(rows, columns=columns) for rows in controls]

    flags = []
    for _, row in sample.iterrows():
        present_in = [control for control in control_tables if not absent_in_negative_control(row, control)]
        flags.append({
            "spike": is_spike(row, spike_species),
            "absent_in_controls": bool(control_tables) and not present_in,
            "enriched": bool(present_in) and all(
                is_enriched(row, sample, control, normalising_spike_species, fold) for control in present_in
            ),
            "low_abundance": bool(is_low_abundance(row, low_abundance_cutoff)),
        })
    return flags
