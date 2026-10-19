def is_spike(row, spike_species):
    return row["species"] in spike_species


def absent_in_negative_control(row, neg_control):
    return row["species"] not in neg_control["species"].values


def is_low_abundance(row, cutoff):
    return row["abundance"] < cutoff


def is_enriched(row, sample, neg_control, normalising_spike_species, fold=25):
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
