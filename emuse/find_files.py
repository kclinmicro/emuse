from pathlib import Path

def find_emu_file(results_dir, sample_name, suffix):
    preprocessing_suffixes = [
        "_downsampled.fastq_",
        "_filtered.fastq_",
        ".porechop_abi.fastq_",
        ".trimmed.trim.fastq_",
        f".trimmed_1.trim.fastq-{sample_name}.trimmed_2.trim.fastq_",
        ".fastq_"
    ]

    files = [
        Path(results_dir)/f"{sample_name}{preprocessing}{suffix}"
        for preprocessing in preprocessing_suffixes
        if (Path(results_dir)/f"{sample_name}{preprocessing}{suffix}").exists()
    ]

    if len(files) != 1:
        raise ValueError(
            f"Expected exactly one {suffix} file for {sample_name}, found: {files}"
        )
    return files[0]
