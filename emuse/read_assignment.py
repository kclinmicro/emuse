import pandas as pd


def read_assignment_summary(path):
    assignment = pd.read_csv(path, sep="\t")
    # Select all columns except the first one
    assignment_filtered = assignment.iloc[:, 1:]
    # Compute mean and median for each column
    assignment_summary = assignment_filtered.agg(['median', 'mean']).T.reset_index()
    # Rename columns
    assignment_summary.columns = ['tax id', 'median probability*', 'mean probability*']
    return assignment_summary
