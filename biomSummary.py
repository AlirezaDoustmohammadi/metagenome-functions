import biom
import pandas as pd
import argparse
import os


def asv_summary():

    # (transpose so that rows correspond to samples and columns to ASVs)
    df_features = df_samples.transpose()
    # Compute summary statistics for each sample
    features_stat = df_features.agg(['sum', 'max'])
    # Count the number of samples that contain each ASV
    asv_sample_counts = (df_features > 0).sum()
    features_stat = features_stat.transpose()
    features_stat['Number of Samples'] = asv_sample_counts.tolist()
    # rename columns
    features_stat = features_stat.rename(columns={'max': 'Maximum abundance in samples',
                                                  'sum': 'Total abundance of ASVs'})
    features_stat = features_stat.sort_values(by=['Total abundance of ASV'])
    # write to csv
    features_stat.to_csv(os.path.basename(args.biom_file).split('.biom')[0] + '.summary.features.csv')
    df_features.to_csv(os.path.basename(args.biom_file).split('.biom')[0] + '.features.csv')


def samples_summary():

    # Compute summary statistics for each sample
    samples_stats = df_samples.agg(['sum'])
    # Count the number of samples that contain each ASV
    feature_counts = (df_samples > 0).sum()
    samples_stats = samples_stats.transpose()
    samples_stats['Number of ASVs'] = feature_counts.tolist()
    # rename column
    samples_stats = samples_stats.rename(columns={'sum': 'Total abundance of ASVs'})
    # sort
    samples_stats = samples_stats.sort_values(by=['Total abundance of ASVs'])
    # write to csv
    samples_stats.to_csv(os.path.basename(args.biom_file).split('.biom')[0] + '.summary.samples.csv')
    df_samples.to_csv(os.path.basename(args.biom_file).split('.biom')[0] + '.samples.csv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summary feature table file')
    parser.add_argument('--biom_file', type=str, required=True, help='Path to the feature table file')

    args = parser.parse_args()

    # Load the table
    table = biom.load_table(args.biom_file)

    # Convert to a pandas DataFrame
    df_samples = table.to_dataframe()

    # call functions
    samples_summary()
    asv_summary()
