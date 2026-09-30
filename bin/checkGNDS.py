#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from fudge import warning as warningModule
from fudge import GNDS_file as GNDS_fileModule

summaryDocString__FUDGE = '''Run FUDGE physics tests on one or more GNDS files.'''

description1 = '''Read one or more GNDS files into Fudge and run all physics tests.
    Sample use: python checkGNDS.py n-001_H_001.xml n-001_H_002.xml ...
    If file n-001_H_001-cov.xml (or -covar.xml) exists, covariances will automatically be read and checked.

    If only one GNDS file is supplied, results will be printed to the console.
    If multiple files were supplied, results will be saved in the outputDir,
    along with a .csv file summarizing the number of warnings by severity level in each file.
'''

__doc__ = description1

parser = argparse.ArgumentParser(description1)
parser.add_argument('gnds', nargs='+',                                  help='GNDS and/or PoPs file(s) to check.')
parser.add_argument('-e', '--ebalance', action='store_true',            help='Include energy balance warnings (takes longer).')
parser.add_argument('--normTolerance', type=float, default=1e-5,        help='Tolerance for warning about unnormalized distributions. Default=1e-5')
parser.add_argument('--threshold', default='Moderate',                  help='Minimum warning threshold, e.g. Pedantic, Minor, Moderate, Severe, Fatal')
parser.add_argument('--outputDir', default='checking_results',          help='Directory for output files (CSV and log files). Default=checking_results')
parser.add_argument('--parallel', type=int, default=1,                  help='Number of parallel processes to use for checking. Default=1 (sequential)')
parser.add_argument('-f', '--failOnException', action='store_true',     help='Fail immediately rather than converting Exception to a warning')
parser.add_argument('-v', '--verbose', action='store_true',             help='Print extra information while checks are running.')


def count_warnings_by_severity(warnings):
    """Count warnings by severity level."""
    counts = {
        warningModule.Level.Pedantic: 0,
        warningModule.Level.Minor: 0,
        warningModule.Level.Moderate: 0,
        warningModule.Level.Severe: 0,
        warningModule.Level.Fatal: 0,
    }

    def count_recursive(warning_list):
        for item in warning_list:
            if isinstance(item, warningModule.Context):
                count_recursive(item.warningList)
            elif hasattr(item, 'level'):
                counts[item.level] = counts.get(item.level, 0) + 1

    count_recursive(warnings.warningList)
    return counts

def check_single_file(fileName, threshold, ebalance, failOnException, normTolerance, verbose, outputDir):
    """Check a single GNDS file and return results."""
    try:
        gnds = GNDS_fileModule.read(fileName)

        warnings = gnds.check(checkEnergyBalance=ebalance, failOnException=failOnException,
                              normTolerance=normTolerance, verbose=verbose)
        filtered, screened = warnings.filter(threshold=threshold)

        # Write warnings to log file
        base_name = Path(fileName).stem
        log_file = Path(outputDir) / f"{base_name}.log"
        with open(log_file, 'w') as f:
            f.write(str(filtered))
            if screened:
                f.write("\n\n  Some warnings were screened\n")
                for key in screened:
                    f.write(f"    {key}: {screened[key]} occurrences\n")

        # Count warnings by severity
        rs_counts = count_warnings_by_severity(warnings)

        covariances = []
        if hasattr(gnds, 'loadCovariances'):
            covariances = gnds.loadCovariances()

        cov_counts = None
        for covarianceSuite in covariances:
            covWarnings = covarianceSuite.check(verbose=verbose)
            filtered_cov, screened_cov = covWarnings.filter(threshold=threshold)

            # Append covariance warnings to log file
            with open(log_file, 'a') as f:
                f.write(f'\n\nChecking covariance file {covarianceSuite.sourcePath}\n')
                f.write(str(filtered_cov))
                if screened_cov:
                    f.write("\n\n  Some covariance warnings were screened\n")
                    for key in screened_cov:
                        f.write(f"    {key}: {screened_cov[key]} occurrences\n")

            # Count covariance warnings by severity (only first covariance suite)
            if cov_counts is None:
                cov_counts = count_warnings_by_severity(covWarnings)

        # Build result row
        row = [fileName]
        row.extend([
            rs_counts[warningModule.Level.Pedantic],
            rs_counts[warningModule.Level.Minor],
            rs_counts[warningModule.Level.Moderate],
            rs_counts[warningModule.Level.Severe],
            rs_counts[warningModule.Level.Fatal]
        ])
        if cov_counts is not None:
            row.extend([
                cov_counts[warningModule.Level.Pedantic],
                cov_counts[warningModule.Level.Minor],
                cov_counts[warningModule.Level.Moderate],
                cov_counts[warningModule.Level.Severe],
                cov_counts[warningModule.Level.Fatal]
            ])
        else:
            row.extend([''] * 5)

        print(f"Completed checking {fileName}")
        return row

    except Exception as e:
        print(f"Error checking {fileName}: {e}")
        return None

if __name__ == '__main__':
    args = parser.parse_args()
    threshold = warningModule.Level.fromString(args.threshold.title())

    # Handle single file case - print to console only
    if len(args.gnds) == 1:
        fileName = args.gnds[0]
        gnds = GNDS_fileModule.read(fileName)

        warnings = gnds.check(checkEnergyBalance=args.ebalance, failOnException=args.failOnException,
                              normTolerance=args.normTolerance, verbose=args.verbose)
        filtered, screened = warnings.filter(threshold=threshold)
        print(filtered)
        if screened:
            print("\n  Some warnings were screened")
            for key in screened:
                print(f"    {key}: {screened[key]} occurrences")

        covariances = []
        if hasattr(gnds, 'loadCovariances'):
            covariances = gnds.loadCovariances()

        for covarianceSuite in covariances:
            print('\nChecking covariance file %s' % covarianceSuite.sourcePath)
            covWarnings = covarianceSuite.check(verbose=args.verbose)
            filtered, screened = covWarnings.filter(threshold=threshold)
            print(filtered)
            if screened:
                print("\n  Some covariance warnings were screened")
                for key in screened:
                    print(f"    {key}: {screened[key]} occurrences")
    else:
        # Multiple files - use outputDir and CSV
        # Create output directory if it doesn't exist
        output_dir = Path(args.outputDir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Setup CSV file
        csv_path = output_dir / 'summary.csv'
        file_exists = csv_path.is_file()
        csv_file = open(csv_path, 'a', newline='')
        csv_writer = csv.writer(csv_file)
        if not file_exists:
            # Write header
            csv_writer.writerow([
                'File',
                'RS_Pedantic', 'RS_Minor', 'RS_Moderate', 'RS_Severe', 'RS_Fatal',
                'Cov_Pedantic', 'Cov_Minor', 'Cov_Moderate', 'Cov_Severe', 'Cov_Fatal'
            ])

        # Process files
        if args.parallel > 1:
            # Parallel processing
            print(f"Processing {len(args.gnds)} files using {args.parallel} parallel processes...")
            with ProcessPoolExecutor(max_workers=args.parallel) as executor:
                # Submit all tasks and maintain order
                futures = [
                    executor.submit(
                        check_single_file,
                        fileName,
                        threshold,
                        args.ebalance,
                        args.failOnException,
                        args.normTolerance,
                        args.verbose,
                        str(output_dir)
                    ) for fileName in args.gnds
                ]

                # Collect results in original order
                for future in futures:
                    result = future.result()
                    if result is not None:
                        csv_writer.writerow(result)
        else:
            # Sequential processing
            for fileName in args.gnds:
                result = check_single_file(
                    fileName,
                    threshold,
                    args.ebalance,
                    args.failOnException,
                    args.normTolerance,
                    args.verbose,
                    str(output_dir)
                )
                if result is not None:
                    csv_writer.writerow(result)

        csv_file.close()
        print(f"\nSummary written to {csv_path}")
        print(f"Log files written to {output_dir}/")
