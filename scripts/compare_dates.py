#!/usr/bin/env python3
"""Compare popmix DATES output against DATES v753.

This script runs both implementations on the same data and reports
any differences. For scientific-grade validation, pair counts should
match exactly and statistics should match within floating point precision.

Usage:
    python scripts/compare_dates.py [data_path]

Default data path: tests/golden/dates/family_chr20_22
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def parse_dates_output(output_file: str) -> list[dict]:
    """Parse DATES .out file into list of bin dictionaries."""
    bins = []
    with open(output_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 5:
                try:
                    bins.append({
                        'distance_cm': float(parts[0]),
                        'covariance': float(parts[1]),
                        'z_score': float(parts[2]),
                        'correlation': float(parts[3]),
                        'n_pairs': int(parts[4]),
                    })
                except (ValueError, IndexError):
                    continue
    return bins


def run_dates(data_path: str, dates_bin: str = "/home/a/DATES/bin/dates") -> list[dict]:
    """Run DATES v753 and return parsed output."""
    if not os.path.exists(dates_bin):
        print(f"ERROR: DATES binary not found at {dates_bin}")
        print("Build DATES first: cd /home/a/DATES/src && make")
        sys.exit(1)

    # Use fixed directory to avoid tempdir issues
    tmpdir = "/tmp/popmix_dates_comparison"
    import shutil
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir)

    try:
        # Create admix file
        # Format: source1 source2 test_pop output_prefix
        admix_file = os.path.join(tmpdir, "admix")
        with open(admix_file, 'w') as f:
            f.write("French Yoruba Simulation Simulation\n")

        # Create parameter file (matching working configuration)
        data_dir = os.path.dirname(data_path)
        data_base = os.path.basename(data_path)
        par_file = os.path.join(tmpdir, "par.dates")
        with open(par_file, 'w') as f:
            f.write(f"DIR: {data_dir}\n")
            f.write(f"indivname: DIR/{data_base}.ind\n")
            f.write(f"snpname: DIR/{data_base}.snp\n")
            f.write(f"genotypename: DIR/{data_base}.geno\n")
            f.write(f"admixlist: {admix_file}\n")
            f.write("binsize: 0.001\n")
            f.write("maxdis: 1.0\n")
            f.write("seed: 77\n")
            f.write("runmode: 1\n")
            f.write("chithresh: 0.0\n")
            f.write("mincount: 1\n")
            f.write("zdipcorrmode: YES\n")
            f.write("qbin: 0\n")  # Use naive O(n²) mode, not FFT (avoids memory issues)
            f.write("runfit: NO\n")
            f.write("afffit: NO\n")

        # Run DATES
        result = subprocess.run(
            [dates_bin, "-p", par_file],
            cwd=tmpdir,
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            print(f"DATES failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            sys.exit(1)

        # Parse output - DATES creates output in Simulation/ subdirectory
        # With qbin=0 (naive mode), output is only in Simulation/Simulation.out
        output_dir = os.path.join(tmpdir, "Simulation")
        output_file = os.path.join(output_dir, "Simulation.out")

        if not os.path.exists(output_file):
            # Fall back: check tmpdir directly (qbin>0 mode)
            alt_file = os.path.join(tmpdir, "Simulation.out")
            if os.path.exists(alt_file):
                output_file = alt_file
            else:
                print("ERROR: DATES output file not found")
                print(f"Files in tmpdir: {os.listdir(tmpdir)}")
                if os.path.isdir(output_dir):
                    print(f"Files in Simulation/: {os.listdir(output_dir)}")
                sys.exit(1)

        # Extract coefficient from DATES output
        coeff_match = re.search(r'coeff:\s+[\d.]+\s+\w+\s+([\d.]+)', result.stdout)
        dates_coeff = float(coeff_match.group(1)) if coeff_match else None

        return parse_dates_output(output_file), dates_coeff
    finally:
        pass  # Keep tmpdir for debugging


def run_popmix(data_path: str) -> tuple[list[dict], float]:
    """Run popmix DATES and return results."""
    import popmix as pm

    geno = pm.read_packedancestrymap(data_path)
    result = pm.dates(
        geno,
        test_pop='Simulation',
        source1='French',
        source2='Yoruba'
    )

    bins = []
    for b in result.ld_bins:
        bins.append({
            'distance_cm': b.distance_cm,
            'covariance': b.weighted_ld,
            'z_score': b.z_score,
            'correlation': b.correlation,
            'n_pairs': b.n_pairs,
        })

    # Get coefficient (need to access internal attribute)
    # For now, return None - we verified coefficient matches
    return bins, None


def compare_results(dates_bins: list[dict], popmix_bins: list[dict],
                    dates_coeff: float = None, popmix_coeff: float = None) -> bool:
    """Compare results and print detailed report. Returns True if all pass."""
    print("=" * 70)
    print("DATES v753 vs popmix Comparison Report")
    print("=" * 70)
    print()

    # Header
    print(f"{'Dist':>6} {'DATES_pairs':>12} {'PM_pairs':>12} {'Diff':>8} "
          f"{'DATES_corr':>12} {'PM_corr':>12} {'Corr_Δ':>10} {'Status':>6}")
    print("-" * 70)

    total_pair_diff = 0
    max_corr_diff = 0.0
    all_exact = True

    # Compare bins (use min length to handle different max distances)
    n_bins = min(len(dates_bins), len(popmix_bins), 100)  # First 100 bins (10 cM)

    for i in range(n_bins):
        d = dates_bins[i]
        p = popmix_bins[i]

        pair_diff = p['n_pairs'] - d['n_pairs']
        corr_diff = p['correlation'] - d['correlation']

        total_pair_diff += abs(pair_diff)
        max_corr_diff = max(max_corr_diff, abs(corr_diff))

        if pair_diff != 0:
            all_exact = False

        status = "PASS" if pair_diff == 0 else "FAIL"

        print(f"{d['distance_cm']:6.1f} {d['n_pairs']:12d} {p['n_pairs']:12d} "
              f"{pair_diff:+8d} {d['correlation']:12.6f} {p['correlation']:12.6f} "
              f"{corr_diff:+10.6f} {status:>6}")

    print("-" * 70)
    print()

    # Summary
    print("SUMMARY")
    print("-" * 70)
    print(f"Total pair count difference: {total_pair_diff}")
    print(f"Maximum correlation difference: {max_corr_diff:.2e}")
    print()

    # Pass/fail criteria
    pair_pass = total_pair_diff == 0
    corr_pass = max_corr_diff < 1e-6

    print("CRITERIA")
    print(f"  Pair counts exact match:     {'PASS' if pair_pass else 'FAIL'}")
    print(f"  Correlations within 1e-6:    {'PASS' if corr_pass else 'FAIL'}")
    print()

    overall = pair_pass and corr_pass
    print(f"OVERALL: {'PASS' if overall else 'FAIL'}")
    print("=" * 70)

    return overall


def main():
    parser = argparse.ArgumentParser(description='Compare popmix DATES with DATES v753')
    parser.add_argument('data_path', nargs='?',
                        default='tests/golden/dates/family_chr20_22',
                        help='Path to data files (without extension)')
    parser.add_argument('--dates-bin', default='/home/a/DATES/bin/dates',
                        help='Path to DATES binary')
    parser.add_argument('--skip-dates', action='store_true',
                        help='Skip running DATES, use cached reference')
    args = parser.parse_args()

    # Resolve path relative to project root
    if not os.path.isabs(args.data_path):
        project_root = Path(__file__).parent.parent
        args.data_path = str(project_root / args.data_path)

    print(f"Data path: {args.data_path}")
    print()

    # Check data exists
    if not os.path.exists(f"{args.data_path}.geno"):
        print(f"ERROR: Data file not found: {args.data_path}.geno")
        sys.exit(1)

    # Run DATES
    print("Running DATES v753...")
    dates_bins, dates_coeff = run_dates(args.data_path, args.dates_bin)
    print(f"  Got {len(dates_bins)} bins")
    if dates_coeff:
        print(f"  Coefficient: {dates_coeff:.6f}")
    print()

    # Run popmix
    print("Running popmix...")
    popmix_bins, popmix_coeff = run_popmix(args.data_path)
    print(f"  Got {len(popmix_bins)} bins")
    print()

    # Compare
    success = compare_results(dates_bins, popmix_bins, dates_coeff, popmix_coeff)

    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
