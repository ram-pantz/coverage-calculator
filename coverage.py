import pandas as pd
from collections import defaultdict
import sys
import os

def calculate_population_coverage(allele_list: list, coverage_data_file: str) -> float:
    """
    Calculates the total population coverage for a given list of HLA alleles
    using data from a CSV file.

    Args:
        allele_list: A list of HLA alleles (e.g., ['A*03:02', 'A*68:01']).
        coverage_data_file: The path to the CSV file containing HLA coverage data.

    Returns:
        The total combined population coverage as a float.
    """
    if not os.path.exists(coverage_data_file):
        print(f"Error: The file '{coverage_data_file}' was not found.", file=sys.stderr)
        return 0.0
    
    try:
        df = pd.read_csv(coverage_data_file)
    except Exception as e:
        print(f"Error reading CSV file: {e}", file=sys.stderr)
        return 0.0

    # Normalize allele names to match the input format
    df['Allele'] = df['Allele'].str.replace('HLA-', '')
    
    # Create a dictionary for quick lookup of allele frequencies and loci
    allele_data = {
        row['Allele'].strip(): {'Locus': row['Locus'].strip(), 'Frequency': row['Frequency']}
        for index, row in df.iterrows()
    }
    
    # Group alleles by locus and sum their frequencies
    locus_frequencies = defaultdict(float)
    for allele in allele_list:
        clean_allele = allele.strip()
        if clean_allele in allele_data:
            data = allele_data[clean_allele]
            locus = data['Locus']
            frequency = data['Frequency']
            locus_frequencies[locus] += frequency
        else:
            print(f"Warning: Allele '{clean_allele}' not found in the data file. Skipping.", file=sys.stderr)
            
    # Calculate population coverage for each locus using Hardy-Weinberg principle
    locus_coverages = {}
    for locus, p in locus_frequencies.items():
        # The formula used is 1 - (1-p)^2, which is equivalent to 2p - p^2
        coverage = 1 - (1 - p)**2
        locus_coverages[locus] = coverage
        print(f"Locus {locus} combined frequency (p): {p:.4f}, Coverage: {coverage*100:.2f}%")
        
    # Calculate the total combined population coverage across all loci
    non_covered_fraction = 1.0
    for coverage in locus_coverages.values():
        non_covered_fraction *= (1 - coverage)
        
    total_coverage = 1 - non_covered_fraction
    
    return total_coverage

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python calculate_coverage.py <path_to_txt>", file=sys.stderr)
        sys.exit(1)

    # Hard-coded path to the CSV file
    csv_file_path = 'HLA_population_coverage.csv'
    txt_file_path = sys.argv[1]
    
    if not os.path.exists(txt_file_path):
        print(f"Error: The allele list file '{txt_file_path}' was not found.", file=sys.stderr)
        sys.exit(1)

    try:
        with open(txt_file_path, 'r') as f:
            alleles_to_analyze = [line.strip() for line in f if line.strip()]
    except Exception as e:
        print(f"Error reading allele list file: {e}", file=sys.stderr)
        sys.exit(1)
    
    if not alleles_to_analyze:
        print("The allele list file is empty.", file=sys.stderr)
        sys.exit(1)

    total_coverage_percentage = calculate_population_coverage(alleles_to_analyze, csv_file_path) * 100
    
    print(f"\nTotal combined population coverage for the given alleles: {total_coverage_percentage:.2f}%")