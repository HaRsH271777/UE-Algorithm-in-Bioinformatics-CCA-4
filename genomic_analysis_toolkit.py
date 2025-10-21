"""
genomic_analysis_toolkit.py
A comprehensive Python module for genomic sequence analysis.

This script provides functions for:
- Part A: GC Content, GC Skew, and k-mer frequency analysis.
- Part B: Hamming distance, exact pattern matching, and fuzzy pattern matching.

Author: [Your Name Here]
Date: 21-Oct-2025
"""

import collections
import re
import math # Needed for CpG calculation

# --- Part A: GC Content and Sequence Analysis ---

def calculate_gc_content(seq):
    """
    Calculates the overall GC content percentage of a DNA sequence.
    
    Ambiguous bases (like 'N') are ignored in the total count,
    leading to a more accurate GC percentage of the known sequence.

    Args:
        seq (str): A DNA sequence string.

    Returns:
        float: The GC content as a percentage (0.0 to 100.0).
               Returns 0.0 for empty sequences or sequences with
               no 'A', 'T', 'C', or 'G' bases.
               
    Performance:
        Time: O(N), where N is the length of the sequence.
              The string .count() method scans the string.
        Space: O(1), as we only store a few count variables.
    """
    if not seq:
        return 0.0

    seq = seq.upper()
    
    g_count = seq.count('G')
    c_count = seq.count('C')
    
    # Total valid bases for this calculation
    total_known_bases = g_count + c_count + seq.count('A') + seq.count('T')
    
    if total_known_bases == 0:
        return 0.0

    gc_percent = (g_count + c_count) / total_known_bases * 100
    return gc_percent

def gc_content_sliding_window(seq, window_size, step):
    """
    Calculates GC content in a sliding window across a sequence.

    Args:
        seq (str): The DNA sequence.
        window_size (int): The size of the sliding window.
        step (int): The number of bases to move the window forward.

    Returns:
        list[float]: A list of GC percentages for each window.
        
    Performance:
        Time: O(N*M), where N is seq length and M is window_size.
        Space: O(N/step), to store the results list.
    """
    seq = seq.upper()
    results = []
    
    for i in range(0, len(seq) - window_size + 1, step):
        window_seq = seq[i:i+window_size]
        window_gc = calculate_gc_content(window_seq)
        results.append(window_gc)
        
    return results

def calculate_gc_skew_sliding_window(seq, window_size, step):
    """
    Calculates GC skew ((G-C)/(G+C)) in a sliding window.
    GC skew is useful for identifying replication origins.

    Args:
        seq (str): The DNA sequence.
        window_size (int): The size of the sliding window.
        step (int): The number of bases to move the window forward.

    Returns:
        list[float]: A list of GC skew values for each window.
                     Returns 0.0 for a window if (G+C) is 0.
                     
    Performance:
        Time: O(N*M), same logic as gc_content_sliding_window.
        Space: O(N/step), to store the results list.
    """
    seq = seq.upper()
    results = []
    
    for i in range(0, len(seq) - window_size + 1, step):
        window_seq = seq[i:i+window_size]
        
        g_count = window_seq.count('G')
        c_count = window_seq.count('C')
        
        gc_sum = g_count + c_count
        
        if gc_sum == 0:
            results.append(0.0)
        else:
            skew = (g_count - c_count) / gc_sum
            results.append(skew)
            
    return results

def calculate_kmer_freq(seq, k):
    """
    Calculates the frequency of all k-mers (e.g., dinucleotides k=2).
    Ignores any k-mers containing ambiguous bases ('N').

    Args:
        seq (str): The DNA sequence.
        k (int): The length of the k-mer (e.g., 2 for dinucleotide).

    Returns:
        collections.Counter: A Counter object mapping k-mers to their counts.
        
    Performance:
        Time: O(N*k), where N is sequence length.
        Space: O(4^k), in the worst case, to store the Counter dict.
    """
    seq = seq.upper()
    kmer_counts = collections.Counter()
    
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        
        if 'N' not in kmer:
            kmer_counts[kmer] += 1
            
    return kmer_counts

def find_cpg_islands(seq, window_size=200, gc_threshold=50.0, oe_threshold=0.6):
    """
    Identifies potential CpG islands in a genomic sequence.
    
    A region is considered a CpG island if it meets:
    1. Length is at least `window_size`.
    2. GC Content is above `gc_threshold`.
    3. Observed-to-Expected CpG ratio is above `oe_threshold`.

    Args:
        seq (str): The DNA sequence.
        window_size (int): The minimum length to check.
        gc_threshold (float): The minimum GC content (e.g., 50.0).
        oe_threshold (float): The minimum Observed/Expected CpG ratio (e.g., 0.6).

    Returns:
        list[tuple(int, int)]: A list of (start, end) tuples for
                               identified CpG islands.
                               
    Performance:
        Time: O(N*M), where N is seq length and M is window_size.
        Space: O(K), where K is the number of islands found.
    """
    seq = seq.upper()
    islands = []
    
    for i in range(len(seq) - window_size + 1):
        window = seq[i:i+window_size]
        
        gc_content = calculate_gc_content(window)
        if gc_content < gc_threshold:
            continue
            
        c_count = window.count('C')
        g_count = window.count('G')
        cg_count = window.count('CG')
        
        if c_count == 0 or g_count == 0:
            continue # Avoid division by zero if no C or G

        # Calculate expected value
        expected_cpg = (c_count * g_count) / window_size
        
        if expected_cpg == 0:
            continue
            
        observed_expected_ratio = cg_count / expected_cpg
            
        if observed_expected_ratio > oe_threshold:
            islands.append((i, i + window_size))

    return islands

# --- Part B: Pattern Matching and Motif Discovery ---

def hamming_distance(seq1, seq2):
    """
    Calculates the Hamming distance between two sequences.
    The distance is the number of positions at which the
    corresponding symbols are different.

    Args:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.

    Returns:
        int: The Hamming distance.

    Raises:
        ValueError: If the sequences are of different lengths.
        
    Performance:
        Time: O(N), where N is the length of the sequences.
        Space: O(1).
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length for Hamming distance.")
        
    distance = 0
    # zip stops at the shortest sequence, but we already checked lengths
    for base1, base2 in zip(seq1.upper(), seq2.upper()):
        if base1 != base2:
            distance += 1
            
    return distance

def find_exact_pattern(seq, pattern):
    """
    Finds all starting indices of an exact pattern, including overlaps.

    Args:
        seq (str): The sequence to search in.
        pattern (str): The pattern to search for.

    Returns:
        list[int]: A list of 0-based starting indices.
        
    Performance:
        Time: O(N*M) in the worst case (e.g., finding 'AAA' in 'AAAAA').
        Space: O(K), where K is the number of matches found.
    """
    seq = seq.upper()
    pattern = pattern.upper()
    
    pattern_escaped = re.escape(pattern)
    
    try:
        # A positive lookahead `(?=...)` allows finding overlapping matches
        indices = [m.start() for m in re.finditer(f'(?={pattern_escaped})', seq)]
    except re.error:
        print(f"Regex error with pattern: {pattern}")
        return []
        
    return indices

def find_fuzzy_pattern(seq, pattern, max_mismatch):
    """
    Finds all starting indices of a pattern with up to N mismatches.

    Args:
        seq (str): The sequence to search in.
        pattern (str): The pattern to search for.
        max_mismatch (int): The maximum number of allowed mismatches
                            (Hamming distance).

    Returns:
        list[int]: A list of 0-based starting indices.
        
    Performance:
        Time: O(N*M), where N is seq length and M is pattern length.
        Space: O(K), where K is the number of matches found.
    """
    seq = seq.upper()
    pattern = pattern.upper()
    indices = []
    pattern_len = len(pattern)
    
    if pattern_len > len(seq):
        return [] # Pattern can't be found if it's longer
        
    for i in range(len(seq) - pattern_len + 1):
        window = seq[i:i+pattern_len]
        
        try:
            dist = hamming_distance(window, pattern)
            if dist <= max_mismatch:
                indices.append(i)
        except ValueError:
            continue 
            
    return indices


# --- Main Test Block ---
if __name__ == "__main__":
    
    print("--- Genomic Analysis Toolkit Tests ---")

    # --- Part A Test Cases ---
    print("\n--- Part A: Sequence Analysis ---")
    seq_a = "AGCTATAGCGCCGATTAGCATGGTATAGTAGAATTC"
    seq_b = "NNNNNNNNNN"
    seq_c = "ATATATATAT"
    seq_d = "" # Empty sequence edge case
    
    print(f"GC Content (seq_a): {calculate_gc_content(seq_a):.2f}%")
    print(f"GC Content (seq_b 'NNN'): {calculate_gc_content(seq_b):.2f}%")
    win_gcs = gc_content_sliding_window(seq_a, 10, 5)
    print(f"Sliding Window GC (10, 5): {[round(g, 2) for g in win_gcs]}")
    win_skews = calculate_gc_skew_sliding_window(seq_a, 10, 5)
    print(f"Sliding Window Skew (10, 5): {[round(s, 2) for s in win_skews]}")
    print(f"Dinucleotide Freq (seq_a): {calculate_kmer_freq(seq_a, 2).most_common(3)}")
    cpg_test_seq = "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG" * 5
    no_cpg_test_seq = "CCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGG" * 5
    print(f"CpG Islands (positive test): {find_cpg_islands(cpg_test_seq, 50, 50, 0.6)}")
    print(f"CpG Islands (negative test): {find_cpg_islands(no_cpg_test_seq, 50, 50, 0.6)}")

    # --- Part B Test Cases ---
    print("\n--- Part B: Pattern Matching ---")
    
    # Q3: Hamming Distance
    try:
        dist = hamming_distance("GATTACA", "GATTTCA")
        print(f"Hamming ('GATTACA', 'GATTTCA'): {dist}") # Should be 1
        print("Testing Hamming with different lengths:")
        dist_fail = hamming_distance("AAA", "AAAA") # Should raise error
    except ValueError as e:
        print(f"Hamming Error (expected): {e}")

    # Q4: Exact Pattern
    pattern_seq = "ATATATACATAT"
    pattern = "ATA"
    print(f"Exact '{pattern}' in '{pattern_seq}': {find_exact_pattern(pattern_seq, pattern)}") # [0, 2, 4, 8]
    print(f"Exact 'AAA' in 'AAAAA': {find_exact_pattern('AAAAA', 'AAA')}") # [0, 1, 2]

    # Q4: Fuzzy Pattern
    fuzzy_seq = "GATTACAT"
    fuzzy_pattern = "GATTTCA"
    print(f"Fuzzy '{fuzzy_pattern}' in '{fuzzy_seq}' (<=1 mismatch): {find_fuzzy_pattern(fuzzy_seq, fuzzy_pattern, 1)}") # [0]
    print(f"Fuzzy '{fuzzy_pattern}' in '{fuzzy_seq}' (<=0 mismatch): {find_fuzzy_pattern(fuzzy_seq, fuzzy_pattern, 0)}") # []