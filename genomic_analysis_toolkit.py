"""
UE - Algorithm in Bioinformatics
CCA4 

Name: Harshvardhan Salunke
TY B.TECH CSE Ai & DS
PRN - 1032230341

"""

import collections
import re
import math # Needed for CpG calculation

# --- Part A---

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
              This is because for each window (approx N/step), 
              we call calculate_gc_content, which is O(M).
              *Optimization Note: This could be O(N) by using a
              sliding window with a deque, adding/removing one 
              base at a time, but this implementation is simpler.*
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
        Time: O(N*k), where N is sequence length. Slicing (seq[i:i+k])
              can take O(k) time, and we do it N-k times.
        Space: O(4^k), in the worst case, to store the Counter dict.
               In practice, it's O(U) where U is unique k-mers found.
    """
    seq = seq.upper()
    kmer_counts = collections.Counter()
    
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        
        # A simple check for 'N'. A more robust check
        # would use a set of allowed characters.
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

    Argumentss:
        seq (str): The DNA sequence.
        window_size (int): The minimum length to check.
        gc_threshold (float): The minimum GC content (e.g., 50.0).
        oe_threshold (float): The minimum Observed/Expected CpG ratio (e.g., 0.6).

    Returns:
        list[tuple(int, int)]: A list of (start, end) tuples for
                               identified CpG islands.
                               
    Performance:
        Time: O(N*M), where N is seq length and M is window_size.
              We loop N-M times, and each window slice involves
              O(M) work for .count() operations.
        Space: O(K), where K is the number of islands found.
    """
    seq = seq.upper()
    islands = []
    
    for i in range(len(seq) - window_size + 1):
        window = seq[i:i+window_size]
        
        # 1. Check GC Content
        gc_content = calculate_gc_content(window)
        if gc_content < gc_threshold:
            continue
            
        # 2. Check Observed/Expected CpG
        # We need to be careful about division by zero
        
        c_count = window.count('C')
        g_count = window.count('G')
        cg_count = window.count('CG')
        
        # Calculate expected value
        # expected = (num_C * num_G) / window_length
        expected_cpg = (c_count * g_count) / window_size
        
        if expected_cpg == 0:
            # If expected is 0, ratio is undefined.
            # If observed (cg_count) is also 0, it's not an island.
            # If observed > 0, this is a "perfect" island (infinite ratio)
            if cg_count > 0:
                 observed_expected_ratio = float('inf')
            else:
                 observed_expected_ratio = 0.0
        else:
            observed_expected_ratio = cg_count / expected_cpg
            
        # 3. Check O/E threshold
        if observed_expected_ratio > oe_threshold:
            islands.append((i, i + window_size))


    return islands


if __name__ == "__main__":
    
    print("--- Genomic Analysis Toolkit Tests ---")

    # --- Part A Test Cases ---
    print("\n--- Part A: ---")
    seq_a = "AGCTATAGCGCCGATTAGCATGGTATAGTAGAATTC"
    seq_b = "NNNNNNNNNN"
    seq_c = "ATATATATAT"
    seq_d = "" # Empty sequence edge case
    
    # Q1: GC Content
    print(f"GC Content (seq_a): {calculate_gc_content(seq_a):.2f}%") # Should be almost 44.4%
    print(f"GC Content (seq_b 'NNN'): {calculate_gc_content(seq_b):.2f}%") # Should be 0.0
    print(f"GC Content (seq_c 'ATAT'): {calculate_gc_content(seq_c):.2f}%") # Should be 0.0
    print(f"GC Content (seq_d ''): {calculate_gc_content(seq_d):.2f}%") # Should be 0.0

    # Q1: GC Sliding Window
    win_gcs = gc_content_sliding_window(seq_a, 10, 5)
    print(f"Sliding Window GC (10, 5): {[round(g, 2) for g in win_gcs]}")
    
    # Q1: GC Skew
    # seq_a[0:10] = 'AGCTATAGCG' (3G, 2C) -> (3-2)/(3+2) = 1/5 = 0.2
    # seq_a[5:15] = 'TAGCGCCGAT' (4G, 2C) -> (4-2)/(4+2) = 2/6 = 0.33
    win_skews = calculate_gc_skew_sliding_window(seq_a, 10, 5)
    print(f"Sliding Window Skew (10, 5): {[round(s, 2) for s in win_skews]}")
    
    # Q2: k-mer Frequencies
    print(f"Dinucleotide Freq (seq_a): {calculate_kmer_freq(seq_a, 2).most_common(3)}")
    print(f"Trinucleotide Freq (seq_a): {calculate_kmer_freq(seq_a, 3).most_common(3)}")

    # Q2: CpG Islands
    # This sequence is high GC and high CpG (CGCGC...)
    cpg_test_seq = "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG" * 5
    # This sequence is high GC but low CpG (CCC...GGG...)
    no_cpg_test_seq = "CCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGG" * 5
    
    print(f"CpG Islands (positive test): {find_cpg_islands(cpg_test_seq, 50, 50, 0.6)}")
    print(f"CpG Islands (negative test): {find_cpg_islands(no_cpg_test_seq, 50, 50, 0.6)}")

