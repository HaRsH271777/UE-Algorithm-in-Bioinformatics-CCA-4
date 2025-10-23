"""
UE - Algorithm in Bioinformatics
CCA3 
Name: Harshvardhan Salunke
TY B.TECH CSE Ai & DS
PRN - 1032230341

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
        Space: O(1).
    """
    if not seq:
        return 0.0

    seq = seq.upper()
    
    g_count = seq.count('G')
    c_count = seq.count('C')
    
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
        Space: O(N/step).
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

    Args:
        seq (str): The DNA sequence.
        window_size (int): The size of the sliding window.
        step (int): The number of bases to move the window forward.

    Returns:
        list[float]: A list of GC skew values for each window.
                     Returns 0.0 for a window if (G+C) is 0.
                     
    Performance:
        Time: O(N*M), same logic as gc_content_sliding_window.
        Space: O(N/step).
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
        Space: O(4^k).
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
            continue

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
        Time: O(N*M) in the worst case.
        Space: O(K), where K is the number of matches found.
    """
    seq = seq.upper()
    pattern = pattern.upper()
    
    pattern_escaped = re.escape(pattern)
    
    try:
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
        max_mismatch (int): The maximum number of allowed mismatches.

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
        return [] 
        
    for i in range(len(seq) - pattern_len + 1):
        window = seq[i:i+pattern_len]
        
        try:
            dist = hamming_distance(window, pattern)
            if dist <= max_mismatch:
                indices.append(i)
        except ValueError:
            continue 
            
    return indices

# --- Part C: Protein Sequence Analysis ---

# Standard DNA codon table. '_' represents a STOP codon.
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def reverse_complement(seq):
    """
    Returns the reverse complement of a DNA sequence.
    Handles 'N' bases by complementing them to 'N'.

    Args:
        seq (str): The DNA sequence.

    Returns:
        str: The reverse complemented sequence.
        
    Performance:
        Time: O(N), where N is the sequence length.
        Space: O(N), to build the new list/string.
    """
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    
    # Map unknown characters to 'N'
    complement_seq = "".join([complement_map.get(base, 'N') for base in seq.upper()])
    reverse_comp = complement_seq[::-1] 
    
    return reverse_comp

def get_reading_frames(seq):
    """
    Gets all 6 reading frames for a DNA sequence (3 forward, 3 reverse).

    Args:
        seq (str): The DNA sequence.

    Returns:
        list[str]: A list of 6 strings, representing frames
                   +1, +2, +3, -1, -2, -3.
                   
    Performance:
        Time: O(N), dominated by the reverse_complement call.
        Space: O(N), to store the 6 frames.
    """
    seq = seq.upper()
    rev_comp = reverse_complement(seq)
    
    frames = []
    # Forward frames
    frames.append(seq)      # Frame +1
    frames.append(seq[1:])  # Frame +2
    frames.append(seq[2:])  # Frame +3
    
    # Reverse frames
    frames.append(rev_comp)      # Frame -1
    frames.append(rev_comp[1:])  # Frame -2
    frames.append(rev_comp[2:])  # Frame -3
    
    return frames

def translate(seq):
    """
    Translates a DNA sequence (representing a single reading frame)
    into a protein sequence.
    
    Translation stops at the first STOP codon.
    Ambiguous codons (containing 'N') are translated as 'X'.

    Args:
        seq (str): A DNA sequence.

    Returns:
        str: The translated amino acid sequence.
        
    Performance:
        Time: O(N), where N is the sequence length.
        Space: O(N/3), to store the resulting protein string.
    """
    seq = seq.upper()
    protein = []
    
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        
        if len(codon) < 3:
            break 
            
        amino_acid = CODON_TABLE.get(codon, 'X')
        
        if amino_acid == '_':
            break # Stop at first STOP codon
            
        protein.append(amino_acid)
        
    return "".join(protein)

def find_all_orfs(seq, start_codons=['ATG'], stop_codons=['TAA', 'TAG', 'TGA']):
    """
    Finds all potential protein-coding Open Reading Frames (ORFs)
    in all 6 reading frames.
    
    An ORF is defined as a sequence starting with a START codon
    and ending with a STOP codon, in the same frame.

    Args:
        seq (str): The DNA sequence.
        start_codons (list[str]): List of start codons.
        stop_codons (list[str]): List of stop codons.

    Returns:
        list[dict]: A list of dictionaries, where each dict contains:
                    'frame' (str): e.g., '+1', '-2'
                    'start' (int): 0-based start index *on the original seq*
                    'end' (int): 0-based end index *on the original seq*
                    'protein' (str): The translated protein sequence.
                    
    Performance:
        Time: O(N^2) in worst case (many starts, few stops).
              Closer to O(N) on typical genomic sequence.
        Space: O(K), where K is the number of ORFs found.
    """
    frames = get_reading_frames(seq)
    frame_names = ['+1', '+2', '+3', '-1', '-2', '-3']
    all_orfs = []
    seq_len = len(seq)
    
    for frame_idx, (frame_name, frame_seq) in enumerate(zip(frame_names, frames)):
        
        # Find all start codons in this frame
        start_indices_in_frame = []
        for start_codon in start_codons:
            for match in re.finditer(f'(?={re.escape(start_codon)})', frame_seq):
                start_indices_in_frame.append(match.start())
        
        start_indices_in_frame.sort()

        for start_pos in start_indices_in_frame:
            protein = []
            
            # Translate from this start codon
            for i in range(start_pos, len(frame_seq), 3):
                codon = frame_seq[i:i+3]
                
                if len(codon) < 3:
                    break 
                    
                amino_acid = CODON_TABLE.get(codon, 'X')
                
                if codon in stop_codons:
                    # Found a valid ORF!
                    protein_seq = "".join(protein)
                    
                    # Map frame coordinates back to original seq coordinates
                    frame_offset = frame_idx % 3
                    
                    if frame_name.startswith('+'):
                        # Forward frame
                        orf_start_on_seq = start_pos + frame_offset
                        orf_end_on_seq = i + 3 + frame_offset
                    else:
                        # Reverse frame
                        rev_start = start_pos + frame_offset
                        rev_end = i + 3 + frame_offset
                        orf_start_on_seq = seq_len - rev_end
                        orf_end_on_seq = seq_len - rev_start
                    
                    all_orfs.append({
                        'frame': frame_name,
                        'start': orf_start_on_seq,
                        'end': orf_end_on_seq,
                        'protein': protein_seq,
                        'length': len(protein_seq)
                    })
                    break # Stop translating this ORF
                
                if amino_acid == 'X':
                    break # Ambiguous codon, invalidates this ORF

                protein.append(amino_acid)
                
    return all_orfs

def find_longest_orf(seq, **kwargs):
    """
    Finds the longest ORF in a DNA sequence.
    
    Args:
        seq (str): The DNA sequence.
        **kwargs: Passed to find_all_orfs (e.g., start_codons)

    Returns:
        dict: The ORF dictionary for the longest protein, 
              or None if no ORFs are found.
    """
    all_orfs = find_all_orfs(seq, **kwargs)
    if not all_orfs:
        return None
        
    longest_orf = max(all_orfs, key=lambda orf: orf['length'])
    return longest_orf


# --- Main Test Block ---
if __name__ == "__main__":
    
    print("--- Genomic Analysis Toolkit Tests ---")

    # --- Part A Test Cases ---
    print("\n--- Part A: Sequence Analysis ---")
    seq_a = "AGCTATAGCGCCGATTAGCATGGTATAGTAGAATTC"
    print(f"GC Content (seq_a): {calculate_gc_content(seq_a):.2f}%")
    win_gcs = gc_content_sliding_window(seq_a, 10, 5)
    print(f"Sliding Window GC (10, 5): {[round(g, 2) for g in win_gcs]}")
    win_skews = calculate_gc_skew_sliding_window(seq_a, 10, 5)
    print(f"Sliding Window Skew (10, 5): {[round(s, 2) for s in win_skews]}")
    print(f"Dinucleotide Freq (seq_a): {calculate_kmer_freq(seq_a, 2).most_common(3)}")

    # --- Part B Test Cases ---
    print("\n--- Part B: Pattern Matching ---")
    try:
        dist = hamming_distance("GATTACA", "GATTTCA")
        print(f"Hamming ('GATTACA', 'GATTTCA'): {dist}")
        dist_fail = hamming_distance("AAA", "AAAA")
    except ValueError as e:
        print(f"Hamming Error (expected): {e}")

    pattern_seq = "ATATATACATAT"
    pattern = "ATA"
    print(f"Exact '{pattern}' in '{pattern_seq}': {find_exact_pattern(pattern_seq, pattern)}")
    fuzzy_seq = "GATTACAT"
    fuzzy_pattern = "GATTTCA"
    print(f"Fuzzy '{fuzzy_pattern}' in '{fuzzy_seq}' (<=1): {find_fuzzy_pattern(fuzzy_seq, fuzzy_pattern, 1)}")

    # --- Part C Test Cases ---
    print("\n--- Part C: Protein Analysis ---")
    
    print(f"Reverse Complement ('ATGCGN'): {reverse_complement('ATGCGN')}") # NCGNNGCAT
    print(f"Translate ('ATGGGTTAA'): {translate('ATGGGTTAA')}") # MG (stops at TAA)
    
    # Test sequence with ORFs in forward and reverse frames
    # +1: ATG GCA TAA  (MA)
    # -2: (rev-comp) TT ATG GGG TAG AA (MG)
    orf_test_seq = "AAATGGCATAATTCACTATGgggTAGAATT" # 'ggg' is part of -2 ORF
    
    print(f"\nTesting ORF finding on: {orf_test_seq}")
    all_orfs = find_all_orfs(orf_test_seq)
    print("All ORFs found:")
    for orf in all_orfs:
        print(f"  - {orf['frame']} (Start: {orf['start']}, End: {orf['end']}): {orf['protein']}")
        
    longest = find_longest_orf(orf_test_seq)
    if longest:
        print(f"\nLongest ORF: {longest['protein']} (from frame {longest['frame']})")
    else:
        print("\nNo ORFs found.")