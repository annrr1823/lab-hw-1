import numpy as np
from collections import defaultdict, deque
import heapq

def complement(base):
    """Return the complementary nucleotide."""
    comp_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return comp_dict.get(base, base)

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    return ''.join(complement(base) for base in reversed(seq))

def find_kmers(seq, k=100):
    """Find all k-mers in a sequence with their positions."""
    kmers = defaultdict(list)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        kmers[kmer].append(i)
    return kmers

def find_anchors(query, reference, k=5, min_match_len=3, max_gap=500):
    """
    Find anchor points between query and reference sequences.
    Returns a list of tuples (q_start, q_end, r_start, r_end, is_reverse)
    """
    # Find k-mers in reference
    ref_kmers = find_kmers(reference, k)
    
    # Also consider reverse complement
    rev_query = reverse_complement(query)
    
    anchors = []
    
    # Find matching k-mers in forward direction
    for i in range(len(query) - k + 1):
        kmer = query[i:i+k]
        if kmer in ref_kmers:
            for j in ref_kmers[kmer]:
                # Extend match in both directions
                q_start, r_start = i, j
                q_end, r_end = i + k, j + k
                
                # Extend left
                ql, rl = q_start - 1, r_start - 1
                while ql >= 0 and rl >= 0 and query[ql] == reference[rl]:
                    ql -= 1
                    rl -= 1
                ql += 1
                rl += 1
                
                # Extend right
                qr, rr = q_end, r_end
                while qr < len(query) and rr < len(reference) and query[qr] == reference[rr]:
                    qr += 1
                    rr += 1
                
                # Calculate extended match length
                match_len = qr - ql
                
                # Add if meeting minimum length requirement
                if match_len >= min_match_len:
                    anchors.append((ql, qr, rl, rr, False))
    
    # Find matching k-mers for reverse complement
    for i in range(len(rev_query) - k + 1):
        kmer = rev_query[i:i+k]
        if kmer in ref_kmers:
            for j in ref_kmers[kmer]:
                # For reverse matches, we need to translate coordinates back to original query
                orig_i = len(query) - (i + k)
                
                # Extend match in both directions
                rev_q_start, r_start = i, j
                rev_q_end, r_end = i + k, j + k
                
                # Extend left
                ql, rl = rev_q_start - 1, r_start - 1
                while ql >= 0 and rl >= 0 and rev_query[ql] == reference[rl]:
                    ql -= 1
                    rl -= 1
                ql += 1
                rl += 1
                
                # Extend right
                qr, rr = rev_q_end, r_end
                while qr < len(rev_query) and rr < len(reference) and rev_query[qr] == reference[rr]:
                    qr += 1
                    rr += 1
                
                # Calculate match length
                match_len = qr - ql
                
                # Convert back to original query coordinates (reversed)
                q_start = len(query) - qr
                q_end = len(query) - ql
                
                # Add if meeting minimum length requirement
                if match_len >= min_match_len:
                    anchors.append((q_start, q_end, rl, rr, True))
    
    # Filter overlapping anchors to keep only the longest ones
    anchors.sort(key=lambda x: (x[0], x[2]))  # Sort by query_start and ref_start
    
    filtered_anchors = []
    for anchor in anchors:
        if not filtered_anchors:
            filtered_anchors.append(anchor)
            continue
            
        last = filtered_anchors[-1]
        # Check if current anchor overlaps with last anchor
        if anchor[0] >= last[1] or anchor[2] >= last[3]:
            filtered_anchors.append(anchor)
        else:
            # If overlapping, keep the longer one
            curr_len = anchor[1] - anchor[0]
            last_len = last[1] - last[0]
            if curr_len > last_len:
                filtered_anchors[-1] = anchor
    
    return filtered_anchors

def score_alignment(q_start1, q_end1, r_start1, r_end1, is_rev1,
                   q_start2, q_end2, r_start2, r_end2, is_rev2,
                   gap_penalty=0.5, inversion_penalty=25):
    """
    Score the alignment transition between two anchors.
    Lower score is better.
    """
    # Basic check that the second anchor comes after the first in query
    if q_start2 < q_end1:
        return float('inf')  # Invalid transition
    
    # Calculate gap sizes
    q_gap = q_start2 - q_end1
    
    # Calculate reference gap size accounting for direction
    if is_rev1 == is_rev2:
        # Same direction
        if is_rev1:  # Both reverse
            r_gap = r_start1 - r_end2 if r_start1 > r_end2 else r_start2 - r_end1
        else:  # Both forward
            r_gap = r_start2 - r_end1 if r_start2 > r_end1 else r_start1 - r_end2
    else:
        # Direction change, apply inversion penalty
        r_gap = abs(r_start2 - r_start1) + abs(r_end2 - r_end1)
    
    # Penalize inversions
    direction_penalty = inversion_penalty if is_rev1 != is_rev2 else 0
    
    # Calculate score (lower is better)
    # We want to minimize gaps and prefer consistent directions
    score = abs(q_gap - r_gap) * gap_penalty + direction_penalty
    
    return score

def chain_anchors(anchors, max_gap=200):
    """
    Chain anchors together to form the best alignment path.
    Uses a variant of the longest path algorithm in a DAG.
    """
    if not anchors:
        return []
    
    # Sort anchors by position in query
    anchors.sort(key=lambda x: x[0])
    
    n = len(anchors)
    # Initialize scores and backtrack pointers
    scores = [0] * n
    prev = [-1] * n
    
    # Calculate best score for each anchor
    for i in range(n):
        print("i:",i,"in n:",n)
        q_start_i, q_end_i, r_start_i, r_end_i, is_rev_i = anchors[i]
        
        # Base score for this anchor (match length)
        scores[i] = q_end_i - q_start_i
        
        # Try connecting to previous anchors
        for j in range(i):
            # print("j:",j,"in i:",i)
            q_start_j, q_end_j, r_start_j, r_end_j, is_rev_j = anchors[j]
            
            # Skip if gap is too large
            if q_start_i - q_end_j > max_gap:
                continue
                
            # Calculate transition score
            alignment_score = score_alignment(
                q_start_j, q_end_j, r_start_j, r_end_j, is_rev_j,
                q_start_i, q_end_i, r_start_i, r_end_i, is_rev_i
            )
            
            # Skip invalid transitions
            if alignment_score == float('inf'):
                continue
                
            # Calculate new score if we connect j to i
            # We want to maximize the score, so we subtract the alignment score (lower is better)
            new_score = scores[j] + (q_end_i - q_start_i) - alignment_score
            
            if new_score > scores[i]:
                scores[i] = new_score
                prev[i] = j
    
    # Find the anchor with the best score
    best_idx = max(range(n), key=lambda i: scores[i])
    
    # Reconstruct the path
    path = []
    while best_idx != -1:
        path.append(anchors[best_idx])
        best_idx = prev[best_idx]
    
    # Reverse the path to get it in the correct order
    path.reverse()
    
    return path


def create_contiguous_alignments(alignments, query, reference):
    """
    Create contiguous, non-overlapping alignments covering the entire query sequence.
    For gaps between anchored alignments, perform local alignment to find the best match.
    Handles both forward and reverse complement alignments.
    Internal processing tracks orientation but final output only includes coordinates.
    """
    if not alignments:
        return [(0, len(query), 0, len(reference))]

    # Sort alignments by query start position
    alignments.sort(key=lambda x: x[0])
    non_overlapping = []
    prev_end = 0
    prev_r_end = 0
    prev_is_rev = False  # Track orientation of previous alignment
    
    # Internal list that tracks orientation
    internal_alignments = []

    for i, align in enumerate(alignments):
        # We need to know if this is a reverse alignment
        if len(align) == 5:  # Includes is_rev flag
            q_start, q_end, r_start, r_end, is_rev = align
        else:  # No explicit is_rev flag, assume forward
            q_start, q_end, r_start, r_end = align
            is_rev = False
            
        # If there's a gap between alignments in query
        if q_start-1 > prev_end + 1:
            gap_query = query[prev_end+1:q_start-1]
            if len(gap_query) > 0:
                best_score = float('inf')
                best_r_start = 0
                best_is_rev = False
                
                # Define search window around expected region in reference
                ref_window_start = max(0, min(r_start, prev_r_end) - 1000)
                ref_window_end = min(len(reference) - len(gap_query), max(r_end, prev_r_end) + 1000)
                step = max(int((ref_window_end-ref_window_start)/5000), 1)
                print("i:",i,"ref_window_start:", ref_window_start, "ref_window_end:", ref_window_end)
                
                # First try forward orientation
                for j in range(ref_window_start, ref_window_end, step):
                    print("forward j:", j)
                    r_segment = reference[j:j+len(gap_query)]
                    if len(r_segment) < len(gap_query):
                        continue
                    score = levenshtein_distance(gap_query, r_segment)
                    if score < best_score:
                        best_score = score
                        best_r_start = j
                        best_is_rev = False
                
                # Then try reverse complement orientation
                rev_gap_query = reverse_complement(gap_query)
                for j in range(ref_window_start, ref_window_end, step):
                    print("reverse j:", j)
                    r_segment = reference[j:j+len(rev_gap_query)]
                    if len(r_segment) < len(rev_gap_query):
                        continue
                    score = levenshtein_distance(rev_gap_query, r_segment)
                    if score < best_score:
                        best_score = score
                        best_r_start = j
                        best_is_rev = True
                
                # Add the gap alignment with orientation (for internal use)
                internal_alignments.append((prev_end+1, q_start-1, best_r_start, 
                                        best_r_start + len(gap_query), best_is_rev))
        
        # Handle the current alignment segment
        adj_start = max(q_start, prev_end)
        if adj_start < q_end:
            if adj_start > q_start:
                # Need to adjust reference position proportionally
                trim_ratio = (adj_start - q_start) / (q_end - q_start) if q_end > q_start else 0
                if is_rev:
                    # For reverse alignments, we need to adjust from the other end
                    adj_r_end = r_end - int((r_end - r_start) * trim_ratio)
                    internal_alignments.append((adj_start, q_end, r_start, adj_r_end, is_rev))
                else:
                    # Forward alignment adjustment
                    adj_r_start = r_start + int((r_end - r_start) * trim_ratio)
                    internal_alignments.append((adj_start, q_end, adj_r_start, r_end, is_rev))
            else:
                # No adjustment needed
                internal_alignments.append((q_start, q_end, r_start, r_end, is_rev))
        
        # Update tracking variables
        prev_end = max(prev_end, q_end)
        if is_rev:
            # For reverse alignments, prev_r_end is actually r_start
            prev_r_end = r_start
        else:
            prev_r_end = r_end
        prev_is_rev = is_rev
    
    # Handle final gap if exists
    if prev_end < len(query)-1:
        final_gap = query[prev_end:len(query)-1]
        if len(final_gap) > 0:
            best_score = float('inf')
            best_r_start = 0
            best_is_rev = False
            
            # Search window based on the last alignment
            ref_window_start = max(0, prev_r_end - 1000)
            ref_window_end = min(len(reference) - len(final_gap), prev_r_end + 1000)
            step = max(int((ref_window_end-ref_window_start)/5000), 1)
            print("final ref_window_start:", ref_window_start, "final ref_window_end:", ref_window_end)
            
            # Try forward orientation
            for j in range(ref_window_start, ref_window_end, step):
                print("final forward j:", j)
                r_segment = reference[j:j+len(final_gap)]
                if len(r_segment) < len(final_gap):
                    continue
                score = levenshtein_distance(final_gap, r_segment)
                if score < best_score:
                    best_score = score
                    best_r_start = j
                    best_is_rev = False
            
            # Try reverse complement orientation
            rev_final_gap = reverse_complement(final_gap)
            for j in range(ref_window_start, ref_window_end, step):
                print("final reverse j:", j)
                r_segment = reference[j:j+len(rev_final_gap)]
                if len(r_segment) < len(rev_final_gap):
                    continue
                score = levenshtein_distance(rev_final_gap, r_segment)
                if score < best_score:
                    best_score = score
                    best_r_start = j
                    best_is_rev = True
            
            # Add the final gap alignment with orientation (for internal use)
            internal_alignments.append((prev_end, len(query), best_r_start, 
                                       best_r_start + len(final_gap), best_is_rev))
    
    # Convert internal alignments to output format without is_rev flag
    for align in internal_alignments:
        q_start, q_end, r_start, r_end = align[0:4]
        non_overlapping.append((q_start, q_end, r_start, r_end))
    
    return non_overlapping

# Update the align_dna_sequences function to process is_rev internally but return only coordinates
def align_dna_sequences(query, reference):
    """
    Align query sequence to reference sequence considering all mutation types.
    Returns a list of contiguous alignment regions as (q_start, q_end, r_start, r_end) tuples.
    Internally tracks orientation (forward/reverse) but final output excludes this information.
    """
    # Find anchor points
    print("start find anchors")
    anchors = find_anchors(query, reference)
    print("end find anchors")
    
    # Chain anchors to find best alignment path
    print("start chain anchors")
    best_path = chain_anchors(anchors)
    print("end chain anchors")
    
    # Keep is_rev flag when converting path to alignment format (for internal use)
    alignments = [(q_start, q_end, r_start, r_end, is_rev) for q_start, q_end, r_start, r_end, is_rev in best_path]
    
    # Merge close alignments using orientation info
    print("start merge")
    merged_alignments = merge_close_alignments(alignments)
    print("end merge")

    # Convert to contiguous non-overlapping segments
    # The function now returns tuples without is_rev flag
    print("start contiguous")
    contiguous_alignments = create_contiguous_alignments(merged_alignments, query, reference)
    print("end contiguous")
    
    return contiguous_alignments

# Update the merge_close_alignments function to handle the is_rev flag internally
def merge_close_alignments(alignments, max_gap=7):
    """
    Merge alignments that are close to each other, respecting orientation.
    Processes alignments with orientation info but only returns coordinate tuples.
    """
    if not alignments:
        return []
    
    # Sort alignments by query start position
    alignments.sort(key=lambda x: x[0])
    
    merged = []
    current = list(alignments[0])
    
    for alignment in alignments[1:]:
        if len(alignment) == 5:  # Includes is_rev flag
            q_start, q_end, r_start, r_end, is_rev = alignment
        else:
            q_start, q_end, r_start, r_end = alignment
            is_rev = False
        
        # Only merge if orientation matches
        if len(current) >= 5 and current[4] == is_rev:
            # Check if the current alignment can be merged with the next one
            if (q_start - current[1] <= max_gap):
                # Same orientation, check reference coordinates appropriately
                if not is_rev and r_start >= current[3] and r_start - current[3] <= max_gap:
                    # Forward direction - extend right
                    current[1] = max(current[1], q_end)
                    current[3] = max(current[3], r_end)
                elif is_rev and r_end <= current[2] and current[2] - r_end <= max_gap:
                    # Reverse direction - extend left in reference
                    current[1] = max(current[1], q_end)
                    current[2] = min(current[2], r_start)
                else:
                    # Cannot merge, add current to merged list and start new current
                    merged.append(tuple(current))
                    current = list(alignment)
            else:
                # Cannot merge, add current to merged list and start new current
                merged.append(tuple(current))
                current = list(alignment)
        else:
            # Different orientation or missing is_rev, don't merge
            merged.append(tuple(current))
            current = list(alignment)
    
    # Add the last alignment
    merged.append(tuple(current))
    
    # We keep the is_rev flag here for internal processing
    # The final output is handled by create_contiguous_alignments
    
    return merged

def levenshtein_distance(s1, s2):
    """
    Calculate the Levenshtein distance between two strings.
    This measures the minimum number of single-character edits needed to change one string into another.
    """
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)
    
    if len(s2) == 0:
        return len(s1)
    
    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            # Calculate insertions, deletions and substitutions
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            
            # Get minimum of the three operations
            current_row.append(min(insertions, deletions, substitutions))
        
        previous_row = current_row
    
    return previous_row[-1]



def load_sequence_from_file(filepath):
    with open(filepath, 'r') as f:
        return f.read().replace('\n', '').strip().upper()


def main():
    reference = load_sequence_from_file(r"D:\作业\大二下\算法\lab 2\refrence.txt")
    query = load_sequence_from_file(r"D:\作业\大二下\算法\lab 2\query.txt")

    alignments = align_dna_sequences(query, reference)

    print("匹配区间元组列表:")
    print(", ".join(str(block) for block in alignments))

    print("len of reference:", len(reference))
    print("len of query:", len(query))

if __name__ == "__main__":
    main()