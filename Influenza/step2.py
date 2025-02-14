DNA_dict = {
    'G': {'G': 2, 'C': -2, 'A': -2, 'T': -2, 'N': -2, '-': -4},
    'C': {'G': -2, 'C': 2, 'A': -2, 'T': -2, 'N': -2, '-': -4},
    'A': {'G': -2, 'C': -2, 'A': 2, 'T': -2, 'N': -2, '-': -4},
    'T': {'G': -2, 'C': -2, 'A': -2, 'T': 2, 'N': -2, '-': -4},
    'N': {'G': -2, 'C': -2, 'A': -2, 'T': -2, 'N': 2, '-': -4},
    '-': {'G': -4, 'C': -4, 'A': -4, 'T': -4, 'N': -4, '-': -4}
}

import numpy as np

def SequenceAlign(seq1, seq2, DNA_dict):
    len1, len2 = len(seq1), len(seq2)
    score_matrix = np.zeros((len1 + 1, len2 + 1))
    traceback_matrix = np.zeros((len1 + 1, len2 + 1), dtype=int)

    for i in range(1, len1 + 1):
        score_matrix[i, 0] = score_matrix[i - 1, 0] + DNA_dict[seq1[i - 1]]['-']
    for j in range(1, len2 + 1):
        score_matrix[0, j] = score_matrix[0, j - 1] + DNA_dict['-'][seq2[j - 1]]

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            match = score_matrix[i - 1, j - 1] + DNA_dict[seq1[i - 1]][seq2[j - 1]]
            delete = score_matrix[i - 1, j] + DNA_dict[seq1[i - 1]]['-']
            insert = score_matrix[i, j - 1] + DNA_dict['-'][seq2[j - 1]]
            score_matrix[i, j] = max(match, delete, insert)
            if score_matrix[i, j] == match:
                traceback_matrix[i, j] = 1
            elif score_matrix[i, j] == delete:
                traceback_matrix[i, j] = 2
            else:
                traceback_matrix[i, j] = 3

    align1, align2 = "", ""
    i, j = len1, len2
    while i > 0 or j > 0:
        if traceback_matrix[i, j] == 1:
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif traceback_matrix[i, j] == 2:
            align1 = seq1[i - 1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    alignment_score = score_matrix[len1, len2]

    return align1, align2, alignment_score

def Profile(sequences):
    alignment_length = len(sequences[0])
    profile_matrix = np.zeros((5, alignment_length), dtype=int)  # 5 rows: A, C, G, T, -

    nucleotide_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-': 4}

    for seq in sequences:
        for i, nucleotide in enumerate(seq):
            if nucleotide in nucleotide_to_index:
                profile_matrix[nucleotide_to_index[nucleotide], i] += 1

    return profile_matrix

def ProfileAlign(sequence, profile_matrix, DNA_dict):
    profile_length = profile_matrix.shape[1]
    seq_length = len(sequence)
    score_matrix = np.zeros((seq_length + 1, profile_length + 1))
    traceback_matrix = np.zeros((seq_length + 1, profile_length + 1), dtype=int)

    for i in range(1, seq_length + 1):
        score_matrix[i, 0] = score_matrix[i - 1, 0] + DNA_dict['-'][sequence[i - 1]]

    for i in range(1, seq_length + 1):
        for j in range(1, profile_length + 1):
            column = profile_matrix[:, j - 1]
            max_base_index = column[:4].argmax()  
            max_base = 'ACGT'[max_base_index]  # Convert index to base
            match = score_matrix[i - 1, j - 1] + DNA_dict[sequence[i - 1]][max_base]
            delete = score_matrix[i - 1, j] + DNA_dict['-'][sequence[i - 1]]
            insert = score_matrix[i, j - 1] + DNA_dict['-'][max_base]
            score_matrix[i, j] = max(match, delete, insert)
            if score_matrix[i, j] == match:
                traceback_matrix[i, j] = 1
            elif score_matrix[i, j] == delete:
                traceback_matrix[i, j] = 2
            else:
                traceback_matrix[i, j] = 3

    aligned_sequence = ""
    i, j = seq_length, profile_length
    while i > 0 and j > 0:
        if traceback_matrix[i, j] == 1:
            aligned_sequence = sequence[i - 1] + aligned_sequence
            i -= 1
            j -= 1
        elif traceback_matrix[i, j] == 2:
            aligned_sequence = sequence[i - 1] + aligned_sequence
            i -= 1
        else:
            aligned_sequence = '-' + aligned_sequence
            j -= 1

    alignment_score = score_matrix[seq_length, profile_length]

    return aligned_sequence, alignment_score

def ProfileMultipleAlignment(seqs, similarityMatrix=DNA_dict):
    n = len(seqs)
    aligned_seq1, aligned_seq2, _ = SequenceAlign(seqs[0], seqs[1], similarityMatrix)
    MSA = [aligned_seq1, aligned_seq2]

    for i in range(2, n):
        profA = Profile(MSA)
        seq_to_add = seqs[i]
        
        # Align seq_to_add to the profile profA
        aligned_seq_to_add, _ = ProfileAlign(seq_to_add, profA, similarityMatrix)
        
        # Ensure all sequences are the same length
        max_len = max(len(seq) for seq in MSA + [aligned_seq_to_add])
        for j in range(len(MSA)):
            if len(MSA[j]) < max_len:
                MSA[j] += '-' * (max_len - len(MSA[j]))
        if len(aligned_seq_to_add) < max_len:
            aligned_seq_to_add += '-' * (max_len - len(aligned_seq_to_add))
        
        MSA.append(aligned_seq_to_add)

    return MSA

def Score(seq1, seq2, DNA_dict):
    score = 0
    for char1, char2 in zip(seq1, seq2):
        if char1 == '-' or char2 == '-':
            score += DNA_dict['-'][char1] + DNA_dict['-'][char2]
        else:
            score += DNA_dict[char1][char2]

    return score

def SimMatrix(sequences, DNA_dict):
    num_sequences = len(sequences)
    matrix = np.zeros((num_sequences, num_sequences))
    
    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            score = Score(sequences[i], sequences[j], DNA_dict)
            matrix[i][j] = score
            matrix[j][i] = score
    
    return matrix
