def kmers(sequence, k):
    
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def k_comparison(seq1, seq2, k):
    
    kmers1 = set(kmers(seq1, k))
    kmers2 = set(kmers(seq2, k))
    return len(kmers1.intersection(kmers2))

def MerMatrix(sequences, k):
    
    import numpy as np
    num_sequences = len(sequences)
    matrix = np.zeros((num_sequences, num_sequences))
    
    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            score = k_comparison(sequences[i], sequences[j], k)
            matrix[i][j] = score
            matrix[j][i] = score
    
    return matrix
