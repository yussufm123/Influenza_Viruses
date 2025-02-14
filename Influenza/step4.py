def WordSeq(sequence, word_length):
    
    return [sequence[i:i+word_length] for i in range(len(sequence) - word_length + 1)]

def words_comparison(seq1, seq2, word_length):
    
    words1 = set(WordSeq(seq1, word_length))
    words2 = set(WordSeq(seq2, word_length))
    return len(words1.intersection(words2))

def LZMatrix(sequences, word_length):
    
    import numpy as np
    num_sequences = len(sequences)
    matrix = np.zeros((num_sequences, num_sequences))
    
    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            score = words_comparison(sequences[i], sequences[j], word_length)
            matrix[i][j] = score
            matrix[j][i] = score
    
    return matrix
