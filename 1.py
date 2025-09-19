def filter_sequences(sequences, quality_scores):
    return [seq for seq, scores in zip(sequences, quality_scores)
            if sum(scores)/len(scores) >= 20 and 50 <= len(seq) <= 200 and seq.count("N")/len(seq) <= 0.1]

n = int(input("Number of sequences: "))
sequences = [input(f"Sequence {i+1}: ") for i in range(n)]
quality_scores = [list(map(int, input(f"Quality scores {i+1}: ").split())) for i in range(n)]

print("Filtered sequences:", filter_sequences(sequences, quality_scores))
