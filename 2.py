seq = input()

counts = {nuc: seq.count(nuc) for nuc in "ATGC"}
rev_comp = "".join({"A":"T","T":"A","G":"C","C":"G"}[b] for b in reversed(seq))
gc_content = (seq.count("G") + seq.count("C")) / len(seq) * 100

print("Nucleotide counts:", counts)
print("Reverse complement:", rev_comp)
print("GC content: {:.2f}%".format(gc_content))
