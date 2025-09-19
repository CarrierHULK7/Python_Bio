# =====================================================
# Master File: Python Programming for Biology (Exam Helper)
# Covers Rosalind-style problems (1–49 + important extras)
# Each block: short 2-line explanation + minimal core code
# =====================================================

# ------------------ DNA String Basics ------------------

# Hamming Distance
# Counts differing symbols between two equal-length strings
def hamming(s, t):
    return sum(1 for a, b in zip(s, t) if a != b)


# Reverse Complement
# Returns reverse complement of DNA string
def reverse_complement(s):
    comp = str.maketrans("ATCG", "TAGC")
    return s.translate(comp)[::-1]


# GC Content
# Computes percentage of G and C in DNA string
def gc_content(s):
    return (s.count("G") + s.count("C")) / len(s) * 100


# Transcription
# Converts DNA string to RNA string
def transcribe(s):
    return s.replace("T", "U")


# Translation
# Converts RNA string into protein using codon table
def translate(rna):
    table = {
        'AUG':'M','UUU':'F','UUC':'F','UUA':'L','UUG':'L','UCU':'S','UCC':'S','UCA':'S','UCG':'S',
        'UAU':'Y','UAC':'Y','UGU':'C','UGC':'C','UGG':'W','CUU':'L','CUC':'L','CUA':'L','CUG':'L',
        'CCU':'P','CCC':'P','CCA':'P','CCG':'P','CAU':'H','CAC':'H','CAA':'Q','CAG':'Q','CGU':'R',
        'CGC':'R','CGA':'R','CGG':'R','AUU':'I','AUC':'I','AUA':'I','ACU':'T','ACC':'T','ACA':'T',
        'ACG':'T','AAU':'N','AAC':'N','AAA':'K','AAG':'K','AGU':'S','AGC':'S','AGA':'R','AGG':'R',
        'GUU':'V','GUC':'V','GUA':'V','GUG':'V','GCU':'A','GCC':'A','GCA':'A','GCG':'A','GAU':'D',
        'GAC':'D','GAA':'E','GAG':'E','GGU':'G','GGC':'G','GGA':'G','GGG':'G','UAA':'*','UAG':'*','UGA':'*'
    }
    protein = []
    for i in range(0, len(rna)-2, 3):
        codon = rna[i:i+3]
        aa = table.get(codon, '')
        if aa == '*': break
        protein.append(aa)
    return ''.join(protein)


# ------------------ Probability & Counting ------------------

# Factorial
# Computes factorial of n
def fact(n):
    res = 1
    for i in range(2, n+1): res *= i
    return res


# nCr (binomial coefficient)
# Computes n choose r
def nCr(n, r):
    if r < 0 or r > n: return 0
    return fact(n) // (fact(r) * fact(n-r))


# Catalan Numbers
# Counts non-crossing perfect matchings (RNA secondary structure)
def catalan(n):
    return fact(2*n) // (fact(n+1) * fact(n))


# Motzkin Numbers
# Counts non-crossing matchings allowing unmatched nodes
def motzkin(n):
    M = [1,1] + [0]*n
    for i in range(2,n+1):
        M[i] = M[i-1] + sum(M[k]*M[i-2-k] for k in range(i-1))
    return M[n]


# Subset Counting mod 1e6
# Counts subsets of {1..n} mod 1,000,000
def count_subsets(n):
    return pow(2,n,1000000)


# ------------------ Dynamic Programming ------------------

# Longest Common Subsequence
# Finds one LCS between two strings
def lcs(s, t):
    n, m = len(s), len(t)
    dp = [[""]*(m+1) for _ in range(n+1)]
    for i in range(n):
        for j in range(m):
            if s[i] == t[j]:
                dp[i+1][j+1] = dp[i][j] + s[i]
            else:
                dp[i+1][j+1] = max(dp[i][j+1], dp[i+1][j], key=len)
    return dp[n][m]


# Longest Increasing Subsequence
# Finds LIS in list of numbers
def lis(seq):
    from bisect import bisect_left
    sub = []
    for x in seq:
        i = bisect_left(sub, x)
        if i == len(sub): sub.append(x)
        else: sub[i] = x
    return len(sub)


# Edit Distance
# Computes minimum edits (insert, delete, substitute) between two strings
def edit_distance(s, t):
    n, m = len(s), len(t)
    dp = [[0]*(m+1) for _ in range(n+1)]
    for i in range(n+1): dp[i][0] = i
    for j in range(m+1): dp[0][j] = j
    for i in range(1, n+1):
        for j in range(1, m+1):
            cost = 0 if s[i-1] == t[j-1] else 1
            dp[i][j] = min(dp[i-1][j]+1, dp[i][j-1]+1, dp[i-1][j-1]+cost)
    return dp[n][m]


# ------------------ Graph & Trees ------------------

# Tree Edge Count
# For n nodes and k edges, missing edges = components - 1
def edges_to_tree(n, edges):
    from collections import defaultdict, deque
    g = defaultdict(list)
    for u,v in edges:
        g[u].append(v); g[v].append(u)
    seen, comps = set(), 0
    for i in range(1,n+1):
        if i not in seen:
            comps += 1
            dq = deque([i]); seen.add(i)
            while dq:
                u = dq.popleft()
                for v in g[u]:
                    if v not in seen:
                        seen.add(v); dq.append(v)
    return comps-1


# Newick Distance (simple version)
# Returns distance between two labels in a rooted tree
def newick_distance(tree, a, b):
    import re
    # parse (very simplified: assumes binary tree)
    tree = re.sub(r"[();]", "", tree)
    nodes = tree.split(",")
    return abs(nodes.index(a) - nodes.index(b))


# ------------------ Utilities ------------------

# FASTA Parser
# Returns dict of {header:sequence}
def parse_fasta(text):
    seqs, label = {}, None
    for line in text.strip().splitlines():
        if line.startswith(">"):
            label = line[1:]; seqs[label] = ""
        else:
            seqs[label] += line.strip()
    return seqs

