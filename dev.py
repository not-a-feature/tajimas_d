
seqs = [
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    "AAAAAAAAAAAAABAAAAAAAABAAAAAAAAAAAAAAAAABAAAAAAA",
    "AAAAACAAAAAAABAAAAAABABAAAAAAABABABAAAAAAAAAAAAA",
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    "AAAAAAAAAAAAABAAAAAAAABAAAAAAAAAAAAAAAAABAAAAAAA",
    "AAAAACAAAAAAABAAAAAABABAAAAAAABABABAAAAAAAAAAAAA",
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
]

def org(sequences) -> int:
    """
    Counts the number of segregating sites.

    Input:
        sequences: list of str, list of sequences.
    Returns:
        seg_sites: int, number of segregating sites.
    """
    seg_sites = 0
    # For each position in sequence
    for i in range(0, len(sequences[0])):
        s = sequences[0][i]    # Take the char of first sequence as reference
        for seq in sequences:  # For each other sequence
            if not seq[i] == s:
                seg_sites += 1  # Add 1 if not equal
                break  # And stop comparing this position
    return seg_sites

def wany(sequences):
    return sum([any([seq[i] != sequences[0][i] for seq in sequences]) for i in range(0, len(sequences[0]))])



# 8
#t_org  = timeit.timeit(lambda:org(seqs), number=10000)
#t_wany = timeit.timeit(lambda:wany(seqs), number=10000)
#print(f"{t_org=}")
#print(f"{t_wany=}")
