class node:
    def __init__(self, kmer):
        self.sequence = kmer
        self.k = len(kmer)
        self.next = []
        self.previous = []


def get_kmer_count_from_sequence(sequence, k=3, cyclic=True):
    """
    Returns dictionary with keys representing all possible kmers in a sequence
    and values counting their occurrence in the sequence.
    """
    # dict to store kmers
    kmers = {}

    # count how many times each occurred in this sequence (treated as cyclic)
    for i in range(0, len(sequence)):
        kmer = sequence[i:i + k]

        # for cyclic sequence get kmers that wrap from end to beginning
        length = len(kmer)
        if cyclic:
            if len(kmer) != k:
                kmer += sequence[:(k - length)]

        # if not cyclic then skip kmers at end of sequence
        else:
            if len(kmer) != k:
                continue

        # count occurrence of this kmer in sequence
        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1

    return kmers

sequences = ['LNNFYPREAKVQWK']
kmers = []
for sequence in sequences:
    kmer = get_kmer_count_from_sequence(sequence,k=3,cyclic=False)

print(kmers)

for kemr in kmers:
    print(kemr)