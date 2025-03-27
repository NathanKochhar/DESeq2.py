'''
make classes
'''


class DESeq_obj:
    def __init__(self, raw_counts, conds, normalized_counts):
        self.raw_counts = raw_counts
        self.conds = conds
        self.normalized_counts = normalized_counts