import regex as re
from typing import List, Tuple, Dict
from at_hook import ATHookManager
from config import REQUIRED_KEYWORDS, GEL, OUTPUT_DIR, DSSP_EXECUTABLE, FASTA, MOTIFS, PATHS, SUFFIXES, ALPHABET, KMERSIZE, COLLAPSE_RC

filtered_manager = ATHookManager().load_from_json("athook_manager_filtered.json")

data = []
for name, hook in filtered_manager.hooks.items():
    if len(hook.e_scores["ME"]) > 0 or len(hook.e_scores["HK"]) > 0:
        for kmer, score in hook.e_scores["ME"].items():
            data.append({"kmer": kmer, "score": score, "type": "ME", "hook_name": hook.name})
        for kmer, score in hook.e_scores["HK"].items():
            data.append({"kmer": kmer, "score": score, "type": "HK", "hook_name": name})

import pandas as pd
df = pd.DataFrame(data)
df_kmers = df.pivot_table(index="kmer", columns=["hook_name", "type"], values="score").reset_index()

df_kmers.columns = ['_'.join(filter(None, col)).strip() for col in df_kmers.columns.values]

def reverse_complement(seq: str) -> str:
    complement = str.maketrans('ACGT', 'TGCA')
    return seq.translate(complement)[::-1]

def drop_reverse_complements(df, kmer_col="kmer"):
    seen = set()
    rows_to_keep = []
    for index, kmer in df[kmer_col].iteritems():
        rev_comp = reverse_complement(kmer)
        if kmer not in seen and rev_comp not in seen:
            rows_to_keep.append(index)
            seen.add(kmer)
            seen.add(rev_comp)
    return df.loc[rows_to_keep]


df_cleaned = drop_reverse_complements(df_kmers, kmer_col="kmer")
