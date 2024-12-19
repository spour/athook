
import regex as re
from typing import List, Tuple, Dict
from at_hook import ATHookManager
from itertools import product
from config import REQUIRED_KEYWORDS, GEL, OUTPUT_DIR, DSSP_EXECUTABLE, FASTA, MOTIFS, PATHS, SUFFIXES, ALPHABET, KMERSIZE, COLLAPSE_RC

def get_kmers(k: int, alphabet: str) -> List[str]:
    return ["".join(kmer) for kmer in product(alphabet, repeat=k)]


if __name__ == "__main__":
    manager = ATHookManager()
    kmers = get_kmers(KMERSIZE, ALPHABET)
    manager.process_fasta(FASTA, n=None, motifs=MOTIFS)
    manager.calculate_all_mean_accessibility(upstream=13, downstream=13)
    for hook_name in manager.hooks.keys():
        manager.process_files(hook_name, paths=PATHS, suffixes=SUFFIXES)
        manager.calculate_e_scores(hook_name, kmers, motifs=MOTIFS, collapse_rc=COLLAPSE_RC)
        
    manager.save_to_json("athook_manager_filtmotif.json")
    df = manager.to_dataframe()
    print("df shape:", df.shape)
