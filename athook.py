
import pandas as pd
from Bio import SeqIO
from typing import Dict, List, Tuple
from tqdm import tqdm
import os
import re
import requests
import time
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
import regex as re
import json
import numpy as np
from glob import glob

from config import REQUIRED_KEYWORDS, GEL, OUTPUT_DIR, DSSP_EXECUTABLE, MOTIFS, PBM_DF

def reverse_complement(seq: str) -> str:
    """reverse complement of a DNA sequence."""
    complement = str.maketrans('ATGC', 'TACG')
    return seq.translate(complement)[::-1]


def find_motif_matches(sequence, motifs=['GRP', 'RGR']):
    """
    find motif matches in a given sequence, return a list of tuples (start, end) for each match.
    """
    matches = []
    pattern = re.compile('|'.join(motifs)) 
    for match in pattern.finditer(sequence):
        matches.append((match.start(), match.end()))
    return matches

def find_grp_matches_with_flanking(sequence, motif='GRP', flank_size=4, required_flank_chars={'R', 'K'}):
    matches = []
    pattern = re.compile(motif)
    
    for match in pattern.finditer(sequence):
        start, end = match.start(), match.end()

        upstream = sequence[max(0, start-flank_size):start]
        downstream = sequence[end:end+flank_size]
        
        if len(upstream) == flank_size and len(downstream) == flank_size:
            if (sum([1 for char in upstream if char in required_flank_chars]) >= 1 and
                sum([1 for char in downstream if char in required_flank_chars]) >= 1):
                matches.append((start-flank_size, end+flank_size))
    
    return matches



class ATHook:
    """
     AT-hook motif within a protein sequence.
    """
    def __init__(self, accession: str, name: str, sequence: str, matches: List[Tuple[int, int]] = None, flag: int = 0):
        self.accession = accession
        self.name = name
        self.sequence = sequence
        self.matches = matches or []
        self.flag = flag  # if the protein is from GEL list
        self.mean_accessibility = []
        self.secondary_structure = []
        self.amino_acid_sequence = []
        self.failed = False
        self.e_scores = {"ME": {}, "HK": {}}
        self.pth_number = None 
        self.associated_files = {}
        self.pth_sequence = None
        self.pbm_df = pd.read_csv(PBM_DF)
        self.pbm_df.dropna(subset=["Protein Sequence"], inplace=True)
    
    def add_file(self, key: str, file_path: str):
        """
        ME or HK file 
        """
        self.associated_files[key] = file_path
        
    def set_pth_number(self, pth_number: str):
        """
        set the pTH number.
        """
        self.pth_number = pth_number
        assert self.pth_number.startswith("pTH"), f"Invalid pTH number: {self.pth_number}"
        if self.pth_number in self.pbm_df["pTH No"].values:
            self.pth_sequence = self.pbm_df[self.pbm_df["pTH No"] == self.pth_number]["Protein Sequence"].values[0]
        else:
            print(f"pTH number {self.pth_number} not found in PBM data.")

    def add_e_scores(self, kmer: str, score: float, me_or_hk: str):
        self.e_scores[me_or_hk][kmer] = score
        

    def fetch_alphafold_structure(self) -> str:
        """
        get  AlphaFold PDB file for the protein from the EBI server.
        """
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        pdb_file_path = f"{OUTPUT_DIR}{self.accession}_alphafold.pdb"

        if os.path.exists(pdb_file_path):
            print(f"Using cached file for {self.accession}")
            return pdb_file_path

        base_url = f"https://alphafold.ebi.ac.uk/files/AF-{self.accession}-F1-model_v4.pdb"

        session = requests.Session()
        retries = Retry(
            total=5,
            backoff_factor=1,
            status_forcelist=[500, 502, 503, 504],
            raise_on_status=False,
        )
        session.mount("https://", HTTPAdapter(max_retries=retries))

        try:
            response = session.get(base_url, timeout=10)
            response.raise_for_status()
            with open(pdb_file_path, "w") as pdb_file:
                pdb_file.write(response.text)
            time.sleep(1)  
            return pdb_file_path
        except requests.exceptions.RequestException as e:
            print(f"Failed to fetch {self.accession}: {e}")
            self.failed = True
            return None

    def parse_pdb_for_residue_data(self, pdb_file_path: str):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file_path)
        model = structure[0]
        dssp = DSSP(model, pdb_file_path, dssp=DSSP_EXECUTABLE)

        residue_data = {
            'residue': [],
            'amino_acid': [],
            'secondary_structure': [],
            'accessibility': []
        }

        for key in dssp.keys():
            residue_data['residue'].append(key[1][1])
            residue_data['amino_acid'].append(dssp[key][1])
            residue_data['secondary_structure'].append(dssp[key][2])
            residue_data['accessibility'].append(dssp[key][3])

        return residue_data

    def calculate_mean_accessibility(self, upstream: int = 5, downstream: int = 5):

        pdb_file_path = self.fetch_alphafold_structure()
        if pdb_file_path is None or self.failed:
            return

        try:
            residue_data = self.parse_pdb_for_residue_data(pdb_file_path)
        except Exception as e:
            print(f"Error parsing {self.accession}: {e}")
            self.failed = True
            return

        for start, end in self.matches:
            start_upstream = max(0, start - upstream)
            end_downstream = end + downstream
            accessibility_slice = residue_data['accessibility'][start_upstream:end_downstream]
            sec_struct_slice = residue_data['secondary_structure'][start_upstream:end_downstream]
            aa_seq_slice = residue_data['amino_acid'][start_upstream:end_downstream]

            if len(accessibility_slice) == 0:
                continue

            mean_acc = sum(accessibility_slice) / len(accessibility_slice)
            self.mean_accessibility.append(mean_acc)
            self.secondary_structure.append(''.join(sec_struct_slice))
            self.amino_acid_sequence.append(''.join(aa_seq_slice))



class ATHookManager:
    """
    Handle multiple ATHook objects.
    """
    def __init__(self):
        self.hooks: Dict[str, ATHook] = {}

    def process_fasta(self, file_path: str, n = None, motifs=MOTIFS):
        i = 0
        for record in SeqIO.parse(file_path, "fasta"):
            i += 1
            if n and i > n:
                break
            full_header = record.description
            header = record.id
            sequence = str(record.seq)

            try:
                name = full_header.split("GN=")[1].split(" ")[0]
            except IndexError:
                continue
            if "Fragment" in full_header:
                continue

            matches = find_motif_matches(sequence, motifs)
            flag = 0

            if any(keyword in full_header for keyword in REQUIRED_KEYWORDS):
                assert len(matches) > 0, f"No GRP matches found for positive control {full_header}"

            keyword_in_gel = next((keyword for keyword in GEL if keyword in full_header), None)
            if keyword_in_gel:
                flag = 1
                if not matches:
                    for tr_record in SeqIO.parse(file_path, "fasta"):
                        tr_full_header = tr_record.description
                        tr_sequence = str(tr_record.seq)
                        if keyword_in_gel in tr_full_header:
                            matches = find_motif_matches(tr_sequence, motifs)
                            if matches:
                                sequence = tr_sequence
                                header = tr_record.id
                                break

            if matches or flag == 1:
                accession = header.split("|")[1]
                self.hooks[header] = ATHook(accession, name, sequence, matches, flag)

    def calculate_all_mean_accessibility(self, upstream: int = 5, downstream: int = 5):
        for hook in tqdm(self.hooks.values(), desc="Calculating mean accessibility"):
            hook.calculate_mean_accessibility(upstream, downstream)

    def to_dataframe(self) -> pd.DataFrame:
        data = []
        for hook in self.hooks.values():
            if hook.failed or not hook.mean_accessibility:
                continue
            for mean_acc, sec_struc, aa_seq in zip(hook.mean_accessibility, hook.secondary_structure, hook.amino_acid_sequence):
                data.append({
                    'mean_acc': mean_acc,
                    'sec_struc': sec_struc,
                    'amino_acid': aa_seq,
                    'accessions': hook.accession,
                    'names': hook.name,
                    'flags': hook.flag
                })
        return pd.DataFrame(data)
    
    def process_files(self, name: str, paths: List[str], suffixes: List[str]):
        name_protein = self.hooks[name].name
        patterns = [os.path.join(path, f"*{name_protein}{suffix}") for path in paths for suffix in suffixes]
        files = [file for pattern in patterns for file in glob(pattern) if "zscores" not in file]
        if not files:
            print(f"No files found for {name}")
            return

        for file in files:
            if "1M-ME" in file:
                self.hooks[name].add_file("ME", file)
            elif "1M-HK" in file:
                self.hooks[name].add_file("HK", file)
            if "pTH" in file:
                pth_number = file.split("pTH")[1].split("_")[0].split(".")[0]
                self.hooks[name].set_pth_number(f"pTH{pth_number}")

    def calculate_e_scores(
        self, name: str, kmers: List[str], motifs: List[str], collapse_rc: bool = True
    ):
        # breakpoint()
        if name not in self.hooks:
            print(f"Error: AT-hook '{name}' not found in hooks.")
            return

        hook = self.hooks[name]

        if "ME" not in hook.associated_files or "HK" not in hook.associated_files:
            print(f"Skipping {hook.name}: Missing associated ME or HK files.")
            return

        if not hasattr(hook, "pth_sequence") or hook.pth_sequence is None or not any(motif in hook.pth_sequence for motif in motifs):
            print(f"Skipping {hook.name}: Does not contain any specified motifs.")
            return


        for kmer in tqdm(kmers, desc=f"Calculating e-scores for {hook.name}"):
            try:
                e_score_me = self.calculate_e_score(hook.associated_files["ME"], kmer, collapse_rc)
                e_score_hk = self.calculate_e_score(hook.associated_files["HK"], kmer, collapse_rc)
                hook.add_e_scores(kmer, e_score_me, "ME")
                hook.add_e_scores(kmer, e_score_hk, "HK")
            except Exception as e:
                print(f"Error calculating e-scores for {kmer} in {hook.name}: {e}")

    
    @staticmethod
    def calculate_e_score(file_path: str, candidate_binding_site: str, collapse_rc: bool,  signal_column: str = 'mean_signal_intensity') -> float:
        data1 = pd.read_csv(file_path, sep="\t")

        if "mean_signal_intensity" not in data1.columns or "pbm_sequence" not in data1.columns:
            raise ValueError(f"File {file_path} must contain 'mean_signal_intensity' and 'pbm_sequence' columns.")

        data = data1.copy()
        
        data = data[data['flag'] != 1]
        
        if signal_column not in data.columns or 'pbm_sequence' not in data.columns:
            raise ValueError(f"DataFrame must contain '{signal_column}' and 'pbm_sequence' columns.")
        
        if collapse_rc:
            candidate_rc = reverse_complement(candidate_binding_site)
            data['type'] = np.where(
                data['pbm_sequence'].str.contains(candidate_binding_site) | data['pbm_sequence'].str.contains(candidate_rc), 'f', 'b')
        else:
            data['type'] = np.where(data['pbm_sequence'].str.contains(candidate_binding_site), 'f', 'b')
        
        data = data.sort_values(signal_column, ascending=False).reset_index(drop=True)
        foreground = data[data['type'] == 'f']
        background = data[data['type'] == 'b']
        if foreground.empty or background.empty:
            return 0

        top_half_foreground = foreground.nlargest(int(len(foreground) * 0.5), signal_column)
        top_half_background = background.nlargest(int(len(background) * 0.5), signal_column)
        sum_fg = top_half_foreground.index.values.sum()
        sum_bg = top_half_background.index.values.sum()
        len_fg = len(top_half_foreground)
        len_bg = len(top_half_background)
        
        if len_fg == 0 or len_bg == 0:
            raise ValueError("Top half of either foreground or background has no data. Cannot calculate enrichment score.")
        
        size_corr = 1 / (len_fg + len_bg)
        wcw_fg = sum_fg / len_fg
        wcw_bg = sum_bg / len_bg
        wcw = wcw_bg - wcw_fg
        ans = size_corr * wcw
        
        return ans
    
    def to_dict(self) -> Dict:
        return {
            name: {
                "accession": hook.accession,
                "name": hook.name,
                "sequence": hook.sequence,
                "matches": hook.matches,
                "flag": hook.flag,
                "mean_accessibility": hook.mean_accessibility,
                "secondary_structure": hook.secondary_structure,
                "amino_acid_sequence": hook.amino_acid_sequence,
                "failed": hook.failed,
                "e_scores": hook.e_scores,
                "pth_number": hook.pth_number,
                "associated_files": hook.associated_files,
                "pth_sequence": hook.pth_sequence,
            }
            for name, hook in self.hooks.items()
        }

    def save_to_json(self, file_path: str):
        with open(file_path, 'w') as f:
            json.dump(self.to_dict(), f, indent=4)
        print(f"ATHookManager saved to {file_path}")

    @staticmethod
    def load_from_json(file_path: str):
        with open(file_path, 'r') as f:
            data = json.load(f)

        manager = ATHookManager()
        for name, hook_data in data.items():
            hook = ATHook(
                accession=hook_data["accession"],
                name=hook_data["name"],
                sequence=hook_data["sequence"],
                matches=hook_data["matches"],
                flag=hook_data["flag"],
            )
            hook.mean_accessibility = hook_data["mean_accessibility"]
            hook.secondary_structure = hook_data["secondary_structure"]
            hook.amino_acid_sequence = hook_data["amino_acid_sequence"]
            hook.failed = hook_data["failed"]
            hook.e_scores = hook_data["e_scores"]
            hook.pth_number = hook_data["pth_number"]
            hook.associated_files = hook_data["associated_files"]
            hook.pth_sequence = hook_data["pth_sequence"]
            manager.hooks[name] = hook

        print(f"ATHookManager loaded from {file_path}")
        return manager
