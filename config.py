
FASTA = "/home/hugheslab1/spour/athook_roundn/uniprotkb_proteome_UP000005640_2024_08_30.fasta"
OUTPUT_DIR = "/home/hugheslab1/spour/athook_roundn/afold/"
DSSP_EXECUTABLE = '/home/hugheslab1/spour/miniforge3/envs/newenv/bin/mkdssp'
PATHS = ["/home/hugheslab1/spour/athook_3_2023/new_athooks", 
         "/home/hugheslab1/spour/AT_Hook_2/all/only_1x_at",
         "/home/hughespub/pbmTest/"
]
SUFFIXES = [".txt", ".AThook.txt", ".DBD.txt"]
MOTIFS = ["GRP", "RGR"]
KMERSIZE = 4
ALPHABET = "ACGT"
COLLAPSE_RC = True
PBM_DF = "/home/hugheslab1/spour/athook_roundn/pbm_sequences.csv"

REQUIRED_KEYWORDS = [
    '|MECP2', 'HMGA1', 'HMGA2', 'Q49A26', 'O00555', '|ELF3', 
    'Q03164', 'CBX2', 'AHDC1', 'FAM171B', 'Q6ZRS2', 'KMT2C', 
    'Q8NFC6', 'Q9UMN6', 'SETBP1', 'PSIP1', "Q13568", "MLL ", 
    "LHX2", "Q96T58", "AKNA_", " RFX5", "O14646", "BAZ2B", 
    "BARX1", "ZBTB24", "CBX2", "Q9NR48", "AHDC1", "Q8TEK3", 
    "Q9BVI0", "PRR12", "AHCTF"
]

GEL = ["ETV5", "IRF5", "CENPX", "FOXO6", "INO80", "SRRM5", "ELF3", "DOT1L", "DNM3A", "MECP2"]

AA = "ACDEFGHIKLMNPQRSTVWY"
