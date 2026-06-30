import gzip
from pathlib import Path
import sqlite3
import sys

def get_genome_name(fasta_path: Path) -> str:
    """Derive a genome accession from the filename, e.g.
    'GCF_000005845.2_ASM584v2_genomic.fna.gz' -> 'GCF_000005845.2'
    """
    name = fasta_path.name
    if name.endswith("_genomic.fna.gz"):
        name = name[: -len("_genomic.fna.gz")]
    parts = name.split("_")
    return "_".join(parts[:2])


def extract_chromosome_hits(fa_files) -> dict:
    """
    Returns a dict: {chromosome_accession: genome_name}
    """
    chrom_dict = {}

    for fa_file in fa_files:
        fa_path = Path(fa_file)
        genome_name = get_genome_name(fa_path)
        
        with gzip.open(fa_path, "rt") as f:
            for line in f:
                if not line.startswith(">"):
                    continue  # skip sequence lines, only headers matter here

                header = line[1:].strip()
                chrom_id = header.split()[0]  # first token = accession
                chrom_dict[chrom_id] = genome_name

    return chrom_dict


def build_chrom_map_db(fa_files, genome_db_path: str, chrom_map_db_path: str) -> None:
    chrom_dict = extract_chromosome_hits(fa_files)

    # Build genome_table
    genome_conn = sqlite3.connect(genome_db_path)
    genome_conn.execute("PRAGMA foreign_keys = ON")
    gcur = genome_conn.cursor()

    gcur.execute("""
        CREATE TABLE IF NOT EXISTS genome_table (
            fname TEXT PRIMARY KEY
        )
    """)
    unique_fnames = set(chrom_dict.values())
    gcur.executemany(
        "INSERT OR IGNORE INTO genome_table (fname) VALUES (?)",
        [(f,) for f in unique_fnames],
    )
    genome_conn.commit()
    genome_conn.close()

    # Build chrom_map_table
    chrom_conn = sqlite3.connect(chrom_map_db_path)
    chrom_conn.execute("PRAGMA foreign_keys = ON")
    ccur = chrom_conn.cursor()

    ccur.execute("""
        CREATE TABLE IF NOT EXISTS chrom_map_table (
            chrom_name TEXT PRIMARY KEY,
            fname TEXT NOT NULL
        )
    """)
    ccur.executemany(
        "INSERT INTO chrom_map_table (chrom_name, fname) VALUES (?, ?)",
        chrom_dict.items(),
    )
    chrom_conn.commit()
    chrom_conn.close()

fa_files = snakemake.input["genome_fa"]
genome_db_path = snakemake.output["genome_db"]
chrommap_db_path = snakemake.output["chrommap_db"]
log = snakemake.log[0]

with open(log, 'w') as log_f:
    sys.stderr = log_f

    build_chrom_map_db(fa_files, genome_db_path, chrommap_db_path)