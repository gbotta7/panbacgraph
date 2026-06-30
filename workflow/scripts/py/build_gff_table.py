import gzip
from glob import glob
from pathlib import Path
import sqlite3
import sys

KEEP_TYPES = {"CDS", "ncRNA", "pseudogene", "mobile_genetic_element", "sequence_feature"}

def parse_gff3(path: Path) -> list[dict]:
    """Yield records from a GFF3 file as dicts."""
    opener = gzip.open(path, "rt") if path.suffix == ".gz" else open(path)
    records = []
    with opener as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = cols
            if ftype not in KEEP_TYPES:
                continue
            phage = any("phage" in tag_val.lower() for tag_val in attrs.split(";"))
            records.append({
                "type": ftype,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "phage": phage,
                "file": str(path.parent.name),
            })
    return records

def build_gff_db(gff3_files, genome_db_path, gff_db_path):
    conn = sqlite3.connect(gff_db_path)
    cur = conn.cursor()

    cur.execute("""
        CREATE TABLE gff_table (
            fname        TEXT NOT NULL,
            feature_type TEXT NOT NULL,
            start_pos    INTEGER NOT NULL,
            end_pos      INTEGER NOT NULL,
            strand       CHAR(1) CHECK (strand IN ('+', '-')),
            phage        BOOLEAN NOT NULL DEFAULT FALSE,
            PRIMARY KEY (fname, start_pos, end_pos, feature_type)
        )
    """)

    cur.execute("PRAGMA journal_mode = WAL")
    cur.execute("PRAGMA synchronous = OFF")

    for path in gff3_files:
        path = Path(path)
        rows = []
        for rec in parse_gff3(path):
            rows.append((
                rec["file"],
                rec["type"],
                rec["start"],
                rec["end"],
                rec["strand"],
                int(rec["phage"]),
            ))
        cur.executemany(
            "INSERT INTO features VALUES (?,?,?,?,?,?)", rows
        )

    # Create B-tree range index
    cur.execute("CREATE INDEX idx_range ON features(file, start, end)")
    conn.commit()
    conn.close()


gff_files = snakemake.input["anno_gff"]
genome_db_path = snakemake.input["genome_db"]
gff_db_path = snakemake.output["gff_db"]
log = snakemake.log[0]

with open(log, 'w') as log_f:
    sys.stderr = log_f

    build_gff_db(gff_files, genome_db_path, gff_db_path)
