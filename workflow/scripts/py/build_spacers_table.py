import sqlite3
import sys
from pathlib import Path


def parse_spacer_file(path: Path):
    """
    Yield spacer sequences from a plain text file, one spacer per line.
    """
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            yield line


def get_sample_name(spacer_path: Path) -> str:
    """
    Derive the sample/genome name from a spacer filename, e.g.
    'GCF_000005845.2_spacers.txt' -> 'GCF_000005845.2'
    """
    name = spacer_path.name
    suffix = "_spacers.txt"
    if name.endswith(suffix):
        name = name[: -len(suffix)]
    return name


def build_spacer_table(spacer_files, genome_db_path: str, spacer_db_path: str) -> None:
    conn = sqlite3.connect(spacer_db_path)
    cur = conn.cursor()

    cur.execute("""
        CREATE TABLE IF NOT EXISTS spacer_table (
            spacer_id       INTEGER PRIMARY KEY AUTOINCREMENT,
            spacer_sequence TEXT NOT NULL,
            source_fname    TEXT NOT NULL
        )
    """)

    cur.execute("PRAGMA journal_mode = WAL")
    cur.execute("PRAGMA synchronous = OFF")

    rows = []
    for spacer_file in spacer_files:
        spacer_path = Path(spacer_file)
        source_fname = get_sample_name(spacer_path)

        for spacer_seq in parse_spacer_file(spacer_path):
            rows.append((spacer_seq, source_fname))

    cur.executemany(
        """
        INSERT INTO spacer_table (spacer_id, spacer_sequence, source_fname)
        VALUES (?, ?, ?)
        """,
        rows,
    )

    # Useful for deduplication / lookup by sequence later
    cur.execute("CREATE INDEX IF NOT EXISTS idx_spacer_sequence ON spacer_table (spacer_sequence)")

    conn.commit()
    conn.close()


spacer_files = snakemake.input["spacers"]
genome_db_path = snakemake.input["genome_db"]
spacers_db_path = snakemake.output["spacers_db"]
log = snakemake.log[0]

with open(log, "w") as log_f:
    sys.stderr = log_f
    build_spacer_table(spacer_files, genome_db_path, spacers_db_path)