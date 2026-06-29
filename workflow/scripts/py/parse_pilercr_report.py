import re
import os
import sys
import gzip

from Bio import SeqIO

crispr_report = snakemake.input["crispr_report"]
in_fa_file = snakemake.input["genome_fa"]
out_fa_file = snakemake.output["mask_genome_fa"]
spacers_file = snakemake.output["spacers"]
arrays_file = snakemake.output["arrays"]
log = snakemake.log[0]

with open(log, 'w') as log_f:
    sys.stderr = log_f

    os.makedirs(os.path.dirname(spacers_file), exist_ok=True)

    try:
        with open(crispr_report, 'r') as f:
            content = f.read()

        # parse spacers
        spacers = re.findall(
            r'\d+\s+\d+\s+[\d.]+\s+\d+\s+\S+\s+[A-Z.]+\s+([A-Z]{10,})',
            content
        )
        with open(spacers_file, 'w') as out:
                for s in spacers:
                    out.write(s + '\n')

        # parse arrays' coordinates
        arrays = []
        in_pos_section = False
        for line in content.splitlines():
            if "SUMMARY BY POSITION" in line:
                in_pos_section = True
                continue
            if not in_pos_section:
                continue
            m = re.match(
                r'^\s*(\d+)\s+(\S+.*?)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)',
                line
            )
            if m:
                seq_id = m.group(2).strip()
                start  = int(m.group(3)) # 1-based
                length = int(m.group(4))
                end = start - 1 + length
                # BED: 0-based start, half-open end
                arrays.append((seq_id, start - 1, end))

        with open(arrays_file, 'w') as bed:
            for seq_id, start, end in arrays:
                bed.write(f"{seq_id}\t{start}\t{end}\n")

        if spacers:
            print(f"Written {len(spacers)} spacers to {spacers_file}", file=sys.stderr)
            print(f"Written {len(arrays)} array regions to {arrays_file}", file=sys.stderr)
        else:
            print(f"No spacers found in {crispr_report}, printed empty spacers file", file=sys.stderr)
            print(f"No CRISPR arrays found in {crispr_report}, printed empty BED file", file=sys.stderr)
        
        # remove arrays from fasta file
        with gzip.open(in_fa_file, "rt") as in_ff:
            record = next(SeqIO.parse(in_ff, "fasta"))
            
        for seq_id, bed_start, bed_end in sorted(arrays, key=lambda x: x[1], reverse=True):
            record.seq = record.seq[:bed_start] + 'N' * (bed_end - bed_start) + record.seq[bed_end:]
        
        with gzip.open(out_fa_file, "wt") as out_ff:
            SeqIO.write(record, out_ff, "fasta")

        print(f"Written the genome without CRISPR arrays to file", file=sys.stderr)
            
    except Exception as e:
        print(f"ERROR processing {crispr_report}: {e}", file=sys.stderr)
        raise