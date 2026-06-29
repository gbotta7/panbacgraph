import os

WORKDIR = "/hlilab/gianfranco/panbacgraph"
configfile: os.path.join(WORKDIR, "workflow/config/config.yaml")

with open(config["sample_list"]) as f:
    sample_list = [l.strip().split("\t")[0] for l in f if l.strip()]

# rule download_dataset:
#     output:
#         asm_zip = temp(os.path.join(WORKDIR, "data/panbac_asm.zip")),
#         genomes_fa = expand(
#             os.path.join(WORKDIR, "data/ncbi_dataset/data/{sample}/{sample}_genomic.fna.gz"),
#             sample=sample_list),
#         accession_list = temp(os.path.join(WORKDIR, "workflow/config/accession_list.txt"))
#     params:
#         genomes_dir = os.path.join(WORKDIR, "data"),
#         accessions = sample_list
#     log:
#         "logs/download_NCBI.log"
#     conda:
#         "../envs/standalone/ncbi.yaml"
#     shell:
#         """
#         rm -rf {params.genomes_dir};
#         mkdir -p {params.genomes_dir};
#         printf '%s\n' {params.accessions} > {output.accession_list};
#         datasets download genome accession \
#             --inputfile {output.accession_list} \
#             --dehydrated \
#             --annotated \
#             --assembly-level complete \
#             --include genome,gff3 \
#             --filename {output.asm_zip} &> {log}
#         unzip {output.asm_zip} -d {params.genomes_dir} &> {log};
#         datasets rehydrate --gzip --directory {params.genomes_dir} &> {log};
#         for acc_dir in {params.genomes_dir}/ncbi_dataset/data/*/; do
#             acc=$(basename "$acc_dir")
#             fna=$(find "$acc_dir" -name "${{acc}}_*_genomic.fna.gz" | head -n1)
#             if [ -n "$fna" ]; then
#                 mv "$fna" "$acc_dir/${{acc}}_genomic.fna.gz"
#             fi
#         done
#         """

        
# rule pilercr_detect_crispr_arrays:
#     input:
#         genome_fa = os.path.join(WORKDIR, "data/ncbi_dataset/data/{sample}/{sample}_genomic.fna.gz")
#     output:
#         crispr_report = os.path.join(WORKDIR, "results/crispr/{sample}_genomic_report.txt")
#     log:
#         "logs/pilercr_crispr_array/{sample}.log"
#     conda:
#         "../envs/standalone/piler.yaml"
#     shell:
#         """
#         tmpfile=$(mktemp /tmp/pilercr_XXXXXX.fna);
#         zcat {input.genome_fa} > "$tmpfile";
#         pilercr -in "$tmpfile" -out {output.crispr_report} 2> {log};
#         rm -f "$tmpfile"
#         """

# rule pilercr_process_report:
#     input:
#         crispr_report = os.path.join(WORKDIR, "results/crispr/{sample}_genomic_report.txt"),
#         genome_fa = os.path.join(WORKDIR, "data/ncbi_dataset/data/{sample}/{sample}_genomic.fna.gz")
#     output:
#         arrays = os.path.join(WORKDIR, "results/spacers/{sample}_arrays.bed"),
#         spacers = os.path.join(WORKDIR, "results/spacers/{sample}_spacers.txt"),
#         mask_genome_fa = os.path.join(WORKDIR, "data/ncbi_dataset/data/{sample}/{sample}_genomic_masked.fna.gz")
#     log:
#         "logs/pilercr_process_report/{sample}.log"
#     conda:
#         "../envs/py/biopython.yaml"
#     script:
#         "../scripts/py/parse_pilercr_report.py"


rule ropebwt3_build_idx:
    input:
        mask_genome_fa = expand(
            os.path.join(WORKDIR, "data/ncbi_dataset/data/{sample}/{sample}_genomic_masked.fna.gz"),
            sample = sample_list)
    output:
        bwt_fmr = os.path.join(WORKDIR, "results/idx/panbac.fmr"),
        bwt_fmd = os.path.join(WORKDIR, "results/idx/panbac.fmd"),
        bwt_ssa = os.path.join(WORKDIR, "results/idx/panbac.fmd.ssa"),
        bwt_len = os.path.join(WORKDIR, "results/idx/panbac.fmd.len.gz")
    log:
        "logs/ropebwt3_idx.log"
    threads:
        24
    conda:
        "../envs/standalone/ropebwt3.yaml"
    shell:
        """
        printf '%s\n' {input.mask_genome_fa} | xargs zcat | \
        ropebwt3 build -t{threads} -bo {output.bwt_fmr} - 2> {log};
        ropebwt3 build -i {output.bwt_fmr} -do {output.bwt_fmd} 2>> {log};
        ropebwt3 ssa -o {output.bwt_ssa} -s8 -t{threads} {output.bwt_fmd} 2>> {log};
        printf '%s\n' {input.mask_genome_fa} | xargs zcat | seqtk comp | cut -f1,2 | gzip > {output.bwt_len} 2>> {log}
        """

rule ropebwt3_search:
    input:
        bwt_fmd = os.path.join(WORKDIR, "results/idx/panbac.fmd"),
        spacers = os.path.join(WORKDIR, "results/spacers/{sample}_spacers.txt")
    output:
        interactions = os.path.join(WORKDIR, "results/interactions/{sample}_spacers_hits.paf")
    log:
        "logs/ropebwt3_search/{sample}.log"
    threads:
        2
    conda:
        "../envs/standalone/ropebwt3.yaml"
    shell:
        """
        if [ -s {input.spacers} ]; then
            ropebwt3 sw -e -N100 -p100 -L -b -t{threads} -m0 -y1 \
                {input.bwt_fmd} {input.spacers} \
                > {output.interactions} 2> {log}
        else
            touch {output.interactions}
            echo "No spacers found, skipping ropebwt3." > {log}
        fi
        """