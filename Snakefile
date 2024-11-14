#!/usr/bin/env python3

from pathlib import Path


def get_fasta(wildcards):
    if wildcards.fasta_file == "all":
        return Path(outdir, "005_combined-fasta-files", "all.fasta")
    else:
        return fasta_dict[wildcards.fasta_file]


# GLOBALS
data_dir = Path("data")
outdir = Path("output")
logdir = Path(outdir, "logs")
kmer_counts = Path(outdir, "020_kmer-counts")
usearch_output = Path(outdir, "030_usearch-clusters")

sequence_identities_to_check = [80, 85, 87, 90, 95, 96, 97, 98, 99]

# containers
bbmap = "docker://quay.io/biocontainers/bbmap:39.01--h92535d8_1"
kmc = "docker://quay.io/biocontainers/kmc:3.2.4--h6dccd9a_2"
r = "docker://ghcr.io/tomharrop/r-containers:r2u_24.04_cv1"
usearch = "docker://quay.io/biocontainers/usearch:12.0_beta--h9ee0642_1"

# MAIN
fasta_files = data_dir.rglob("*.fasta")

fasta_dict = {}
alignment_dict = {}

# which bait sizes do we want to analyse?
bait_sizes = ["80", "120"]

# get the fasta files
for fasta_file in fasta_files:
    my_stem = fasta_file.stem
    if not (my_stem.startswith(("_", "."))):
        if (
            "Angiosperms353_alignments_PerezEscobar2024_Diurideae_Pterostylidinae"
            in fasta_file.parts
        ):
            alignment_dict[my_stem] = fasta_file
        else:
            fasta_dict[my_stem] = fasta_file


# add the combined alignment file to the fasta dict
fasta_dict[
    "Angiosperms353_alignments_PerezEscobar2024_Diurideae_Pterostylidinae"
] = Path(
    outdir,
    "010_alignment-files",
    "Angiosperms353_alignments_PerezEscobar2024_Diurideae_Pterostylidinae.fasta",
)

all_fastas = sorted(set(fasta_dict.keys())) + ["all"]


rule collate_counts:
    input:
        unclustered_kmer_counts=expand(
            Path(
                kmer_counts,
                "{fasta_file}.k{{k}}",
                "kmer_counts.json",
            ),
            fasta_file=all_fastas,
        ),
        clustered_kmer_counts=expand(
            Path(
                usearch_output,
                "{fasta_file}.k{{k}}",
                "clusters.id{sequence_identity}.count.csv",
            ),
            fasta_file=all_fastas,
            sequence_identity=sequence_identities_to_check,
        ),
    output:
        csv=Path(outdir, "040_collated-kmer-counts", "counts.k{k}.csv"),
    log:
        Path(logdir, "collate_counts.k{k}.log"),
    container:
        r
    script:
        "src/collate_counts.R"


rule count_clustered_kmers:
    input:
        Path(
            usearch_output,
            "{fasta_file}.k{k}",
            "clusters.id{sequence_identity}.fasta",
        ),
    output:
        Path(
            usearch_output,
            "{fasta_file}.k{k}",
            "clusters.id{sequence_identity}.count.csv",
        ),
    shell:
        "printf '%s,%s,%s\\n' "
        "{wildcards.fasta_file} "
        "{wildcards.sequence_identity} "
        "$( grep -c '^>' {input} ) "
        "> {output} "


rule usearch:
    input:
        counts=Path(
            kmer_counts,
            "{fasta_file}.k{k}",
            "kmer_counts.txt",
        ),
    output:
        Path(
            usearch_output,
            "{fasta_file}.k{k}",
            "clusters.id{sequence_identity}.fasta",
        ),
    log:
        Path(
            logdir,
            "usearch",
            "{fasta_file}.k{k}.id{sequence_identity}.log",
        ),
    threads: 32
    resources:
        mem_mb=int(6e3),
        time=360,
    container:
        usearch
    shadow:
        "minimal"
    shell:
        'awk \'{{print ">seq" NR "_" $2 "\\n" $1}}\' {input.counts} > kmer_counts.fasta '
        "&& "
        "usearch "
        "-threads {threads} "
        "-cluster_fast kmer_counts.fasta "
        "-id 0.{wildcards.sequence_identity} "
        "-centroids clusters.fasta "
        "&> {log} "
        "&& "
        "mv clusters.fasta {output}"


rule count_kmers:
    input:
        fasta_file=get_fasta,
    output:
        kmc_pre=Path(
            kmer_counts,
            "{fasta_file}.k{k}",
            "kmer_counts.kmc_pre",
        ),
        kmc_suf=Path(
            kmer_counts,
            "{fasta_file}.k{k}",
            "kmer_counts.kmc_suf",
        ),
        counts=Path(
            kmer_counts,
            "{fasta_file}.k{k}",
            "kmer_counts.txt",
        ),
        stats=Path(
            kmer_counts,
            "{fasta_file}.k{k}",
            "kmer_counts.json",
        ),
    params:
        outdir=lambda wildcards, output: Path(output.kmc_pre).parent,
    log:
        Path(
            logdir,
            "count_kmers",
            "{fasta_file}.k{k}.log",
        ),
    threads: 2
    resources:
        mem_mb=int(2e3),
    shadow:
        "minimal"
    container:
        kmc
    shell:
        "kmc "
        "-k{wildcards.k} "
        "-ci1 "
        '"-m$(( {resources.mem_mb}/1000 ))" '
        "-fm "
        "-jkmer_counts.json "
        "-t{threads} "
        "{input.fasta_file} "
        "kmer_counts "
        "$(mktemp -d) "
        "&> {log} "
        "&& "
        "kmc_tools "
        "transform "
        "kmer_counts "
        "dump "
        "kmer_counts.txt "
        "&>> {log} "
        "&& "
        "mv kmer_counts.* {params.outdir}/ "


rule combine_fastas:
    input:
        fasta_dict.values(),
    output:
        Path(outdir, "005_combined-fasta-files", "all.fasta"),
    shell:
        "cat {input} > {output}"


rule combine_alignments:
    input:
        expand(
            Path(outdir, "010_alignment-files", "{alignment}.fasta"),
            alignment=alignment_dict.keys(),
        ),
    output:
        fasta_dict[
            "Angiosperms353_alignments_PerezEscobar2024_Diurideae_Pterostylidinae"
        ],
    shell:
        "cat {input} > {output}"


rule dotdashxton:
    input:
        fasta=lambda wildcards: alignment_dict[wildcards.alignment],
    output:
        fasta=Path(outdir, "010_alignment-files", "{alignment}.fasta"),
    log:
        Path(logdir, "dotdashxton", "{alignment}.log"),
    threads: 1
    resources:
        time=1,
    container:
        bbmap
    shell:
        "reformat.sh "
        "in={input.fasta} "
        "out={output.fasta} "
        "dotdashxton=t "
        "2> {log}"


# target
rule target:
    default_target: True
    input:
        expand(
            Path(outdir, "040_collated-kmer-counts", "counts.k{k}.csv"),
            k=bait_sizes,
        ),
