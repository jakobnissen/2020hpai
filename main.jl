cd("/Users/jakobnissen/Documents/ssi/projects/2020/hpai")

using Pkg
Pkg.activate(".")

import Conda
using BioSequences
using BioAlignments
using FASTX
using Plots
using CodecZlib

SEGMENTS = ["HA", "NA", "MP", "NP", "PB1", "PB2", "PA", "NS"]
maybedir(path::AbstractString) = isdir(path) || mkdir(path)

## Step 1: Move inputs to raw directory
include("scripts/01_gather_inputs.jl")
println("Completed step 1 - Move input")

## Step 2: QC reads - look at these reports manually.
include("scripts/02_read_qc.jl")
basenames = readdir("results/fastqc")
println("Completed step 2 - QC reads")

## Step 3: Split ref by segment
include("scripts/03_split_ref_by_segment.jl")
println("Completed step 3 - Split references")

## Step 4: KMA align
include("scripts/04_kma.jl")
println("Completed step 4 - Align to references")

## Step 5: Remove primers and cat consensus together
include("scripts/05_consensus.jl")
println("Completed step 5 - Create consensus seq")

## Step 6: Get depths plot
include("scripts/06_depths.jl")
println("Completed step 6 - Plot depths")

## Step 7: Compare consensus with the ones already uploaded
include("scripts/07_compare_consensus.jl")
println("Completed step 7 - Compare consensus seq")

## Step 8: Translate the longest ORF to AA for each segment
# Relative to ref/relative.fna (which is 14535-1 annotated by an online tool)
GENE_ORFS = [
    ("HA", "HA", [17:1720]),
    ("M1", "MP", [14:772]),
    ("M2", "MP", [14:39, 728:995]),
    ("NA", "NA", [9:1421]), # 1:1419 for N5
    ("NP", "NP", [34:1530]),
    ("NS1", "NS", [15:668]),
    ("NEP", "NS", [15:44, 517:852]),
    ("PA", "PA", [13:2163]),
    ("PAX", "PA", [13:582, 584:772]),
    ("PB1", "PB1", [13:2286]),
    ("PB1-F2", "PB1", [221:379]),
    ("PB2", "PB2", [16:2295])
]

GENES = [i[1] for i in GENE_ORFS]

include("scripts/08_translate.jl")
println("Completed step 8 - Translate consensus to AA")

## Step 9: Pairwise differences
include("scripts/09_pairwise_diffs.jl")
println("Completed step 9 - Get pairwise identities")

## Step 10: Initial, crude phylogeny - we do this in order to determine which
# other sequences we might want in our real phylogeny
include("scripts/10_phylo.jl")
println("Completed step 10 - Initial phylogeny")

## Step 11: Process the second round of references including more N5 references
include("scripts/11_process_ref2.jl")
println("Completed step 11 - Process second round of references")

## Step 12: Second round of phylogeny
include("scripts/12_phylo2.jl")
println("Completed step 12 - Second round of phylogeny")