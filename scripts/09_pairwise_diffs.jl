# Get pairwise diffs for all the H5N8 (all except the falcon)
function matches(seqa::NucleotideSeq, seqb::NucleotideSeq)
    seqa == seqb && return (length(seqa), length(seqa))
    model = AffineGapScoreModel(match=1, mismatch=-2, gap_open=-5, gap_extend=-1) # same as I used for KMA
    aln = alignment(pairalign(GlobalAlignment(), seqa, seqb, model))
    return count_matches(aln), length(aln)
end

println("Differences ranges on nt/aa level")
for included in [true, false]
    println("With H5N5 sample? ", included, '\n')

    println("Nucleotide")
    for segment in SEGMENTS
        (names, seqs) = open(FASTA.Reader, "results/consensus/set1/$segment.fna") do reader
            names, seqs = [], []
            for record in reader
                name = FASTA.identifier(record)
                if !included && name == "A/peregrine_falcon/Denmark/13776-1/2020-10-30(H5N5)"
                    continue
                end
                push!(names, name)
                push!(seqs, FASTA.sequence(LongDNASeq, record))
            end
            (names, seqs)
        end
        identities = Float64[]
        for i in 1:length(seqs)-1
            for j in i+1:length(seqs)
                m, len = matches(seqs[i], seqs[j])
                identity = m / len
                push!(identities, identity)
            end
        end
        println(segment, " ", round(100*minimum(identities), digits=3), " ", round(100*maximum(identities), digits=3))
    end
    println("")

    # Also for AAs
    println("Amino acid - discount NA if H5N5 is present")
    for gene in GENES
        seqs = []
        for basename in basenames
            !included && basename == "13776_S1" && continue
            seq = open(FASTA.Reader, "results/translated/set1/$(basename)_$gene.faa") do reader
                FASTA.sequence(LongAminoAcidSeq, iterate(reader)[1])
            end
            push!(seqs, seq)
        end

        identities = Float64[]
        for i in 1:length(seqs)-1
            for j in i+1:length(seqs)
                m = sum(i == j for (i, j) in zip(seqs[i], seqs[j]))
                identity = m / length(seqs[i])
                push!(identities, identity)
            end
        end
        println(gene, " ", round(100*minimum(identities), digits=3), " ", round(100*maximum(identities), digits=3))
    end
    println("")
end