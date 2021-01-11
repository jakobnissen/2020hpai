const alignment_model = AffineGapScoreModel(match=1, mismatch=-2, gap_open=-5, gap_extend=-1) # same as I used for KMA

# Relative to ref/relative.fna (which is 14535-1 annotated by a tool)
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

function get_map(ref::LongDNASeq, query::LongDNASeq)
    aln = alignment(pairalign(GlobalAlignment(), query, ref, alignment_model))
    refgap = LongDNASeq([r for (q,r) in aln])
    qgap = LongDNASeq([q for (q,r) in aln])
    
    map = Vector{Union{Nothing, DNA}}(undef, length(refgap))
    refind = 0
    for i in eachindex(refgap)
        if refgap[i] != DNA_Gap
            refind += 1
        else
            continue
        end
        if qgap[i] == DNA_Gap
            map[refind] = nothing
        else
            map[refind] = qgap[i]
        end
    end
    return map
end

function get_orf(ref::LongDNASeq, query::LongDNASeq, reforfs::Vector{<:UnitRange})
    map = get_map(ref, query)
    seqbuffer = DNA[]
    for reforf in reforfs
        for i in reforf
            if map[i] === nothing
                error("Pos $i translates to gap in ref seq")
            end
            push!(seqbuffer, map[i])
        end
    end
    return LongDNASeq(seqbuffer)
end

function get_aa(ref::LongDNASeq, query::LongDNASeq, reforfs::Vector{<:UnitRange})
    orf = get_orf(ref, query, reforfs)
    aas = BioSequences.translate(orf)
    if last(aas) == AA_Term
        aas = aas[1:end-1]
    end
    return aas
end

maybedir("results/translated")
maybedir("results/translated/set1")

refseqs = Dict()
for segment in SEGMENTS
    open("ref/relative.fna.gz") do file
        reader = FASTA.Reader(GzipDecompressorStream(file))
        for record in reader
            refseqs[FASTA.identifier(record)] = FASTA.sequence(LongDNASeq, record)
        end
    end
end

for segment in SEGMENTS
    for basename in basenames
        record = open(FASTA.Reader, "results/consensus/consensus/$(basename)_$(segment).fna") do reader
            iterate(reader)[1]
        end
        seq = FASTA.sequence(LongDNASeq, record)
        
        for (gene, smt, orfs) in GENE_ORFS
            if segment == smt
                aa = if occursin("13776", basename) && segment == "NA"
                    refs = replace(replace("ATGAATCCAAATCAGAAAATAATAACAATTGGCTCAATATCATTAGGATTGGTTGTTTTCAACATTCTTC
                    TTCATGTTGCATCAATAGTCCTAGGGATAATATCAGTGACCAAAGACCACGAAGCATACACATGCAACAC
                    AACTGAGGTGTACAATGAGACAGTGAGGGTGGAAACAGTAACTATCCCTGTCAACAACACTATTTACATA
                    GAAAGGGAATTGACCCATGAACCAGAATTTCTTAACAACACGGAGCCTCTCTGTGAGGTATCAGGATTTG
                    CCATTGTTTCCAAAGACAATGGGATCAGAATAGGTTCAAGAGGGCATGTCTTCGTCATAAGAGAGCCTTT
                    TGTGGCTTGTGGTCCTTCAGAGTGCAGGACATTTTTCTTAACTCAAGGCGCTCTGTTGAATGATAAGCAT
                    TCAAACAATACAGTAAAAGACAGGAGTCCCTATCGAGCTTTAATGAGCGTGCCATTGGGATCCTCTCCCA
                    ATGCTTACCAAGCCAAATTTGAGTCCGTCGGATGGTCTGCTACAGCCTGCCATGACGGGAAGGAGTGGAT
                    GGCTATTGGGGTGAGTGGTGCAGACGATGATGCCTATGCCGTCATCCATTATGGAGGGATACCAACAGAC
                    GTAGTGAGGTCATGGAGAAAGCAAATACTGAGGACACAAGAGTCTTCATGTGTTTGCATGAAAGGAGAGT
                    GTTATTGGGTAATGACAGACGGTCCAGCAAATAACCAAGCAAGTTACAAAATCTTTAAATCACAGAAAGG
                    ACTAGTTGTAGATGAAAAGGAAATTTCATTTCAAGGTGGACACATAGAGGAGTGTTCCTGCTATCCCAAT
                    ATGGGGAAAGTGGAATGCGCCTGCAGGGATAATTGGAACGGAATGAATAGGCCAATTCTAACATTCGATG
                    AAAACCTTGAATATGAGGTTGGTTATCTGTGTGCTGGGATTCCAACAGATACTCCGCGAGTTCAAGACAG
                    CAGTTTCACTGGCTCATGCACCAATGCTGTTGGGGGAAGTGGAACAAATAATTATGGAGTGAAGGGATTT
                    GGTTTTAGGCAAGGGACCAGTGTATGGGCAGGAAGGACAATAAGCACTTCTTCTCGGAGTGGCTTCGAGG
                    TCTTACTAATAGAAAATGGATGGGTTAGGCCAAGTAAAACCATTAGCAAAAAAGTTGAGGTTTTGAATAA
                    TAAGAACTGGTCAGGATACAGTGGGTCCTTCACCATCCCCACTACGATGACAAGCAAGAGCTGTCTGGTT
                    CCATGCTTTTGGCTGGAAATGATTCGAGGTAAACCGGAAGAGAGAACTAGCATTTGGACCTCAAGTAGTT
                    CCACTGTGTTTTGTGGTGTTTCTAGTGAGGTCCCAGGATGGTCCTGGGATGATGGAGCAATTCTACCATT
                    CGACATCGATAAGATGTAG", " "=>""), "\n"=>"")
                    ref = LongDNASeq(refs)
                    get_aa(ref, seq, [1:1419])
                else
                    ref = refseqs[segment]
                    get_aa(ref, seq, orfs)
                end
                open(FASTA.Writer, "results/translated/set1/$(basename)_$(gene).faa") do writer
                    newrecord = FASTA.Record(FASTA.identifier(record), aa)
                    write(writer, newrecord)
                end
            end
        end
    end
end

for gene in GENES
    open(FASTA.Writer, "results/translated/set1/$gene.faa") do writer
        for basename in basenames
            record, _ = open(iterate, FASTA.Reader, "results/translated/set1/$(basename)_$gene.faa")
            write(writer, record)
        end
    end
end