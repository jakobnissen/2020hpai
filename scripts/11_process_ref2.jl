maybedir("results/phylo/set2")

# Add set1 first
for segment in TREE_SEGMENTS
    maybedir("results/phylo/set2/$segment")
    cp("results/phylo/set1/$segment/cat.fna", "results/phylo/set2/$segment/cat.fna", force=true)
end

# Then add set2
bysegment = open(F"ref/set2/ref.fna.gz") do file
    reader = FASTA.Reader(GzipDecompressorStream(file))
    records = Dict(s => FASTA.Record[] for s in TREE_SEGMENTS)
    for record in reader
        identifier, segment, htype, ntype = extract_header_data(record)
        key = segment == "NA" ? ("N" * ntype) : segment
        seq = extract_sequence_data(record)
        newrecord = FASTA.Record("$identifier(H$(htype)N$(ntype))", seq)
        push!(records[key], newrecord)
    end
    records
end

# And N5
n5_header_regex = r"^(A/[^\|]+)\|A_/_(H\d+N5)\|NA\|.+"
open(FASTA.Reader, "ref/set2/n5.fna.gz") do file
    reader = FASTA.Reader(GzipDecompressorStream(file))
    for record in reader
        header = String(record.data[first(record.identifier):max(last(record.identifier), last(record.description))])
        m = match(n5_header_regex, header)
        m === nothing && error(header)
        identifier, hn = m[1], m[2]
        seq = extract_sequence_data(record)
        newrecord = FASTA.Record("$identifier($hn)", seq)
        push!(bysegment["N5"], newrecord)
    end
end

# Cat these onto the set2 ones we copied at the top
for segment in TREE_SEGMENTS
    open(FASTA.Writer, "results/phylo/set2/$segment/cat.fna", append=true) do writer
        foreach(record -> write(writer, record), bysegment[segment])
    end
end
