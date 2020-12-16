function extract_header_data(record::FASTA.Record)
    header_regex = r"^(A/[^\|]+)\|  ([\d\w]+) \| A / H(\d+)N(\d+).+"
    header = String(record.data[first(record.identifier):last(record.description)])
    m = match(header_regex, header)
    isnothing(m) && throw(ArgumentError("Did not match header regex: " * header))

    return Tuple([String(strip(m[i])) for i in 1:4])
end

function extract_sequence_data(record::FASTA.Record)
    seqstr = String(record.data[record.sequence])

    # Special cases:
    if seqstr[1:18] == "achickenpolandhnmp"
        seqstr = seqstr[19:end]
    end

    # Replace X in sequence
    seqstr = replace(seqstr, "x"=>"n")
    return LongDNASeq(seqstr)
end

# Split by segment and rename
bysegment = open("ref/set1/set1.fna.gz") do file
    reader = FASTA.Reader(GzipDecompressorStream(file))
    bysegment = Dict()
    for record in reader
        identifier, segment, htype, ntype = extract_header_data(record)
        newname = "$identifier(H$(htype)N$(ntype))"
        seq = extract_sequence_data(record)
        newrecord = FASTA.Record(newname, seq)
        push!(get!(bysegment, segment, []), newrecord)
    end
    return bysegment
end

# One segment in each sample
@assert length(Set(map(length, values(bysegment)))) == 1

# Print to files
maybedir("results/ref")
maybedir("results/ref/set1")
for (segment, records) in bysegment
    open(FASTA.Writer, "results/ref/set1/$segment.fna") do writer
        foreach(x -> write(writer, x), records)
    end
end
