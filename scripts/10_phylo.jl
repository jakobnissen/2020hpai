maybedir("results/phylo")
maybedir("results/phylo/set1")
maybedir("results/log/phylo/set1")

# Step one: Get the references
TREE_SEGMENTS = let
    arr = copy(SEGMENTS)
    deleteat!(arr, findfirst(isequal("NA"), arr))
    append!(arr, ["N5", "N8"])
end

# Add all except N5, N8
records = Dict()
for segment in TREE_SEGMENTS
    in(segment, ["N5", "N8"]) && continue
    records[segment] = open(collect, FASTA.Reader, "results/ref/set1/$segment.fna")
    append!(records[segment], open(FASTA.Reader, "results/consensus/set1/$segment.fna") do reader
        records = []
        for record in reader
            header = String(record.data[first(record.identifier):max(last(record.identifier), last(record.description))])
            if occursin("/Denmark/13776-1", header) || occursin("Denmark/14138-1", header)
                continue
            end
            push!(records, record)
        end
        records
    end
    )
end

# Add N5, N8
for n in ["N5", "N8"]
    records[n] = FASTA.Record[]
end

# Add N5/N8 from ref
open(FASTA.Reader, "results/ref/set1/NA.fna") do reader
    for record in reader
        header = String(record.data[first(record.identifier):max(last(record.identifier), last(record.description))])
        m = match(r"^A/.+\((H\d+)(N\d+)\)$", header)
        m === nothing && error(header)
        in(m[2], ["N5", "N8"]) || continue
        push!(records[m[2]], record)
    end
end

# Add own sequences
open(FASTA.Reader, "results/consensus/set1/NA.fna") do reader
    for record in reader
        header = String(record.data[first(record.identifier):max(last(record.identifier), last(record.description))])

        # These are already in our references
        if occursin("/Denmark/13776-1", header) || occursin("Denmark/14138-1", header)
            continue
        end
        
        # Special case: This one was annotated without NA (which is N8)
        if header == "A/barnacle_goose/Denmark/14599-1/2020-11-07(H5)"
            push!(records["N8"], record)
            continue
        end

        m = match(r"^A/.+\((H\d+)(N\d+)\)$", header)
        m === nothing && error(header)
        in(m[2], ["N5", "N8"]) || continue
        push!(records[m[2]], record)
    end
end

# Write 'em down
for segment in TREE_SEGMENTS
    maybedir("results/phylo/set1/$segment")
    catpath = "results/phylo/set1/$segment/cat.fna"
    open(FASTA.Writer, catpath) do writer
        foreach(record -> write(writer, record), records[segment])
    end
end

Threads.@threads for segment in TREE_SEGMENTS
    # Step two: MAFFT align
    catpath = "results/phylo/set1/$segment/cat.fna"
    alnpath = "results/phylo/set1/$segment/cat.aln.fna"
    run(pipeline(`mafft --thread 2 $catpath`, stdout=alnpath))

    # Step three: IQTREE
    iqdir = "results/phylo/set1/$segment"
    iqfiles = readdir(iqdir)
    for file in iqfiles
        startswith(file, "iq") && rm(joinpath(iqdir, file))
    end
    iqpath = "$iqdir/iq"
    stdout = "results/log/phylo/set1/$segment.log"
    pipe = pipeline(`iqtree -safe -s $alnpath -pre $iqpath -nt 2 -m HKY+G2 -nm 2500 -bb 1000`, stdout=stdout)
    run(pipe)
end