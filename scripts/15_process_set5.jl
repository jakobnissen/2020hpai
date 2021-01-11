REMOVE_SET5 = Set([
    # Too far away
    "A/chicken/Bulgaria/77_20VIR1727/2020(H5N2)"
])

maybedir("results/phylo/set5")

# Add set 4
bysegment = Dict(s => FASTA.Record[] for s in TREE_SEGMENTS)
for segment in TREE_SEGMENTS
    maybedir("results/phylo/set5/$segment")
    open(FASTA.Reader, "results/phylo/set4/$segment/cat.fna") do reader
        for record in reader
            header = getheader(record)
            startswith(header, ">") && error(header)
            if in(header, REMOVE_SET5)
                println("REMOVED *",header,"*")
                continue
            end
            push!(bysegment[segment], record)
        end
    end
end

## Add in new segments
# Looks like this: A/pigeon/Italy/17VIR9842/2017 HA  2017-10-23 | A / H5N8
# or this MW026118.1 Influenza A virus (A/mute swan/Denmark/19192-1/2016(H5N8)) [ ... ]
# or this A_pheasant_Denmark_12106-3_2018-08-27(H5N6)HA
# Should look like A/chicken/Czech Republic/1175-1/2020(H5N8)
for segment in SEGMENTS
    smt = segment == "NA" ? "N8" : segment
    existing = Set(getheader(r) for r in bysegment[smt])
    open(FASTA.Reader, "ref/set5/Add_$(segment)_seq_16_12_20.fasta") do reader
        for record in reader
            header = getheader(record)
            #occursin(segment, header) || (println(header); error())
            identifier = ""
            m = match(r"^A/", header)
            if m !== nothing
                before, after = split(header, segment, limit=2)
                identifier = String(strip(before))
                subtype = String(last(split(header)))
            end
            m2 = match(r"^[\w\d\.]+ Influenza A virus \((A/[^\(]+)\((H\d+N\d+)\)\)", header)
            if m2 !== nothing
                identifier, subtype = m2.captures
            end
            isempty(identifier) && (println(header); error())

            newheader = "$identifier($subtype)"
            in(newheader, existing) && continue

            push!(bysegment[smt], FASTA.Record(newheader, FASTA.sequence(LongDNASeq, record)))
        end
    end
end

# Write
for segment in TREE_SEGMENTS
    segment == "N5" && continue # we dont have any of these
    maybedir("results/phylo/set5/$segment")
    open("results/phylo/set5/$segment/cat.fna", "w") do file
        foreach(bysegment[segment]) do record
            newheader = getheader(record) |> (s -> replace(s, ' '=>'_')) |> (s -> replace(s, ')'=>']')) |>
            (s -> replace(s, '('=>'['))
            println(file, ">", newheader)
            println(file, FASTA.sequence(LongDNASeq, record))
        end
    end
end

maybedir("results/log/phylo/set5")

Threads.@threads for segment in TREE_SEGMENTS
    segment == "N5" && continue
    # Step 1: MAFFT align
    catpath = "results/phylo/set5/$segment/cat.fna"
    alnpath = "results/phylo/set5/$segment/cat.aln.fna"
    run(pipeline(`mafft --thread 2 $catpath`, stdout=alnpath))
end

Threads.@threads for segment in TREE_SEGMENTS
    segment == "N5" && continue
    alnpath = "results/phylo/set5/$segment/cat.aln.fna"
    # Step 2: IQTREE
    iqdir = "results/phylo/set5/$segment"
    iqfiles = readdir(iqdir)
    for file in iqfiles
        startswith(file, "iq") && rm(joinpath(iqdir, file))
    end
    iqpath = "$iqdir/iq"
    stdout = "results/log/phylo/set5/$segment.log"
    pipe = pipeline(`iqtree -safe -s $alnpath -pre $iqpath -nt 2 -m HKY+G2 -nm 2500 -bb 1000`, stdout=stdout)
    run(pipe)
end