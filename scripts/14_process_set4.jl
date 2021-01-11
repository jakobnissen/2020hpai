REMOVE_SET4 = Set([
    # These are too far away!
    "A/Whooper swan/Xinjiang/1/2020(H5N6)",
    "A/Mute swan/Xinjiang/5/2020(H5N6)",
    "A/Whooper swan/Xinjiang/13/2020(H5N6)",

    # Not sure about these - check tree. I think they are too far in N5 tree
    "A/mallard_duck/Netherlands/83/2008(H12N5)",
    "A/Anas_platyrhynchos/Belgium/398_H189467/2017(H10N5)",
    "A/mallard/Jiangxi/8341/2004(H6N5)",
    "A/duck/Hunan/748/2005(H6N5)",
])

function getheader(record::FASTA.Record)
    ind = first(record.identifier):max(last(record.identifier), last(record.description))
    String(record.data[ind])
end

## Add set4
bysegment = Dict(s => FASTA.Record[] for s in TREE_SEGMENTS)
for segment in TREE_SEGMENTS
    #maybedir("results/phylo/set3/$segment")
    open(FASTA.Reader, "results/phylo/set3/$segment/cat.fna") do reader
        for record in reader
            header = getheader(record)
            if in(header, REMOVE_SET4)
                continue
            end
            push!(bysegment[segment], record)
        end
    end
end

## Add in new segments
open("ref/set4/set4.fna.gz") do file
    reader = FASTA.Reader(GzipDecompressorStream(file))

    for record in reader
        header = getheader(record)
        # A/Great Black-backed Gull/Netherlands/1/2017  HA  2017-12-18 | A / H5N6
        fields = split(header, "  ")
        identifier = String(strip(fields[1]))
        segment = String(fields[2])
        kind = String(last(split(header)))
        if segment == "NA"
            segment = String(kind[3:4])
        end
        newheader = "$identifier($kind)"
        newrecord = FASTA.Record(newheader, extract_sequence_data(record))

    end
end

# Write
for segment in TREE_SEGMENTS
    maybedir("results/phylo/set4/$segment")
    open(FASTA.Writer, "results/phylo/set4/$segment/cat.fna") do writer
        foreach(record -> write(writer, record), bysegment[segment])
    end
end

maybedir("results/log/phylo/set4")

Threads.@threads for segment in TREE_SEGMENTS
    # Step 1: MAFFT align
    catpath = "results/phylo/set4/$segment/cat.fna"
    alnpath = "results/phylo/set4/$segment/cat.aln.fna"
    run(pipeline(`mafft --thread 2 $catpath`, stdout=alnpath))

    # Step 2: IQTREE
    iqdir = "results/phylo/set4/$segment"
    iqfiles = readdir(iqdir)
    for file in iqfiles
        startswith(file, "iq") && rm(joinpath(iqdir, file))
    end
    iqpath = "$iqdir/iq"
    stdout = "results/log/phylo/set4/$segment.log"
    pipe = pipeline(`iqtree -safe -s $alnpath -pre $iqpath -nt 2 -m HKY+G2 -nm 2500 -bb 1000`, stdout=stdout)
    run(pipe)
end

###
segmnts = Dict()
for segment in SEGMENTS
    segmnts[segment] = []
    for basename in basenames
        open(FASTA.Reader, "results/consensus/set1/$(basename)_$(segment).fna") do reader
            rec, _ = iterate(reader)
            push!(segmnts[segment], (getheader(rec), FASTA.sequence(LongDNASeq, rec)))
        end
    end
end
