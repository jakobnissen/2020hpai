REMOVE_SET3 = Set([
    "A/Whooper swan/Xinjiang/6/2020(H5N6)",
    "A/Whooper swan/Xinjiang/7/2020(H5N6)",
    "A/Whooper swan/Xinjiang/8/2020(H5N6)",
    "A/Whooper swan/Xinjiang/9/2020(H5N6)",
    "A/Whooper swan/Xinjiang/10/2020(H5N6)",
    "A/Whooper swan/Xinjiang/11/2020(H5N6)",
    "A/Whooper swan/Xinjiang/12/2020(H5N6)",
    "A/Whooper swan/Xinjiang/2/2020(H5N6)",
    "A/Whooper swan/Xinjiang/3/2020(H5N6)",
    "A/Mute swan/Xinjiang/4/2020(H5N6)",

    "A/chicken/Poland/004/2020(H5N8)",
    "A/chicken/Poland/003/2020(H5N8)",
    "A/turkey/Poland/027/2020(H5N8)",
    "A/domestic_goose/Poland/028/2020(H5N8)",
    "A/chicken/Poland/054/2020(H5N8)",
    "A/turkey/Poland/079/2020(H5N8)",
    "A/laying_hen/Poland/095/2020(H5N8)",
    "A/turkey/Poland/096/2020(H5N8)",
    "A/turkey/Poland/182/2020(H5N8)",
    "A/domestic_duck/Poland/219/2020(H5N8)",
    "A/domestic_duck/Poland/221/2020(H5N8)",
    
    "A/domestic_duck/Poland/223/2020(H5N8)",
    "A/domestic_duck/Poland/229/2020(H5N8)",
    "A/domestic_duck/Poland/230/2020(H5N8)",
    "A/domestic_duck/Poland/237/2020(H5N8)",
    "A/turkey/Poland/366/2020(H5N8)",
    "A/domestic_goose/Poland/274/2020(H5N8)",
    "A/domestic_duck/Poland/285/2020(H5N8)",
    "A/domestic_duck/Poland/263/2020(H5N8)",
    "A/domestic_duck/Poland/271/2020(H5N8)",

    "A/Eurasian Wigeon/Netherlands/4/2020(H5N1)",
    "A/Eurasian Wigeon/Netherlands/5/2020(H5N1)",
    "A/Eurasian Wigeon/Netherlands/7/2020(H5N8)",
    "A/greylag goose/Netherlands/20016582-004/2020(H5N1)",
    "A/eurasian teal/Netherlands/20016896-013/2020(H5N1)",
    "A/barnacle goose/Netherlands/20016935-002/2020(H5N8)",
    "A/chicken/Netherlands/20017639-001/2020(H5N8)",
    "A/chicken/Netherlands/20017694-004/2020(H5N8)",

    "A/guinea fowl/Germany-NW/AI01184/2020(H5N8)"
    ])

function getheader(record::FASTA.Record)
    ind = first(record.identifier):max(last(record.identifier), last(record.description))
    String(record.data[ind])
end

maybedir("results/phylo/set3")

# Add set1 first
bysegment = Dict(s => FASTA.Record[] for s in TREE_SEGMENTS)
for segment in TREE_SEGMENTS
    maybedir("results/phylo/set2/$segment")
    open(FASTA.Reader, "results/phylo/set2/$segment/cat.fna") do reader
        for record in reader
            header = getheader(record)
            if in(header, REMOVE_SET3)
                continue
            end
            push!(bysegment[segment], record)
        end
    end
end

# Add the Guangdong one
# A/Goose/Guangdong/1/96 |  PB2 | A / H5N1 | 1996-01-01 | EPI_ISL_1254 | A/Goose/Guangdong/1/96 | 5779
open("ref/set3/guangdong.fna.gz") do file
    reader = FASTA.Reader(GzipDecompressorStream(file))
    for record in reader
        header = getheader(record)
        header_regex = r"^(A/[^\|]+)\|  ([\d\w]+) \|"
        m = match(header_regex, header)
        m === nothing && error(header)
        identifier = String(strip(m[1]))
        segment = m[2]
        newrecord = FASTA.Record("$identifier(H5N1)", extract_sequence_data(record))
        if segment != "NA"
            push!(bysegment[segment], newrecord)
        end
    end
end

# Write
for segment in TREE_SEGMENTS
    maybedir("results/phylo/set3/$segment")
    open(FASTA.Writer, "results/phylo/set3/$segment/cat.fna") do writer
        foreach(record -> write(writer, record), bysegment[segment])
    end
end

maybedir("results/log/phylo/set3")

Threads.@threads for segment in TREE_SEGMENTS
    # Step 1: MAFFT align
    catpath = "results/phylo/set3/$segment/cat.fna"
    alnpath = "results/phylo/set3/$segment/cat.aln.fna"
    run(pipeline(`mafft --thread 2 $catpath`, stdout=alnpath))

    # Step 2: IQTREE
    iqdir = "results/phylo/set3/$segment"
    iqfiles = readdir(iqdir)
    for file in iqfiles
        startswith(file, "iq") && rm(joinpath(iqdir, file))
    end
    iqpath = "$iqdir/iq"
    stdout = "results/log/phylo/set3/$segment.log"
    pipe = pipeline(`iqtree -safe -s $alnpath -pre $iqpath -nt 2 -m HKY+G2 -nm 2500 -bb 1000`, stdout=stdout)
    run(pipe)
end
