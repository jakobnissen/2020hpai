maybedir("results/phylo/set6")

# Load single new HA
new_record = open(x -> first(iterate(x)), FASTA.Reader, "ref/set6/HA.fna")

function rename_record(rec::FASTA.Record)
    patterns = [
        "(H5N5)(H5N5)"=>s"_H5N5",
        "[H5N5][H5N5]"=>s"_H5N5",
        "[H5N8][H5N8]"=>s"_H5N8",
        "__"=>"_",
    ]
    header = foldl(replace, patterns, init=getheader(rec))
    header = replace(header, r"\[(H\dN\d)\]"=>s"_\1")
    return FASTA.Record(header, FASTA.sequence(rec))
end

cp("results/phylo/set4/N5/cat.fna", "results/phylo/set5/N5/cat.fna", force=true)

for segment in TREE_SEGMENTS
    dst = "results/phylo/set6/$segment/cat.fna"
    src = "results/phylo/set5/$segment/cat.fna"
    reader = open(FASTA.Reader, src)
    writer = open(FASTA.Writer, dst)

    for record in reader
        write(writer, rename_record(record))
    end

    if segment == "HA"
        write(writer, FASTA.Record("A/Duck/Egypt/SMG4/2019_H5N8", FASTA.sequence(new_record)))
    end
    close(reader)
    close(writer)
end

Threads.@threads for segment in TREE_SEGMENTS
    segment == "N5" && continue
    # Step 1: MAFFT align
    catpath = "results/phylo/set6/$segment/cat.fna"
    alnpath = "results/phylo/set6/$segment/cat.aln.fna"
    run(pipeline(`mafft --thread 2 $catpath`, stdout=alnpath))
end

maybedir( "results/log/phylo/set6")

Threads.@threads for segment in TREE_SEGMENTS
    segment == "N5" && continue
    alnpath = "results/phylo/set6/$segment/cat.aln.fna"
    # Step 2: IQTREE
    iqdir = "results/phylo/set6/$segment"
    iqfiles = readdir(iqdir)
    for file in iqfiles
        startswith(file, "iq") && rm(joinpath(iqdir, file))
    end
    iqpath = "$iqdir/iq"
    stdout = "results/log/phylo/set6/$segment.log"
    pipe = pipeline(`iqtree -safe -s $alnpath -pre $iqpath -nt 2 -m HKY+G2 -nm 2500 -bb 1000`, stdout=stdout)
    run(pipe)
end