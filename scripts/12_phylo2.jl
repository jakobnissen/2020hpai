maybedir("results/log/phylo/set2")

Threads.@threads for segment in TREE_SEGMENTS
    # Step 1: MAFFT align
    catpath = "results/phylo/set2/$segment/cat.fna"
    alnpath = "results/phylo/set2/$segment/cat.aln.fna"
    run(pipeline(`mafft --thread 2 $catpath`, stdout=alnpath))

    # Step 2: IQTREE
    iqdir = "results/phylo/set2/$segment"
    iqfiles = readdir(iqdir)
    for file in iqfiles
        startswith(file, "iq") && rm(joinpath(iqdir, file))
    end
    iqpath = "$iqdir/iq"
    stdout = "results/log/phylo/set2/$segment.log"
    pipe = pipeline(`iqtree -safe -s $alnpath -pre $iqpath -nt 2 -m HKY+G2 -nm 2500 -bb 1000`, stdout=stdout)
    run(pipe)
end