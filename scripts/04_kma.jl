maybedir("results/aln")
maybedir("results/aln/set1")
maybedir("results/log/aln1")

# Index set1
Threads.@threads for segment in SEGMENTS
    path = "results/ref/set1/$segment.fna"
    run(`kma index -k 14 -i $path -o $path`)
end

# Sparse mapping to find best templates
Threads.@threads for basename in basenames
    for segment in SEGMENTS
        fw = "results/trim/$basename.pair1.truncated.gz"
        rv = "results/trim/$basename.pair2.truncated.gz"
        out = "results/aln/set1/$(basename)_$(segment)"
        db = "results/ref/set1/$segment.fna"
        stderr = "results/log/aln1/$basename.spa"
        pipe = pipeline(`kma -ipe $fw $rv -o $out -t_db $db -t 1 -Sparse`, stderr=stderr)
        run(pipe)
    end
end

# Get best templates for each
bestindex = Dict()
for segment in SEGMENTS
    bestindex[segment] = Dict()
    for basename in basenames
        open("results/aln/set1/$(basename)_$(segment).spa") do file
            readline(file) # skip header
            fields = readline(file)
            bestnum = parse(Int, split(fields, '\t')[2])
            bestindex[segment][basename] = bestnum
        end
    end
end

# Now do second round of mapping
for basename in basenames
    Threads.@threads for segment in SEGMENTS
        fw = "results/trim/$basename.pair1.truncated.gz"
        rv = "results/trim/$basename.pair2.truncated.gz"
        out = "results/aln/set1/$(basename)_$(segment)"
        db = "results/ref/set1/$segment.fna"
        index = bestindex[segment][basename]
        run(`kma -ipe $fw $rv -o $out -t_db $db
        -t 2 -k 8 -gapopen -5 -nf -matrix -Mt1 $index`)
    end
end