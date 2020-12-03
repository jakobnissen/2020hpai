maybedir("results/plots")
maybedir("results/plots/depths")

function get_depths(matrixpath::String)
    result = Int[]
    decompressor = GzipDecompressorStream(open(matrixpath))
    header = readline(decompressor)
    for line in eachline(decompressor)
        isempty(strip(line)) && continue
        fields = split(line)
        fields[1] == "-" && continue
        depth = sum(map(x -> parse(Int, x), fields[2:5]))
        push!(result, depth)
    end
    close(decompressor)
    return result
end

for basename in basenames
    isdir("results/plots/depths") || mkdir("results/plots/depths")
    plt = plot()
    for segment in SEGMENTS
        depths = log10.(get_depths("results/aln/set1/$(basename)_$(segment).mat.gz"))
        xs = range(0.0, 1.0; length=length(depths))
        plot!(plt, xs, depths, label=segment, ylims=(0, 5), legend=:bottomright)
    end
    savefig(plt, "results/plots/depths/$basename.pdf")
end
