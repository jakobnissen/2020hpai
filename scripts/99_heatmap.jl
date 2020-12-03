using Plots
using FASTX

mkdir("results/plots/heatmaps")

# Begin by aligning just the consensus sequences
SEGMENTS = ["HA", "NA", "MP", "NP", "PB1", "PB2", "PA", "NS"]
Threads.@threads for segment in SEGMENTS
    inpath = "results/consensus/set1/$segment.fna"
    alnpath = "results/plots/heatmaps/$segment.aln.fna"
    run(pipeline(`linsi --thread 2 $inpath`, stdout=alnpath))
end

# Get the alignment from the Phylo
function distance_matrix(alnpath::AbstractString)
    (sequences, names) = open(FASTA.Reader, alnpath) do reader
        names, sequences = [], []
        for record in reader
            push!(names, FASTA.identifier(record))
            push!(sequences, FASTA.sequence(record))
        end
        sequences, names
    end

    # Reorder cols here based on the plots I've observed
    order = [2:7 ; [10, 13, 14, 16,
    9, 15, 11,
    8, 12, 1]]
    names = names[order]
    sequences = sequences[order]
    #

    res = zeros(Int, (length(names), length(names)))
    for i in 1:length(names)-1
        for j in i+1:length(names)
            d = 0
            for (x, y) in zip(sequences[i], sequences[j])
                if x !== y && x !== DNA_Gap && y !== DNA_Gap
                    d += 1
                end
            end
            res[i, j] = d
            res[j, i] = d
        end
    end
    return res, names
end

function make_heatmap(matrix::Matrix{<:Number}, names)
    ticks = (eachindex(names), names)
    plt = heatmap(matrix, yticks=ticks, c=palette(:inferno, maximum(matrix)))
    return plt
end

for segment in SEGMENTS
    # The falcon is quite different in these genes. So we remove them from the heatmap.
    (matrix, names) = if !(segment in ["NA", "PA", "NP", "PB2"])
        matrix, names = distance_matrix("results/plots/heatmaps/$segment.aln.fna")
        (matrix, names)
    else
        matrix, names = distance_matrix("results/plots/heatmaps/$segment.aln.fna")
        bad_index = argmax(vec(sum(matrix, dims=2)))
        indices = [i for i in eachindex(names) if i != bad_index]
        matrix = matrix[:, indices][indices, :]
        names = names[indices]
        (matrix, names)
    end
    println(segment, " ", maximum(matrix))
    plt = make_heatmap(matrix, names)
    savefig(plt, "results/plots/heatmaps/$segment.png")
end