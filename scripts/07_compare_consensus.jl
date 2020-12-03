# These have already been uploaded to GISAID
duplicates = ["A/barnacle_goose/Denmark/14138-1/2020-11-04(H5N8)", "A/peregrine_falcon/Denmark/13776-1/2020-10-30(H5N5)"]

function getheader(record::FASTA.Record)
    mn = first(record.identifier)
    mx = max(last(record.identifier), last(record.description))
    return String(record.data[mn:mx])
end

ref_dict = Dict(s=>Dict() for s in SEGMENTS)
own_dict = Dict(s=>Dict() for s in SEGMENTS)
for segment in SEGMENTS
    open(FASTA.Reader, "results/ref/set1/$segment.fna") do reader
        for record in reader
            if in(getheader(record), duplicates)
                ref_dict[segment][getheader(record)] = FASTA.sequence(record)
            end
        end
    end
end

for segment in SEGMENTS
    open(FASTA.Reader, "results/consensus/set1/$segment.fna") do reader
        for record in reader
            if in(getheader(record), duplicates)
                own_dict[segment][getheader(record)] = FASTA.sequence(record)
            end
        end
    end
end

# Now compare!
println("Differences from uploaded consensus")
model = AffineGapScoreModel(match=1, mismatch=-2, gap_open=-5, gap_extend=-1) # same as KMA
for segment in SEGMENTS
    for sample in keys(ref_dict[segment])
        seq_a = ref_dict[segment][sample]
        seq_b = own_dict[segment][sample]

        aln = alignment(pairalign(GlobalAlignment(), seq_a, seq_b, model))
        println(count_mismatches(aln), " ", segment, " ", sample)
    end
end