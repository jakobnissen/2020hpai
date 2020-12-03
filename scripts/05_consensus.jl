function is_possibility(needle::NucleotideSeq, haystack::NucleotideSeq)
    @inbounds for i in eachindex(needle)
        iscompatible(needle[i], haystack[i]) || return false
    end
    return true
end

function slide!(primer::NucleotideSeq, seq::NucleotideSeq, minlen::Int)
    length(seq) < length(primer) && return 0
    primer = copy(primer)
    seq = copy(seq)[1:length(primer)]
    length(primer) == length(seq) || throw(ArgumentError("Must match in lengths"))
    for overlap in length(primer):-1:minlen
        is_possibility(primer, seq) && return overlap
        popfirst!(primer)
        pop!(seq)
    end
    return 0
end

function getheader(record::FASTA.Record)
    mn = first(record.identifier)
    mx = max(last(record.identifier), last(record.description))
    return String(record.data[mn:mx])
end

function load_primers(path::String)
    open(FASTA.Reader, path) do reader
        map(record -> (getheader(record), FASTA.sequence(LongDNASeq, record)), reader)
    end
end

function load_consensus(path::String)
    open(FASTA.Reader, path) do reader
        record, _ = iterate(reader)
        seq = FASTA.sequence(LongDNASeq, record)
        header = getheader(record)
        iterate(reader) === nothing || error("Multiple records in file $path")
        (header, seq)
    end
end

function remove_primers(seq::NucleotideSeq, primers::Vector{<:Tuple{String, NucleotideSeq}}, minlength::Int)
    overlaps = [slide!(p, seq, minlength) for (h,p) in primers]
    if maximum(overlaps) > 0
        i = argmax(overlaps)
        overlap = overlaps[i]
        header = primers[i][1]
        println(stdout, "Primer $header found with overlap $overlap")
        seq = seq[1+overlap:end]
    end
    return seq
end

function remove_primers(primerpath::String, consensuspath::String, output::String, minlength::Int)
    primers = load_primers(primerpath)
    header, seq = load_consensus(consensuspath)

    if !isempty(primers)
        seq = remove_primers(seq, primers, minlength)
        seq = remove_primers(reverse_complement(seq), primers, minlength)
        seq = reverse_complement(seq)
    end

    open(FASTA.Writer, output) do writer
        write(writer, FASTA.Record(header, seq))
    end 
end

renames = Dict(
        "13776_S1" => "A/peregrine falcon/Denmark/13776-1/2020-10-30(H5N5)",
        "14819-6_S13" => "A/chicken/Denmark/14819-6/2020-11-15(H5N8)",
        "14534-1-01_S4" => "A/barnacle goose/Denmark/14534-1/2020-11-04(H5N8)",
        "14139-3-01_S7" => "A/barnacle goose/Denmark/14139-3/2020-11-04(H5N8)",
        "14537-1-01_S7" => "A/barnacle goose/Denmark/14537-1/2020-11-04(H5N8)",
        "14139-2-01_S6" => "A/barnacle goose/Denmark/14139-2/2020-11-04(H5N8)",
        "14536-1-01_S6" => "A/barnacle goose/Denmark/14536-1/2020-11-05(H5N8)",
        "14139-1-01_S5" => "A/barnacle goose/Denmark/14139-1/2020-11-04(H5N8)",
        "14599-1-01_S10" => "A/barnacle goose/Denmark/14599-1/2020-11-07(H5)",
        "14600-1-01_S11" => "A/common_buzzard/Denmark/14600-1/2020-11-06(H5N8)",
        "14596-1-01_S9" => "A/peregrine falcon/Denmark/14596-1/2020-11-06(H5N8)",
        "14138-1-01_S4" => "A/barnacle goose/Denmark/14138-1/2020-11-04(H5N8)",
        "14600-2-01_S12" => "A/barnacle goose/Denmark/14600-2/2020-11-07(H5N8)",
        "14535-1-01_S5" => "A/graylag_goose/Denmark/14535-1/2020-11-09(H5N8)",
        "14139-4-01_S8" => "A/black-headed gull/Denmark/14139-4/2020-11-04(H5N8)",
        "14538-1-01_S8" => "A/barnacle goose/Denmark/14538-1/2020-11-05(H5N8)"
)

function rename_ownseq(basename::AbstractString, rec::FASTA.Record)
    newheader = replace(renames[basename], " "=>"_")
    seq = FASTA.sequence(rec)
    return FASTA.Record(newheader, seq)
end

primers = "ref/primers.fna"
maybedir("results/consensus")
maybedir("results/consensus/set1")

# Remove primers from segments
for segment in SEGMENTS
    for basename in basenames
        consensus = "results/aln/set1/$(basename)_$(segment).fsa"
        output = "/tmp/$(basename)_$(segment).fna"
        remove_primers(primers, consensus, output, 7)

        record = open(FASTA.Reader, output) do reader
            iterate(reader)[1]
        end
        renamed = rename_ownseq(basename, record)
        open(FASTA.Writer, "results/consensus/set1/$(basename)_$(segment).fna") do writer
            write(writer, renamed)
        end
    end
end

# Cat consensus together
for segment in SEGMENTS
    records = []
    for basename in basenames
        record = open(FASTA.Reader, "results/consensus/set1/$(basename)_$(segment).fna") do reader
            record = iterate(reader)[1]
        end
        push!(records, record)
    end
    open(FASTA.Writer, "results/consensus/set1/$segment.fna") do writer
        for record in records
            write(writer, record)
        end
    end
end