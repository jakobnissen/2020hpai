maybdir("results/trim")
fastq_regex = r"S-(.*)_L001_R[12]_001.fastq.gz"

files = readdir("results/raw_reads")
pairs = Dict()
for file in files
    m = match(fastq_regex, file)
    basename = m[1]
    push!(get!(pairs, basename, []), file)
end
@assert all(x -> length(x) == 2, values(pairs))
foreach(sort!, values(pairs))

# Trimming
Threads.@threads for (basename, (fw, rv)) in collect(pairs)
    fw = "results/raw_reads/$fw"
    rv = "results/raw_reads/$rv"
    bn = "results/trim/$basename"
    command = `AdapterRemoval --file1 $fw --file2 $rv --basename $bn --minlength 20 --trimns 
              --trimqualities --minquality 20 --qualitybase 33 --gzip --trimwindows 5 --threads 1`
    run(command)
end

# FastQC report
# 14139-4-01_S8 looks a little bad to be honest. But let's continue
mkdir("results/fastqc")
Threads.@threads for (basename, pair) in collect(pairs)
    mkdir("results/fastqc/$basename")
    fw = "results/trim/$(basename).pair1.truncated.gz"
    rv = "results/trim/$(basename).pair2.truncated.gz"
    run(`fastqc -t 1 $fw $rv -o results/fastqc/$basename`)
end

