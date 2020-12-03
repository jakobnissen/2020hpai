mkdir("results/raw_reads")
fastq_regex = r"S-(.*)_L001_R[12]_001.fastq.gz"

# We just put all the samples together as one
for ngs in readdir("raw")
    files = readdir("raw/$ngs")
    fastqs = filter(name -> !isnothing(match(fastq_regex, name)), files)
    for fastq in fastqs
        # Skip these samples - they're too weak to be of any use.
        startswith(fastq, "S-12947") && continue
        startswith(fastq, "S-14320") && continue
        run(`ln -s $(pwd())/raw/$ngs/$fastq results/raw_reads/$fastq`)
    end
end