require("fastaread.jl")
using FastaRead

function read_fasta_file(filename::String)
    fr = FastaReader(filename)

    n = 0
    for (name, seq) in fr
        n += 1
        println("num=$n name=$name seq=$(ascii(seq))")
    end
    println("read $n entries")
end

read_fasta_file("test.fasta.gz")
