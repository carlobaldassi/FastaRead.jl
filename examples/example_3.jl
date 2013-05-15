# example 3: read with readall
#            default return type (ASCIIString)

module FastaReadExample3
using FastaRead

function read_fasta_file(filename::String)
    fr = FastaReader(filename)

    out = readall(fr)
    n = 0
    for (name, seq) in out
        n += 1
        println("num=$n name=$name seq=$seq")
    end
    # here of course n == fr.num_parsed == length(out)
    println("read $n entries")
end

current_dir = dirname(Base.source_path())
read_fasta_file(joinpath(current_dir, "test.fasta.gz"))
end
