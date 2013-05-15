# example 1: read with iterators,
#            default return type (ASCIIString)

module FastaReadExample1
using FastaRead

function read_fasta_file(filename::String)
    fr = FastaReader(filename)

    for (name, seq) in fr
        println("num=$(fr.num_parsed) name=$name seq=$seq")
    end
    println("read $(fr.num_parsed) entries")
end

current_dir = dirname(Base.source_path())
read_fasta_file(joinpath(current_dir, "test.fasta.gz"))
end
