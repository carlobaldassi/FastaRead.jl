# example 2: read with readentry
#            default return type (ASCIIString)

module FastaReadExample2
using FastaRead

function read_fasta_file(filename::String)
    fr = FastaReader(filename)

    while !eof(fr)
        name, seq = readentry(fr)
        println("num=$(fr.num_parsed) name=$name seq=$seq")
    end
    println("read $(fr.num_parsed) entries")
end

current_dir = dirname(Base.source_path())
read_fasta_file(joinpath(current_dir, "test.fasta.gz"))
end
