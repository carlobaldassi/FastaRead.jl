# example 4: read with iterators,
#            custom return type

module FastaReadExample4
using FastaRead

function read_fasta_file(filename::String)
    # We set the reader so that it returns a Vector{Char};
    # If you have BioSeq package, you can use more
    # specialized uints, e.g. Vector{AminoAcid}
    fr = FastaReader{Vector{Char}}(filename)

    for (name, seq) in fr
        println("num=$(fr.num_parsed) name=$name seq=$seq")
    end
    println("read $(fr.num_parsed) entries")
end

current_dir = dirname(Base.source_path())
read_fasta_file(joinpath(current_dir, "test.fasta.gz"))
end
