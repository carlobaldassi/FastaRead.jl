using FastaRead

let
    out_ascii = {
        ("C0Z3L5_CAEEL/75-127",
         ".................................s-CTLPR...QIG." *
         "....T...GP.....Y..................RI.......PRW." *
         ".YYN............PVRGRCE...LF..Y.W...S.........." *
         ".......G.C...C.G...N...G.......NNFQTFQT....CQST" *
         "C-e......."),
        ("A8X0E5_CAEBR/76-128",
         ".................................s-CTLPR...QIG." *
         "....T...GP.....Y..................RI.......PRW." *
         ".YYN............PVRGRCE...LF..Y.W...S.........." *
         ".......G.C...C.G...N...G.......NNFQTFQT....CQST" *
         "C-e......."),
        ("C0Z3L4_CAEEL/75-127",
         ".................................s-CTLPR...QIG." *
         "....T...GP.....Y..................RI.......PRW." *
         ".YYN............PVRGRCE...LF..Y.W...S.........." *
         ".......G.C...C.G...N...G.......NNFQTFQT....CQST" *
         "C-e......."),
        ("Q9XWX5_CAEEL/75-127",
         ".................................s-CTLPR...QIG." *
         "....T...GP.....Y..................RI.......PRW." *
         ".YYN............PVRGRCE...LF..Y.W...S.........." *
         ".......G.C...C.G...N...G.......NNFQTFQT....CQST" *
         "C-e.......")}
    out_uint8 = map(x->(x[1],convert(Vector{Uint8}, x[2])), out_ascii)
    out_char = map(x->(x[1],convert(Vector{Char}, x[2])), out_ascii)

    function test_fastaread(T::Type, out)
        fr = FastaReader{T}(joinpath(dirname(Base.source_path()), "test.fasta.gz"))

        Test.@test readall(fr) == out

        rewind(fr)

        for (name, seq) in fr
            Test.@test name == out[fr.num_parsed][1]
            Test.@test seq == out[fr.num_parsed][2]
        end

        rewind(fr)
        while !eof(fr)
            name, seq = readentry(fr)
            Test.@test name == out[fr.num_parsed][1]
            Test.@test seq == out[fr.num_parsed][2]
        end

        close(fr)
    end

    test_fastaread(Vector{Uint8}, out_uint8)
    test_fastaread(Vector{Char}, out_char)
    test_fastaread(ASCIIString, out_ascii)
    test_fastaread(UTF8String, out_ascii)
end
