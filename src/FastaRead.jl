module FastaRead
export
    FastaReader,
    readentry,
    rewind

using GZip

import Base.start, Base.done, Base.next, Base.readall,
       Base.close, Base.show, Base.eof

const fasta_buffer_size = 4096

type FastaReader{T}
    # public read-only
    filename::String
    num_parsed::Int        # number of parsed entries so far
    # private
    f::IO
    is_eof::Bool           # did we reach end of file?
    rbuffer::Vector{Uint8} # read buffer
    rbuf_sz::Int           # read buffer size
    rbuf_pos::Int          # read buffer cursor
    lbuffer::Vector{Uint8} # line buffer
    lbuf_sz::Int           # line buffer size
    mbuffer::Vector{Uint8} # multi-line buffer
    mbuf_sz::Int           # multi-line buffer size
    function FastaReader(filename::String)
        new(filename, 0, gzopen(filename), false,
            Array(Uint8, fasta_buffer_size), 0, 0,
            Array(Uint8, fasta_buffer_size), 0,
            Array(Uint8, fasta_buffer_size), 0)
    end
end

FastaReader(filename::String) = FastaReader{ASCIIString}(filename)

close(fr::FastaReader) = close(fr.f)
rewind(fr::FastaReader) = (seek(fr.f, 0); fr.is_eof = false; fr.num_parsed = 0; nothing)

function read_chunk(fr::FastaReader)
    if fr.is_eof
        return
    end
    ret = gzread(fr.f, pointer(fr.rbuffer), fasta_buffer_size)
    if ret == -1
        error("gzread failed")
    end
    fr.rbuf_sz = ret
    fr.rbuf_pos = 1
    if ret == 0
        fr.is_eof = true
    end
    return
end

function readline(fr::FastaReader)
    fr.lbuf_sz = 0
    found = false
    while !fr.is_eof
        if fr.rbuf_pos == 0
            read_chunk(fr::FastaReader)
        end
        i = fr.rbuf_pos
        while i <= fr.rbuf_sz
            if fr.rbuffer[i] == '\n'
                found = true
                break
            end
            i += 1
        end
        i -= 1
        chunk_len = i - fr.rbuf_pos + 1
        free_sbuf = length(fr.lbuffer) - fr.lbuf_sz
        gap = chunk_len - free_sbuf
        if gap > 0
            resize!(fr.lbuffer, length(fr.lbuffer) + gap)
        end

        #fr.lbuffer[fr.lbuf_sz + (1:chunk_len)] = fr.rbuffer[fr.rbuf_pos:i]
        copy!(fr.lbuffer, fr.lbuf_sz + 1, fr.rbuffer, fr.rbuf_pos, chunk_len)
        fr.lbuf_sz += chunk_len

        i += 2
        if i > fr.rbuf_sz
            i = 0
        end
        fr.rbuf_pos = i
        if found
            break
        end
    end
    return
end

function start(fr::FastaReader)
    rewind(fr)
    readline(fr)
    if fr.lbuf_sz == 0
        error("empty fasta file")
    end
    return
end
done(fr::FastaReader, x::Nothing) = fr.is_eof
function _next_step(fr::FastaReader)
    if fr.lbuffer[1] != '>'
        error("invalid fasta file: seq name does not start with '>'")
    end
    if length(fr.lbuffer) == 1
        error("invalid fasta file: empty seq name")
    end
    name = ascii(fr.lbuffer[2:fr.lbuf_sz])
    fr.mbuf_sz = 0
    while true
        readline(fr)
        if fr.is_eof || fr.lbuffer[1] == '>'
            break
        end
        gap = fr.lbuf_sz - (length(fr.mbuffer) - fr.mbuf_sz)
        if gap > 0
            resize!(fr.mbuffer, length(fr.mbuffer) + gap)
        end
        #fr.mbuffer[fr.mbuf_sz + (1:fr.lbuf_sz)] = fr.lbuffer[1:fr.lbuf_sz]
        copy!(fr.mbuffer, fr.mbuf_sz + 1, fr.lbuffer, 1, fr.lbuf_sz)
        fr.mbuf_sz += fr.lbuf_sz
    end
    return name
end
function _next(fr::FastaReader{Vector{Uint8}})
    name = _next_step(fr)
    fr.num_parsed += 1
    return (name, fr.mbuffer[1:fr.mbuf_sz])
end
function _next(fr::FastaReader{ASCIIString})
    name = _next_step(fr)
    out_str = ccall(:jl_pchar_to_string, ByteString, (Ptr{Uint8},Int), fr.mbuffer, fr.mbuf_sz)
    fr.num_parsed += 1
    return (name, out_str)
end
function _next{T}(fr::FastaReader{T})
    name = _next_step(fr)
    fr.num_parsed += 1
    return (name, convert(T, fr.mbuffer[1:fr.mbuf_sz]))
end

next(fr::FastaReader, x::Nothing) = (_next(fr), nothing)

function readall(fr::FastaReader)
    ret = Any[]
    for item in fr
        push!(ret, item)
    end
    return ret
end

function readentry(fr::FastaReader)
    fr.is_eof && throw(EOFError())
    if fr.num_parsed == 0
        readline(fr)
        if fr.lbuf_sz == 0
            error("empty fasta file")
        end
    end
    item, _ = next(fr, nothing)
    return item
end

eof(fr::FastaReader) = fr.is_eof

function show{T}(io::IO, fr::FastaReader{T})
    print(io, "FastaReader(filename=\"$(fr.filename)\", out_type=$T, eof=$(fr.is_eof))")
end

end # module FastaRead
