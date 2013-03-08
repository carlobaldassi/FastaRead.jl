module FastaRead
export
    FastaReader,
    rewind

using GZip

import Base.start, Base.done, Base.next, Base.readall,
       Base.close


const fasta_buffer_size = 4096

type FastaReader
    f::IO
    rbuffer::Vector{Uint8} # read buffer
    rbuf_sz::Int           # read buffer size
    rbuf_pos::Int          # read buffer cursor
    is_eof::Bool           # did we reach end of file?
    lbuffer::Vector{Uint8} # line buffer
    lbuf_sz::Int           # line buffer size
    mbuffer::Vector{Uint8} # multi-line buffer
    mbuf_sz::Int           # multi-line buffer size
    function FastaReader(filename::String)
        new(gzopen(filename), Array(Uint8, fasta_buffer_size), 0, 0, false,
            Array(Uint8, fasta_buffer_size), 0,
            Array(Uint8, fasta_buffer_size), 0)
    end
end

close(fr::FastaReader) = close(fr.f)
rewind(fr::FastaReader) = seek(fr.f, 0)

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
    fr.is_eof = false
    readline(fr)
    if fr.lbuf_sz == 0
        error("empty fasta file")
    end
    return 0
end
done(fr::FastaReader, lnum) = fr.is_eof
function next(fr::FastaReader, lnum)
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
    return (name, fr.mbuffer[1:fr.mbuf_sz]), lnum + 1
end

function readall(fr::FastaReader)
    ret = Any[]
    for item in fr
        push!(ret, item)
    end
    return ret
end

end # module FastaRead
