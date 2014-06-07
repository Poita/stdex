module stdex.io;

import std.stdio;
import std.algorithm;

struct Lines
{
    this(File file, size_t bufferSize)
    {
        m_file = file;
        m_chunks = file.byChunk(bufferSize);
        peek();
    }

    @property bool empty() const
    {
        return m_empty;
    }

    @property char[] front()
    {
        return m_front;
    }

    void popFront()
    {
        peek();
    }

    void peek()
    {
        if (m_chunks.empty)
        {
            m_empty = true;
        }
        else 
        {
            size_t i = 0;
            bool done = false;
            do
            {
                char[] chunk = cast(char[]) m_chunks.front[m_chunkOffset..$];
                size_t n = 0;
                while (n < chunk.length && chunk[n] != '\n')
                    ++n;
                
                m_frontBuffer.length = max(m_frontBuffer.length, i + n);
                m_frontBuffer[i..i+n] = chunk[0..n];
                i += n;

                if (n < chunk.length)
                {
                    ++n;
                    done = true;
                }
                m_chunkOffset += n;
                if (n == chunk.length)
                {
                    m_chunks.popFront();
                    m_chunkOffset = 0;
                    if (m_chunks.empty)
                        break;
                }
            }
            while (!done);
            m_front = m_frontBuffer[0..i];
        }
    }

    bool m_empty = false;
    char[] m_front;
    char[] m_frontBuffer;
    char[] m_buffer;
    private File m_file;
    File.ByChunk m_chunks;
    size_t m_chunkOffset = 0;
}

Lines byLineBuffered(File file, size_t bufferSize = 4096)
{
    return Lines(file, bufferSize);
}