module stdex.range;

import std.algorithm;
import std.array;
import std.functional;
import std.range;
import std.typecons;
import std.conv;

struct Indices(size_t N)
{
public:
    static assert(N != 0, "0 dimensions not allowed");

    this(size_t[] dimensions)
    {
        m_dimensions[] = dimensions[];
        m_indices[] = 0;
    }

    @property ref const(size_t[N]) front() const
    {
        return m_indices;
    }

    @property bool empty() const
    {
        return m_empty;
    }

    void popFront()
    {
        size_t i = 0;
        while(i != N)
        {
            m_indices[i]++;
            if (m_indices[i] != m_dimensions[i])
                break;
            m_indices[i] = 0;
            ++i;
        }
        m_empty = i == N;
    }

    // TODO: forward, bidirectional, random access

private:
    size_t[N] m_dimensions;
    size_t[N] m_indices = 0;
    bool m_empty = false;
}

Indices!2 indices(size_t a, size_t b)
{
    size_t[2] dimensions = void;
    dimensions[0] = a;
    dimensions[1] = b;
    return Indices!2(dimensions[]);
}

struct GroupSlices(alias _pred = "a == b", Range)
{
    static assert(isForwardRange!Range, "Range must be a forward range.");
    static assert(hasSlicing!Range, "Range must support slicing.");

public:
    alias binaryFun!_pred pred;
    alias typeof(Range.init[0..1]) Slice;

    this(Range r)
    {
        m_input = r;
        peek();
    }

    @property Slice front()
    {
        return m_front;
    }

    @property bool empty()
    {
        return m_input.empty;
    }

    void popFront()
    {
        m_input = m_next;
        if (!m_input.empty)
            peek();
    }

    @property auto save()
    {
        return typeof(this)(m_input.save);
    }

    // TODO: back, popBack

private:
    void peek()
    {
        m_next = m_input.save;
        size_t n = 1;
        m_next.popFront();
        while (!m_next.empty && pred(m_input.front, m_next.front))
        {
            m_next.popFront();
            ++n;
        }
        m_front = m_input[0..n];
    }

    Range m_input;
    Range m_next;
    Slice m_front;
}

GroupSlices!(pred, Range) groupSlices(alias pred = "a == b", Range)(Range r)
{
    return GroupSlices!(pred, Range)(r);
}

auto mapZip(alias f, Range)(Range r)
    if (isInputRange!Range && 
        is(typeof(unaryFun!f(r.front))))
{
    alias unaryFun!f fun;
    return r.map!(e => tuple(e, fun(e)));
}

unittest
{
    import std.stdio;
    writeln("testing stdex.range.groupSlices");

    int[] a = [1, 2, 3, 4, 5];
    assert(equal(a.groupSlices(), [[1], [2], [3], [4], [5]]));
    assert(equal(a.groupSlices!"a%2==b%2"(), [[1], [2], [3], [4], [5]]));
    assert(equal(a.groupSlices!"(a<3)==(b<3)"(), [[1, 2], [3, 4, 5]]));
}

struct MultiArrayIndices
{
public:
    this(Dimensions)(Dimensions divs)
    {
        size_t n = divs.length;
        m_divs = new size_t[n];
        copy(divs, m_divs);
        m_front = new size_t[n];
        m_front[] = 0;
        m_empty = false;
    }

    @property auto front()
    {
        return m_front;
    }

    @property bool empty()
    {
        return m_empty;
    }

    void popFront()
    {
        size_t i = 0;
        while(true)
        {
            m_front[i]++;
            if (m_front[i] == m_divs[i])
            {
                m_front[i] = 0;
                ++i;
                if (i == m_front.length)
                {
                    m_empty = true;
                    break;
                }
            }
            else
                break;
        }
    }

private:
    size_t[] m_divs;
    size_t[] m_front;
    bool m_empty;
}

unittest
{
    import std.stdio;
    writeln("testing stdex.range.MultiArrayIndeces");
    size_t[][] ans = [[0, 0, 0], [0, 1, 0], [0, 2, 0], [0, 0, 1], [0, 1, 1], [0, 2, 1]];
    assert(equal(MultiArrayIndices([1, 3, 2]), ans));
}

struct Generate(alias f)
{
public:
    @property auto front()
    {
        if (!m_haveFront)
            m_front = f();
        return m_front;
    }

    enum bool empty = false;

    void popFront()
    {
        m_haveFront = false;
    }

private:
    typeof(f()) m_front;
    bool m_haveFront = false;
}

auto generate(alias f)()
{
    return Generate!f();
}