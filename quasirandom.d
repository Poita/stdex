module stdex.quasirandom;

import stdex.number;

debug
{
    import std.string;
}

/++
    Gets the nth Van der Corput number of a particular base.
  +/
double vanDerCorput(size_t n, uint base) pure
{
    double x = 0.0;
    double denum = 1.0;
    double baseFloat = base;

    while (n != 0)
    {
        denum *= baseFloat;
        x += (n % base) / denum;
        n /= base;
    }
    return x;
}

/++
    Infinite sequence of Van der Corput numbers of a given base.
    Iterating this sequence is faster than generating them individually.
  +/
struct VanDerCorput
{
public:
    this(uint base)
    {
        m_base = base;
    }

    @property double front() const pure
    {
        return vanDerCorput(m_offset, m_base);
    }

    void popFront() pure
    {
        ++m_offset;
    }

    enum bool empty = false; // Infinite range

    double opIndex(size_t i) const pure
    {
        return vanDerCorput(m_offset + i, m_base);
    }

private:
    size_t m_offset = 0;
    uint m_base = 2;
}

/++
    Writes the nth Halton vector to output.
  +/
void halton(size_t n, double[] output) pure
{
    size_t d = output.length;
    debug assert(d < smallPrimes.length, format("%d-dimensional space not supported.", d));
    foreach (i; 0..d)
    {
        uint base = smallPrimes[i];
        output[i] = vanDerCorput(n, base);
    }
}

/++
    The Halton sequence.
  +/
struct Halton(size_t Order)
{
public:
    this(size_t offset)
    {
        m_offset = offset;
    }

    @property double[Order] front() const pure
    {
        double[Order] output = void;
        halton(m_offset, output[]);
        return output;
    }

    void popFront() pure
    {
        ++m_offset;
    }

    enum bool empty = false; // Infinite range

    double[Order] opIndex(size_t i) const pure
    {
        double[Order] output = void;
        halton(m_offset + i, output[]);
        return output;
    }

private:
    size_t m_offset = 0;   
}

version(unittest)
{
    import std.algorithm;
    import std.math;
    import std.range;
    import std.stdio;
    import std.string;
}

unittest
{
    import std.stdio;
    writeln("testing stdex.quasirandom");

    immutable double[] vdc2 = [ 0.0,        1.0 / 2,    1.0 / 4,    3.0 / 4,    1.0 / 8,
                                5.0 / 8,    3.0 / 8,    7.0 / 8,    1.0 / 16,   9.0 / 16 ];
    immutable double[] vdc3 = [ 0.0,        1.0 / 3,    2.0 / 3,    1.0 / 9,    4.0 / 9,
                                7.0 / 9,    2.0 / 9,    5.0 / 9,    8.0 / 9,    1.0 / 27 ];
    immutable double[] vdc4 = [ 0.0,        1.0 / 4,    1.0 / 2,    3.0 / 4,    1.0 / 16,
                                5.0 / 16,   9.0 / 16,   13.0 / 16,  1.0 / 8,    3.0 / 8 ];
    immutable double[] vdc5 = [ 0.0,        1.0 / 5,    2.0 / 5,    3.0 / 5,    4.0 / 5,
                                1.0 / 25,   6.0 / 25,   11.0 / 25,  16.0 / 25,  21.0 / 25 ];
    immutable double[][6] vdc = [null, null, vdc2, vdc3, vdc4, vdc5];

    static bool approx(double a, double b) { return abs(a-b) < 0.0000001; }

    // Van der Corput testing
    foreach (n; 0..10)
        foreach (base; 2..6)
            assert(approx(vanDerCorput(n, base), vdc[base][n]),
                   format("vanDerCorput(%d, %d) == %f != %f", n, base, vanDerCorput(n, base), vdc[base][n]));

    foreach (base; 2..6)
        assert(VanDerCorput(base).take(10).equal!approx(vdc[base]),
               format("VanDerCorput(%d) incorrect", base));

    // Halton testing
    assert(Halton!2().take(10).map!("a[0]")().equal!approx(vdc[2]), "Halton incorrect.");
    assert(Halton!2().take(10).map!("a[1]")().equal!approx(vdc[3]), "Halton incorrect.");
}