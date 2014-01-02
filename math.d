module stdex.math;

import std.algorithm;
import std.range;

T sqr(T)(T x)
{
    return x * x;
}

T dot(T)(in T[] a, in T[] b)
{
    // TODO: ranges
    // TODO: vectorize
    T r = 0;
    foreach (i; 0..a.length)
        r += a[i] * b[i];
    return r;
}

V dot(K, V)(in V[K] a, in V[K] b)
{
    V r = 0;
    foreach (k, v; a)
    {
        const(V)* bv = k in b;
        if (bv)
            r += v * (*bv);
    }
    return r;
}

V dot(K, V)(in V[K] a, in V[] b)
{
    V r = 0;
    foreach (k, v; a)
        if (k < b.length)
            r += v * b[k];
    return r;
}

V dot(K, V)(in V[] a, in V[K] b)
{
    V r = 0;
    foreach (k, v; b)
        if (k < a.length)
            r += a[k] * v;
    return r;
}

T clamp(T)(in T x, in T a, in T b)
{
    return min(max(x, a), b);
}

auto sum(Range)(Range r, ElementType!Range seed = 0)
{
    return reduce!("a + b")(seed, r);
}