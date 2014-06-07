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
    if (b.length < a.length)
        return dot(b, a);
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
        if (k < b.length && k >= 0)
            r += v * b[k];
    return r;
}

V dot(K, V)(in V[] a, in V[K] b)
{
    V r = 0;
    foreach (k, v; b)
        if (k < a.length && k >= 0)
            r += a[k] * v;
    return r;
}

T clamp(T)(in T x, in T a, in T b)
{
    return min(max(x, a), b);
}

T clamp01(T)(in T x)
{
    return min(max(x, 0), 1);
}

U lerp(T, U)(in T t, in U from, in U to)
{
    return (1 - t) * from + t * to;
}

U lerpClamped(T, U)(in T t, in U from, in U to)
{
    return lerp(clamp01(t), from, to);
}

T inverseLerp(T)(in T x, in T from, in T to)
{
    assert(to != from, "Cannot inverse lerp between equal values.");
    return (x - from) / (to - from);
}

T inverseLerpClamp(T)(in T x, in T from, in T to)
{
    return inverseLerp(x, from, to).clamp01();
}

auto sum(Range)(Range r, ElementType!Range seed = 0)
{
    return reduce!("a + b")(seed, r);
}