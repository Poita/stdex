module stdex.optimize;

import std.algorithm;
import std.range;
import std.math;

import stdex.math;
import stdex.range;

/++
    Finds the maximum of the unimodal function $(D f) within the domain between
    $(D left) and $(D right). $(D f) must strictly increase from $(D left) then
    strictly decrease towards $(D right) for ternary search to work.
  +/
T ternarySearch(alias f, T)(T left, T right, T precision)
{
    assert(left < right, "Left/right are wrong way round.");
    while (right - left > precision)
    {
        T leftThird = (2 * left + right) / 3;
        T rightThird = (left + 2 * right) / 3;
        if (f(leftThird) < f(rightThird))
        {
            left = leftThird;
        }
        else
        {
            right = rightThird;
        }
    }
    return (left + right) / 2;
}

unittest
{
    import std.stdio;
    writeln("testing stdex.optimize.ternarySearch");
    assert(abs(ternarySearch!(x => -(1 - x) * (1 - x))(-10.0f, 10.0f, 0.0001f) - 1.0f) < 0.001f);
}

/++
   Finds the minimum of function $(D f) by sampling between $(D lower) and
   $(D upper) with the $(I i)th dimension subdivided $(D divs[i]) times.
   The sample where $(D f) is minimum is returned.
  +/
ElementType!V1[] gridSearch(alias f, V1, V2, V3)(V1 lower, V2 upper, V3 divs)
{
    alias ElementType!V1 E;

    size_t n = walkLength(lower.save);
    E[] working = new E[n];
    E[] best = new E[n];
    copy(lower, best);
    auto bestValue = f(best);

    foreach (size_t[] idx; MultiArrayIndices(divs))
    {
        foreach (i; 0..n)
            working[i] = cast(E)(lower[i] + idx[i] * (upper[i] - lower[i]) / (divs[i] - 1));

        auto workingValue = f(working);
        if (workingValue < bestValue)
        {
            bestValue = workingValue;
            copy(working, best);
        }
    }

    return best;
}

unittest
{
    import std.stdio;
    writeln("testing stdex.optimize.gridSearch");

    static float f0(float[] v)
    {
        return (v[0] * v[0]) + (v[1] * v[1]);
    }

    float[] a0 = gridSearch!f0([-10.0f, -10.0f], [20.0f, 20.0f], [31, 31]);
    assert(a0 == [0.0f, 0.0f]);
}
