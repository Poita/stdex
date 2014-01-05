module stdex.algorithm;

import std.algorithm;
import std.array;
import std.functional;
import std.range;

void forEach(alias func, Range)(Range range)
{
	alias unaryFun!func f;
	foreach (e; range)
		f(e);
}

auto sorted(alias pred = "a < b", InputRange)(InputRange range)
{
    static assert(isInputRange!InputRange, "input must be an input range.");
    static assert(!isInfinite!InputRange, "input must not be an infinite range.");

    return range.array.sort!pred();
}

auto argMin(alias f, alias less = "a < b", InputRange)(InputRange args)
{
    static assert(isInputRange!InputRange, "args must be an input range.");
    static assert(!isInfinite!InputRange, "args must not be an infinite range.");

    assert(!args.empty, "Cannot find argmin of empty range.");
    auto best = args.front;
    auto bestVal = unaryFun!f(best);
    args.popFront();
    while (!args.empty)
    {
        auto testVal = unaryFun!f(args.front);
        if (binaryFun!less(testVal, bestVal))
        {
            bestVal = testVal;
            best = args.front;
        }
        args.popFront();
    }
    return best;
}

auto argMax(alias f, alias more = "a > b", InputRange)(InputRange args)
{
    return argMin!(f, more, InputRange)(args);
}

ulong[ElementType!Range] histogram(Range)(Range input)
{
    static assert(isInputRange!Range, "input must be an input range.");
    static assert(!isInfinite!Range, "input must not be an infinite range.");

    ulong[ElementType!Range] histo;
    foreach (e; input)
        histo[e]++;
    return histo;
}

unittest
{
    import std.stdio;
    writeln("testing stdex.algorithm.histogram");

    auto h1 = [1, 1, 2, 2, 2, 3, 4, 5, 5].histogram();
    assert(h1[1] == 2);
    assert(h1[2] == 3);
    assert(h1[3] == 1);
    assert(h1[4] == 1);
    assert(h1[5] == 2);
}

ulong fnv1(Range)(Range x)
{
    ulong h = 14695981039346656037UL;
    foreach (b; x)
        h = (h * 1099511628211UL) ^ b;
    return h;
}
