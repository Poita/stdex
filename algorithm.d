module stdex.algorithm;

import std.algorithm;
import std.array;
import std.functional;
import std.range;
import std.traits;

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

auto longestCommonSubsequence(R1, R2)(R1 a, R2 b)
{
    // TODO: use hirschberg's algorithm
    alias T = CommonType!(Unqual!(ElementType!R1), Unqual!(ElementType!R1));
    auto L = new uint[][](a.length + 1, b.length + 1);
 
    foreach (i; 0..a.length)
    {
        import std.stdio;
        writeln(i);
        foreach (j; 0..b.length)
            L[i + 1][j + 1] = (a[i] == b[j]) ? (1 + L[i][j]) :
                              max(L[i + 1][j], L[i][j + 1]);
    }
 
    Unqual!T[] result;
    for (auto i = a.length, j = b.length; i > 0 && j > 0; ) {
        if (a[i - 1] == b[j - 1]) {
            result ~= a[i - 1];
            i--;
            j--;
        } else
            if (L[i][j - 1] < L[i - 1][j] || L[i][j - 1] == L[i - 1][j] && a[i - 1] < b[j - 1])
                i--;
            else
                j--;
    }
 
    result.reverse();
    return result;
}

///
unittest
{
    assert(longestCommonSubsequence("hello"d, "hero"d) == "heo"d);
}

bool hasSubsequence(R1, R2)(R1 a, R2 b)
{
    while (!a.empty && !b.empty)
    {
        if (a.front == b.front)
            b.popFront();
        a.popFront();
    }
    return b.empty;
}

///
unittest
{
    assert("hello".hasSubsequence("hlo"));
    assert(!"hello".hasSubsequence("hle"));
}