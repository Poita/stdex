/++
    Module containing utility functions that are of general use in all
    programs.

    This module provides:
    $(UL
        $(LI Formatted assert convenience function.)
    )

    Note:
        None

    See_Also:
        Nothing

    Copyright: Copyright 2011-2012
    License:   $(WEB www.boost.org/LICENSE_1_0.txt, Boost License 1.0).
    Authors:   Peter J Alexander
  +/

module stdex.util;

import std.typecons;
import std.typetuple;

/++
    Throws an AssertError if test is false. The args can be used to provide a
    format string.
  +/
void assertf(T, string file = __FILE__, int line = __LINE__, Args...)
            (lazy T test, lazy Args args)
{
    //version(assert)
    {
        import std.string : format;
        import core.exception : AssertError;
        if (!test)
        {
            throw new AssertError(format(args), file, line);
        }
    }
}

/++
    Produces an array of key-value pair tuples from an associative array.
    The order of items produces is the same order produced by foreach
    iteration, which is unspecified.
  +/
Tuple!(Key, Value)[] keyValueArray(Key, Value)(Value[Key] aa)
{
    alias Tuple!(Key, Value) Element;
    Element[] arr = new Element[aa.length];
    size_t i = 0;
    foreach (key, value; aa)
        arr[i++] = Element(key, value);
    return arr;
}

/++
    Interleaves two type tuples.
  +/
template Interleave(A...)
{
    template and(B...)
        if (A.length == B.length)
    {
        static if (A.length == 1)
        {
            alias TypeTuple!(A[0], B[0]) and;
        }
        else
        {
            alias TypeTuple!(A[0], B[0], Interleave!(A[1..$]).and!(B[1..$])) and;
        }
    }
}

/++
    Constructs a tuple with named members.

    namedTuple!("a", "b")(a, b) is equivalent to
    Tuple!(typeof(a), "a", typeof(b), "b")(a, b).
  +/
template namedTuple(names...)
{
    auto namedTuple(Args...)(Args args)
        if (Args.length == names.length &&
            args.length > 0)
    {
        return Tuple!(Interleave!(Args).and!(names))(args);
    }
}