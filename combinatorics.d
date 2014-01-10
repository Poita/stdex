// Written in the D programming language.

/**
This module defines combinatorics stuff.

Source: $(PHOBOSSRC std/_combinatorics.d)

Macros:

WIKI = Phobos/StdCombinatorics

Copyright: Copyright by authors 2012-.

License: $(WEB boost.org/LICENSE_1_0.txt, Boost License 1.0).

Authors: $(WEB poita.org, Peter Alexander).
 */
module stdex.combinatorics;

import core.bitop;
import std.algorithm;
import std.array;
import std.conv;
import std.exception;
import std.functional;
import std.range;
import std.traits;

version (unittest)
{
    import std.bigint;
    import std.conv;
    import std.stdio;
    import std.string;
    import std.typetuple;
}

/**
   Computes the factorial of a non-negative integer.

   The factorial of $(D n) is defined to be the product of all integers
   between 1 and $(D n) (inclusive). The factorial of 0 is 1.
 */
Int factorial(Int = ulong)(uint n)
{
    // TODO: use lookup tables
    // TODO: handle floats
    if (n < 2)
        return 1;
    Int result = cast(Int) n;
    while (--n != 1)
        result *= n; // TODO: handle overflow
    return result;
}

unittest
{
    foreach (T; TypeTuple!(short, ushort, int, uint, long, ulong))
    {
        assert(equal(iota(8).map!(factorial!T)(), [1, 1, 2, 6, 24, 120, 720, 5040]));
    }

    assert(factorial!ulong(20) == 2_432_902_008_176_640_000UL);
    // TODO: test bigint
}

/**
   Computes the multinomial coefficient of a set of integers.
 */
Int multinomial(Int, Range)(Range ks)
{
    // TODO: use better algo to avoid intermediate overflow
    // http://home.comcast.net/~tamivox/dave/multinomial/index.html
    alias reduce!((a, b) => a + b) sum;
    alias reduce!((a, b) => cast(Int)(a * b)) product;
    alias factorial!Int fact;

    uint n = sum(0u, ks);
    return to!Int(fact(n) / product(cast(Int) 1, map!fact(ks)));
}

/**
   Computes the binomial coefficient.
 */
Int binomial(Int)(uint n, uint k)
{
    // TODO: see links
    if (k > n)
        return 0;
    k = min(k, n - k);
    Int c = 1;
    foreach (uint i; 1..k+1)
    {
        // TODO: possible intermediate overflow
        c *= (n - (k - i));
        c /= i;
    }
    return c;
}

unittest
{
    int[][] pascal = [
        [1],
        [1,  1],
        [1,  2,  1],
        [1,  3,  3,  1],
        [1,  4,  6,  4,  1],
        [1,  5,  10, 10, 5,  1],
        [1,  6,  15, 20, 15, 6,  1]
    ];
    foreach (T; TypeTuple!(short, ushort, int, uint, long, ulong))
    {
        alias binomial!T bc;
        foreach (n; 0..pascal.length)
            foreach (k; 0..pascal[n].length)
                assert(bc(cast(uint)n, cast(uint)k) == pascal[n][k]);
        assert(bc(26, 3) == 2600);
        //assert(bc(31, 4) == 31465);
        // TODO: this overflows, need overflow-free algo
    }
}

/**
   Gets the $(D n)th Bell number.
 */
Int bellNumber(Int)(uint n)
{
    // TODO: overflow
    if (n < bellNumbersCache.length)
        return bellNumbersCache[n].to!Int();
    // TODO: BigInt algorithm
    // TODO: real number algorithm
    return 0;
}

// Cache of all Bell numbers below ulong.max
private immutable ulong[] bellNumbersCache = [ 1UL,
                                               1UL,
                                               2UL,
                                               5UL,
                                              15UL,
                                              52UL,
                                             203UL,
                                             877UL,
                                           4_140UL,
                                          21_147UL,
                                         115_975UL,
                                         678_570UL,
                                       4_213_597UL,
                                      27_644_437UL,
                                     190_899_322UL,
                                   1_382_958_545UL,
                                  10_480_142_147UL,
                                  82_864_869_804UL,
                                 682_076_806_159UL,
                               5_832_742_205_057UL,
                              51_724_158_235_372UL,
                             474_869_816_156_751UL,
                           4_506_715_738_447_323UL,
                          44_152_005_855_084_346UL,
                         445_958_869_294_805_289UL,
                       4_638_590_332_229_999_353UL ];

/**
   Gets the $(D n)th Catalan number.
 */
Int catalanNumber(Int)(uint n)
{
    // TODO: overflow
    if (n < catalanNumberCache.length)
        return catalanNumberCache[n].to!Int();
    // TODO: BigInt algorithm
    // TODO: real number algorithm
    return 0;
}

// Cache of all the Catalan numbers below ulong.max
private immutable ulong[] catalanNumberCache = [ 1UL,
                                                 1UL,
                                                 2UL,
                                                 5UL,
                                                14UL,
                                                42UL,
                                               132UL,
                                               429UL,
                                             1_430UL,
                                             4_862UL,
                                            16_796UL,
                                            58_786UL,
                                           208_012UL,
                                           742_900UL,
                                         2_674_440UL,
                                         9_694_845UL,
                                        35_357_670UL,
                                       129_644_790UL,
                                       477_638_700UL,
                                     1_767_263_190UL,
                                     6_564_120_420UL,
                                    24_466_267_020UL,
                                    91_482_563_640UL,
                                   343_059_613_650UL,
                                 1_289_904_147_324UL,
                                 4_861_946_401_452UL,
                                18_367_353_072_152UL,
                                69_533_550_916_004UL,
                               263_747_951_750_360UL,
                             1_002_242_216_651_368UL,
                             3_814_986_502_092_304UL,
                            14_544_636_039_226_909UL ];

/**
   Next permutation.
 */
bool nextPermutation(alias pred = ((a, b) => a < b), Range)(Range r)
    if (isBidirectionalRange!Range &&
        hasSwappableElements!Range &&
        is(typeof(binaryFun!pred(r.front, r.front)) : bool))
{
    // TODO: optimize for array
    if (!r.empty)
    {
        auto iter1 = r.save;
        auto iter2 = r.save;
        iter1.popBack();
        for (size_t n = 1; !iter1.empty; ++n)
        { 
            alias binaryFun!pred comp;
            if (comp(iter1.back, iter2.back))
            {
                auto iter3 = r.save;
                while (!comp(iter1.back, iter3.back))
                    iter3.popBack();
                swap(iter1.back, iter3.back);
                static if (hasSlicing!Range && hasLength!Range)
                    reverse(r[r.length - n .. r.length]);
                else
                    reverse(retro(r).takeExactly(n));
                return true;
            }
            iter1.popBack();
            iter2.popBack();
        }
        reverse(r);
    }
    return false;
} // TODO: add tests

/**
Previous permutation.
 */
bool previousPermutation(alias pred = ((a, b) => a < b), Range)(Range r)
    if (isBidirectionalRange!Range &&
        hasSwappableElements!Range &&
        is(typeof(binaryFun!pred(r.front, r.front)) : bool))
{
    alias binaryFun!pred comp;
    return nextPermutation!((a, b) => comp(b, a), Range)(r);
} // TODO: add tests

/**
    Count the number of permutations of a range under an ordering predicate.
    Multisets are handled correctly.
 */
Int countPermutations(Int = ulong, alias pred = ((a, b) => a < b), Range)(Range r)
{
    alias binaryFun!pred order;
    alias binaryFun!((a, b) => !(order(a, b) || order(b, a))) equivalent;
    alias unaryFun!(a => a[1]) groupCount;

    // TODO: check order
    auto set = array(r);
    sort!order(set);
    return multinomial!Int(set.group!equivalent().map!groupCount());
}

unittest
{
    assert(countPermutations(cast(int[]) []) == 1);
    assert(countPermutations([1]) == 1);
    assert(countPermutations([1, 2]) == 2);
    assert(countPermutations([1, 2, 3]) == 6);
    assert(countPermutations([3, 1, 2]) == 6);
    assert(countPermutations([1, 1, 1]) == 1);
    assert(countPermutations([1, 2, 1]) == 3);
    assert(countPermutations([2, 1, 2, 1]) == 6);
    // TODO: larger test
    // TODO: test more range types
}

/**
   Returns a range of the permutations of $(D r).
 */
auto permutations(alias pred = ((a, b) => a < b), Range)(Range r)
    if (isInputRange!Range &&
        is(typeof(binaryFun!pred(r.front, r.front)) : bool))
{
    return Permutations!(pred, Unqual!(ElementType!Range))(r);
}

private struct Permutations(alias pred = ((a, b) => a < b), Value)
{
    alias binaryFun!pred comp;

    this(Range)(Range source)
    {
        _front = source.save.array();
        _back = _front.dup;
        _length = countPermutations!(size_t, pred)(_front);
        previousPermutation(_back);
        _onLast = equal(_front, _back);
        _done = false;
    }

    @property auto front() { return _front.dup; }
    @property auto back() { return _back; }
    @property bool empty() const { return _done; }
    @property size_t length() const { return _length; }

    void popFront()
    {
        assert(!_done);
        nextPermutation!pred(_front);
        _done = _onLast;
        _onLast = equal(_front, _back);
        --_length;
    }

    void popBack()
    {
        assert(!_done);
        previousPermutation!pred(_back);
        _done = _onLast;
        _onLast = equal(_front, _back);
        --_length;
    }

    @property auto save()
    {
        Permutations copy = void;
        copy._front = _front.dup;
        copy._back = _back.dup;
        copy._length = _length;
        copy._onLast = _onLast;
        copy._done = _done;
        return copy;
    }

    // TODO: random access
    // TODO: different orders

    private Value[] _front;
    private Value[] _back;
    private size_t _length;
    private bool _onLast = false;
    private bool _done = false;
}

unittest
{
    alias filter!"true" F;
    int[][] perms123 = [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]];
    
    assert(equal!equal(permutations([1, 2, 3]), perms123));
    assert(equal!equal(permutations([1, 2, 3]).retro(), perms123.retro()));
    assert(equal!equal(permutations(chain([1, 2], [3])), perms123));
    assert(equal!equal(permutations(chain([1, 2], [3])).retro(), perms123.retro()));
    assert(permutations(iota(1)).filter!"true"().walkLength() == 1);
    assert(permutations(iota(2)).filter!"true"().walkLength() == 2);
    assert(permutations(iota(3)).filter!"true"().walkLength() == 6);
    assert(permutations(iota(4)).filter!"true"().walkLength() == 24);
    assert(permutations(iota(5)).filter!"true"().walkLength() == 120);
    assert(permutations(iota(6)).filter!"true"().walkLength() == 720);

    int[][] perms112 = [[1, 1, 2], [1, 2, 1], [2, 1, 1]];
    assert(equal!equal(permutations([1, 1, 2]), perms112));
    assert(equal!equal(permutations([1, 1, 2]).retro(), perms112.retro()));

    int[][] perms111 = [[1, 1, 1]];
    assert(equal!equal(permutations([1, 1, 1]), perms111));
    assert(equal!equal(permutations([1, 1, 1]).retro(), perms111));

    int[] emptyArray = [];
    assert(equal!equal(permutations(emptyArray), [emptyArray]));

    // TODO: test more range types
    // TODO: try larger permutations
}

private struct Subset(Range)
    if (isInputRange!Range)
{
    // TODO: add bidir support
    // TODO: add random access support ???
    this(Range source, ulong mask) // TODO: support larger subsets
    {
        _current = source;
        _mask = mask;
        _length = popcnt(cast(uint) mask) + popcnt(cast(uint)(mask >> 32));
        static if (hasLength!Range)
            assert(_length <= _current.length, "Mask invalid for range.");
        moveNext();
    }

    @property auto front() { return _current.front; }
    @property bool empty() const { return _length == 0; }
    @property size_t length() const { return _length; }

    void popFront()
    {
        _mask >>= 1;
        _current.popFront();
        --_length;
        moveNext();
    }

    static if (isForwardRange!Range)
        @property Subset save()
        {
            Subset copy = void;
            copy._current = _current.save;
            copy._mask = _mask;
            copy._length = _length;
            return copy;
            // TODO: tests
        }

    private void moveNext()
    {
        while ((_mask & 1) == 0 && _mask != 0 && !_current.empty)
        {
            // TODO: count trailing 0's and use popFrontN
            _mask >>= 1;
            _current.popFront();
        }
        assert(!(_current.empty && _mask != 0), "Mask invalid for range.");
    }

    private Range _current;
    private ulong _mask = 0;
    private size_t _length = 0;
}

/**
   Returns the subset of $(D range) defined by $(D mask). The $(D n)th element
   of $(D range) will be present in the subset if the ($ n)th bit is set in
   $(D mask).
  */
auto subset(Range)(Range range, ulong mask)
    if (isInputRange!Range)
{
    return Subset!Range(range, mask);
}

unittest 
{
    dstring a = "abc"d;
    assert(equal(subset(a, 0), ""));
    assert(equal(subset(a, 1), "a"));
    assert(equal(subset(a, 2), "b"));
    assert(equal(subset(a, 3), "ab"));
    assert(equal(subset(a, 4), "c"));
    assert(equal(subset(a, 5), "ac"));
    assert(equal(subset(a, 6), "bc"));
    assert(equal(subset(a, 7), "abc"));

    assert(equal(subset(iota(64), 0xFFFF_FFFF_FFFF_FFFF), iota(64)));
    assert(equal(subset(iota(64), 0x0000_0000_FFFF_FFFF), iota(32)));
    assert(equal(subset(iota(64), 0xFFFF_FFFF_0000_0000), iota(32, 64)));
    assert(equal(subset(iota(64), 0x0000_FFFF_FFFF_0000), iota(16, 48)));
    assert(equal(subset(iota(64), 0x5555_5555_5555_5555), iota(0, 64, 2)));
    assert(equal(subset(iota(64), 0xAAAA_AAAA_AAAA_AAAA), iota(1, 64, 2)));
    assert(equal(subset(iota(64), (1UL << 17) | (1UL << 63)), [17, 63]));

    // TODO: test more range types
}

private struct PowerSet(Range)
    if (isForwardRange!Range)
{
    // TODO: different orderings (colex etc.)

    enum infinite = isInfinite!Range;

    this(Range source)
    {
        _source = source;
        static if (!infinite)
        {
            auto numElements = walkLength(_source.save);
            enum maxElements = 8 * size_t.sizeof - 1;
            assert(numElements <= maxElements, "Source range is too large for power set.");
            _length = (cast(size_t) 1) << numElements;
        }
    }

    @property Subset!Range front() { return subset(_source.save, _mask); }
    @property PowerSet!Range save() { return this; }
    Subset!Range opIndex(size_t i) { return subset(_source.save, _mask + i); }

    static if (infinite)
    {
        enum bool empty = false;
        void popFront() { ++_mask; }
    }
    else
    {
        @property bool empty() const { return _length == 0; }
        void popFront() { ++_mask; --_length; }
        @property auto back() { return subset(_source.save, _mask + _length - 1); }
        void popBack() { --_length; }
        @property size_t length() const { return _length; }

        PowerSet!Range opSlice(size_t from, size_t to)
        {
            assert(from >= 0, "Slice out of bounds.");
            assert(to - from <= _length, "Slice out of bounds.");

            PowerSet slice = void;
            slice._source = _source;
            slice._length = to - from;
            slice._mask = _mask + from;
            return slice;
        }
    }

    private Range _source;
    static if (!infinite)
        private size_t _length = 0;
    private size_t _mask = 0; // TODO: support > 63 elements
}

/**
   Returns a range of all the subsets of $(D r), i.e. the $(I power set) of $(D r).
 */
PowerSet!Range subsets(Range)(Range r)
    if (isForwardRange!Range)
{
    return PowerSet!Range(r);
}

unittest
{
    auto p123 = [[], [1], [2], [1, 2], [3], [1, 3], [2, 3], [1, 2, 3]];
    auto p = subsets([1, 2, 3]);
    assert(equal!equal(p, p123));
    assert(equal!equal(retro(p), retro(p123)));
    assert(equal!equal(radial(p), radial(p123)));
    assert(equal!equal(filter!"true"(p), p123));
    assert(p.length == 8);
    assert(equal(p[0], p123[0]));
    assert(equal(p[5], p123[5]));
    assert(equal(p[7], p123[7]));
    assert(equal!equal(p[0..8], p123));
    assert(equal!equal(p[1..6], p123[1..6]));
    assert(equal!equal(p[1..5][2..4], p123[3..5]));
    assert(equal!equal(subsets(cycle([1, 2, 3])).take(8), p123));
    assert(equal!equal(subsets(iota(1, 4)), p123));

    auto p111 = [[], [1], [1], [1, 1], [1], [1, 1], [1, 1], [1, 1, 1]];
    assert(equal!equal(subsets([1, 1, 1]), p111));

    int[] emptySet = [];
    assert(equal!equal(subsets(emptySet), [emptySet]));
    assert(equal!equal(subsets([1]), [[], [1]]));
}

private struct KSubsets(Range)
    if (isForwardRange!Range)
{
    // TODO: orderings
    // TODO: length
    // TODO: random access??
    // TODO: bidirectional?
    enum infinite = isInfinite!Range;

    this(Range range, size_t k)
    {
        _source = range.save;
        _mask = ((cast(size_t) 1) << k) - 1; // lowest k bits set
        static if (!infinite)
        {
            size_t n = walkLength(range);
            _length = binomial!size_t(cast(uint)n, cast(uint)k);
        }
    }

    @property auto front() { return subset(_source.save, _mask); }

    void popFront()
    {
        // See: http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
        // TODO: this can be faster if processor has ctz support. See link.
        static if (!infinite)
        {
            if (--_length == 0)
                return;
        }
        else if (_mask == 0)
            return;

        size_t v = _mask;
        size_t t = (v | (v - 1)) + 1;  
        _mask = t | ((((t & -t) / (v & -v)) >> 1) - 1); 
    }

    static if (infinite)
    {
        enum bool empty = false;
    }
    else
    {
        @property bool empty() const { return _length == 0; }
        @property size_t length() const { return _length; }
    }

    Range _source;
    size_t _mask;
    static if (!infinite)
    {
        size_t _length;
    }
}

/**
   Returns a range of the $(D k)-length subsets of $(D range).
 */
KSubsets!Range kSubsets(Range)(Range range, size_t k)
{
    return KSubsets!Range(range, k);
}

unittest
{
    int[] emptySet = [];
    int[][] ks0_0 = [[]];
    int[][] ks1_1 = [[1]];
    int[][] ks12_1 = [[1], [2]];
    int[][] ks12_2 = [[1, 2]];
    int[][] ks123_1 = [[1], [2], [3]];
    int[][] ks123_2 = [[1, 2], [1, 3], [2, 3]];
    int[][] ks123_3 = [[1, 2, 3]];
    int[][] ks1234_1 = [[1], [2], [3], [4]];
    int[][] ks1234_2 = [[1, 2], [1, 3], [2, 3], [1, 4], [2, 4], [3, 4]];
    int[][] ks1234_3 = [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]];
    int[][] ks1234_4 = [[1, 2, 3, 4]];

    assert(equal!equal(kSubsets(emptySet, 0), ks0_0));
    assert(equal!equal(kSubsets([1], 1), ks1_1));
    assert(equal!equal(kSubsets([1, 2], 1), ks12_1));
    assert(equal!equal(kSubsets([1, 2], 2), ks12_2));
    assert(equal!equal(kSubsets([1, 2, 3], 1), ks123_1));
    assert(equal!equal(kSubsets([1, 2, 3], 2), ks123_2));
    assert(equal!equal(kSubsets([1, 2, 3], 3), ks123_3));
    assert(equal!equal(kSubsets([1, 2, 3, 4], 1), ks1234_1));
    assert(equal!equal(kSubsets([1, 2, 3, 4], 2), ks1234_2));
    assert(equal!equal(kSubsets([1, 2, 3, 4], 3), ks1234_3));
    assert(equal!equal(kSubsets([1, 2, 3, 4], 4), ks1234_4));

    assert(equal!equal(kSubsets(cycle([1, 2, 3, 4]), 0).take(1), ks0_0));
    assert(equal!equal(kSubsets(cycle([1, 2, 3, 4]), 1).take(4), ks1234_1));
    assert(equal!equal(kSubsets(cycle([1, 2, 3, 4]), 2).take(6), ks1234_2));
    assert(equal!equal(kSubsets(cycle([1, 2, 3, 4]), 3).take(4), ks1234_3));
    assert(equal!equal(kSubsets(cycle([1, 2, 3, 4]), 4).take(1), ks1234_4));

    dstring alphabet = "abcdefghijklmnopqrstuvwxyz"d;
    assert(kSubsets(alphabet, 3).length == 2600);
    assert(equal!equal(kSubsets(alphabet, 3).drop(2597), ["vyz"d, "wyz"d, "xyz"d]));

    //TODO: more range types
}

/**
   Return the number of integer partitions of $(D n).
 */
Int countIntegerPartitions(Int)(size_t n)
{
    // TODO: what to do about large n?
    
    // Lookup table.
    static Int[] lut = [1];
    if (n < lut.length && lut[n] != 0)
        return lut[n];

    Int s = 0;
    long j = cast(long)(n - 1);
    uint k = 2;
    while (j >= 0)
    {
        Int t = countIntegerPartitions!Int(cast(size_t) j);
        if ((k / 2) % 2 == 0)
            s -= t;
        else
            s += t;
        j -= (k % 2 == 0) ? k / 2 : k;
        ++k;
    }
    lut.length = max(lut.length, n + 1);
    lut[n] = s;
    return s;
}

unittest
{
    static immutable int[] first50 = [
        1,      1,      2,      3,      5,      7,      11,     15,
        22,     30,     42,     56,     77,     101,    135,    176,
        231,    297,    385,    490,    627,    792,    1002,   1255,
        1575,   1958,   2436,   3010,   3718,   4565,   5604,   6842,
        8349,   10143,  12310,  14883,  17977,  21637,  26015,  31185,
        37338,  44583,  53174,  63261,  75175,  89134,  105558, 124754,
        147273, 173525
    ]; 

    foreach (T; TypeTuple!(int, uint, long, ulong, size_t, BigInt))
        foreach (i; 0..50)
            assert(countIntegerPartitions!T(i) == first50[i]);

    assert(countIntegerPartitions!size_t(100) == 190_569_292);
    assert(countIntegerPartitions!BigInt(1000) == BigInt("24061467864032622473692149727991"));
}

struct IntegerPartitions(Int)
    if (isIntegral!Int)
{
    this(Int n)
    {
        _size = 1;
        _lastNon1 = 0;
        _buffer.length = n == 0 ? 1 : cast(size_t) n;
        _buffer[0] = n;
        _length = countIntegerPartitions!size_t(n);
    }

    @property const(Int)[] front() const { return _buffer[0.._size]; }
    @property bool empty() const { return _size == 0; }

    void popFront()
    {
        --_length;
        if (_buffer[0] == 1 || _buffer[0] == 0)
        {
            _size = 0;
            return;
        }

        assert(_size != 0);

        size_t i = _lastNon1;
        assert(_buffer[i] != 1);
        if (--_buffer[i] == 1)
            --_lastNon1;
        Int rem = cast(Int)(_size - i);
        ++i;
        for (; rem != 0; rem -= _buffer[i], ++i)
        {
            _buffer[i] = min(_buffer[i-1], rem);
            if (_buffer[i] != 1)
                _lastNon1 = i;
        }
        _size = i;
    }

    @property IntegerPartitions save() const
    {
        IntegerPartitions copy;
        copy._buffer = _buffer.dup;
        copy._size = _size;
        copy._lastNon1 = _lastNon1;
        return copy;
    }

    @property size_t length() const { return _length; }

    // TODO: bidir
    // TODO: random
    // TODO: orders

    private Int[] _buffer;
    private size_t _size;
    private size_t _lastNon1;
    private size_t _length;
}

/**
   Returns a range of integer partitions of $(D n).
 */
auto integerPartitions(Int)(Int n)
{
    return IntegerPartitions!Int(n);
}

unittest
{
    int[][] p0 = [[0]];
    int[][] p1 = [[1]];
    int[][] p2 = [[2], [1, 1]];
    int[][] p3 = [[3], [2, 1], [1, 1, 1]];
    int[][] p4 = [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]];
    int[][] p5 = [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]];
    assert(equal!equal(integerPartitions(0), p0));
    assert(equal!equal(integerPartitions(1), p1));
    assert(equal!equal(integerPartitions(2), p2));
    assert(equal!equal(integerPartitions(3), p3));
    assert(equal!equal(integerPartitions(4), p4));
    assert(equal!equal(integerPartitions(5), p5));
}

/**
   Return the number of $(D k)-length integer partitions of $(D n).
 */
Int countIntegerPartitions(Int)(size_t n, size_t k)
{
    // TODO: what to do with large values?
    if (n == 0) return 1;
    if (n < k) return 0;
    if (k == 1 || k == n) return 1;

    alias countIntegerPartitions p;
    return p(n - 1, k - 1) + p(n - k, k);
}

unittest
{
    assert(countIntegerPartitions!size_t(0, 1) == 1);
    assert(countIntegerPartitions!size_t(1, 1) == 1);
    assert(countIntegerPartitions!size_t(2, 1) == 1);
    assert(countIntegerPartitions!size_t(2, 2) == 1);
    assert(countIntegerPartitions!size_t(3, 1) == 1);
    assert(countIntegerPartitions!size_t(3, 2) == 1);
    assert(countIntegerPartitions!size_t(3, 3) == 1);
    assert(countIntegerPartitions!size_t(4, 1) == 1);
    assert(countIntegerPartitions!size_t(4, 2) == 2);
    assert(countIntegerPartitions!size_t(4, 3) == 1);
    assert(countIntegerPartitions!size_t(4, 4) == 1);
    assert(countIntegerPartitions!size_t(5, 1) == 1);
    assert(countIntegerPartitions!size_t(5, 2) == 2);
    assert(countIntegerPartitions!size_t(5, 3) == 2);
    assert(countIntegerPartitions!size_t(5, 4) == 1);
    assert(countIntegerPartitions!size_t(5, 5) == 1);
    assert(countIntegerPartitions!size_t(11, 4) == 11);
}

struct IntegerKPartitions(Int)
    if (isIntegral!Int)
{
    this(Int n, Int k)
    {
        assert(k <= n || n == 0);
        assert(k != 0);
        _buffer.length = k;
        _buffer[0] = cast(Int)(n + 1 - k);
        _buffer[1..$] = 1;
        _length = countIntegerPartitions!size_t(n, k);
        _n = n;
    }

    @property const(Int)[] front() const { return _buffer; }
    @property bool empty() const { return _empty; }
    @property size_t length() const { return _length; }

    void popFront()
    {
        --_length;
        size_t i = 1;
        for (; i < _buffer.length; ++i)
            if (_buffer[i] < _buffer[0] - 1)
                break;
        _empty = (i == _buffer.length);
        if (_empty)
            return;
        _buffer[1..i+1] = _buffer[i] + 1;
        _buffer[0] = _n - reduce!((a, b) => a + b)(_buffer[1..$]);
    }

    @property IntegerKPartitions save() const
    {
        IntegerKPartitions copy;
        copy._buffer = _buffer.dup;
        copy._empty = _empty;
        copy._n = _n;
        return copy;
    }


    // TODO: bidir
    // TODO: random
    // TODO: orders

    private Int[] _buffer;
    private size_t _length;
    private bool _empty = false;
    private Int _n;
}

/**
   Returns a range of integer partitions of $(D n) into exactly $(D k) parts,
   where ($(D k <= n)).
 */
auto integerPartitions(Int)(Int n, Int k)
{
    return IntegerKPartitions!Int(n, k);
}

unittest
{
    int[][] p0_1 = [[0]];
    int[][] p1_1 = [[1]];
    int[][] p2_1 = [[2]];
    int[][] p2_2 = [[1, 1]];
    int[][] p3_1 = [[3]];
    int[][] p3_2 = [[2, 1]];
    int[][] p3_3 = [[1, 1, 1]];
    int[][] p4_1 = [[4]];
    int[][] p4_2 = [[3, 1], [2, 2]];
    int[][] p4_3 = [[2, 1, 1]];
    int[][] p4_4 = [[1, 1, 1, 1]];
    int[][] p5_1 = [[5]];
    int[][] p5_2 = [[4, 1], [3, 2]];
    int[][] p5_3 = [[3, 1, 1], [2, 2, 1]];
    int[][] p5_4 = [[2, 1, 1, 1]];
    int[][] p5_5 = [[1, 1, 1, 1, 1]];
    int[][] p11_4 = [[8, 1, 1, 1], [7, 2, 1, 1], [6, 3, 1, 1], [5, 4, 1, 1],
                     [6, 2, 2, 1], [5, 3, 2, 1], [4, 4, 2, 1], [4, 3, 3, 1],
                     [5, 2, 2, 2], [4, 3, 2, 2], [3, 3, 3, 2]];
    assert(equal!equal(integerPartitions(0, 1), p0_1));
    assert(equal!equal(integerPartitions(1, 1), p1_1));
    assert(equal!equal(integerPartitions(2, 1), p2_1));
    assert(equal!equal(integerPartitions(2, 2), p2_2));
    assert(equal!equal(integerPartitions(3, 1), p3_1));
    assert(equal!equal(integerPartitions(3, 2), p3_2));
    assert(equal!equal(integerPartitions(3, 3), p3_3));
    assert(equal!equal(integerPartitions(4, 1), p4_1));
    assert(equal!equal(integerPartitions(4, 2), p4_2));
    assert(equal!equal(integerPartitions(4, 3), p4_3));
    assert(equal!equal(integerPartitions(4, 4), p4_4));
    assert(equal!equal(integerPartitions(5, 1), p5_1));
    assert(equal!equal(integerPartitions(5, 2), p5_2));
    assert(equal!equal(integerPartitions(5, 3), p5_3));
    assert(equal!equal(integerPartitions(5, 4), p5_4));
    assert(equal!equal(integerPartitions(5, 5), p5_5));
    assert(equal!equal(integerPartitions(11, 4), p11_4));
}

private struct SetPartitions(Range)
    if(isRandomAccessRange!Range &&
       !isInfinite!Range)
{
    this(Range r)
    {
        size_t n = walkLength(r.save);
        _ks.length = n;
        _ms.length = n;
        _masks.length = n;
        _ks[] = 0;
        _ms[] = 0;
        _masks[1..$] = 0;
        _masks[0] = (1 << n) - 1;
        _set = r;
        _length = bellNumber!size_t(cast(uint)n);
    }

    @property auto front()
    {
        return _masks.filter!(a => a != 0)().map!(a => subset(_set, a))();
    }

    @property bool empty() const { return _empty; }
    @property size_t length() const { return _length; }

    void popFront()
    {
        --_length;
        size_t n = _ks.length;
        for (size_t i = n - 1; i > 0; --i)
        {
            if (_ks[i] <= _ms[i - 1])
            {
                auto imask = 1 << i;
                _masks[_ks[i]] ^= imask;
                ++_ks[i];
                _masks[_ks[i]] ^= imask;

                _ms[i] = max(_ms[i], _ks[i]);
                for (size_t j = i + 1; j < n; ++j)
                {
                    auto jmask = 1 << j;
                    _masks[_ks[j]] ^= jmask;
                    _ks[j] = _ks[0];
                    _masks[_ks[j]] ^= jmask; // TODO: is this right!?
                    _ms[j] = _ms[i];
                }
                return;
            }
        }
        _empty = true;
        assert(_length == 0, "Reached end, but length is not zero.");
    }

    // TODO back, popBack
    // TODO opIndex

    private Range _set;
    private size_t[] _ks;
    private size_t[] _ms;
    private size_t[] _masks;
    private size_t _length;
    private bool _empty = false;
}

/**
   Returns a range of the set partitions of the range $(D r). Each set
   partition is a range of ranges, so the result of this function is a range
   of range of ranges.
*/
auto setPartitions(Range)(Range r)
{
    return SetPartitions!Range(r);
}

unittest
{
    dstring[][] abcd_parts = [["abcd"],
                              ["abc", "d"],
                              ["abd", "c"],
                              ["ab", "cd"],
                              ["ab", "c", "d"],
                              ["acd", "b"],
                              ["ac", "bd"],
                              ["ac", "b", "d"],
                              ["ad", "bc"],
                              ["a", "bcd"],
                              ["a", "bc", "d"],
                              ["ad", "b", "c"],
                              ["a", "bd", "c"],
                              ["a", "b", "cd"],
                              ["a", "b", "c", "d"]];
    alias binaryFun!((a, b) => equal!equal(a, b)) equalRoR;
    assert(equal!equalRoR(setPartitions("abcd"d), abcd_parts));
}

private struct Word(Range)
    if (isRandomAccessRange!Range &&
        hasLength!Range &&
        !isInfinite!Range)
{
    this(Range alphabet, size_t index, size_t divisor, size_t length)
    {
        _alphabet = alphabet;
        _index = index;
        _divisor = divisor;
        _length = length;
    }

    @property auto front() { return _alphabet[(_index / _divisor) % _alphabet.length]; }
    @property bool empty() const { return _length == 0; }
    @property size_t length() const { return _length; }

    void popFront()
    {
        _divisor /= _alphabet.length;
        --_length;
    }

    // TODO: opIndex, back, save

    private Range _alphabet;
    private size_t _index;
    private size_t _divisor;
    private size_t _length;
}

// TODO: this has a bug with integer overflow
private struct Words(Range)
    if (isRandomAccessRange!Range &&
        hasLength!Range &&
        !isInfinite!Range)
{
    this(Range alphabet)
    {
        assert(!alphabet.empty, "words over an empty alphabet not supported.");
        _alphabet = alphabet.save;
    }

    @property auto front() { return Word!Range(_alphabet.save, _index - _offset, _divisor, _wordLength); }
    @property enum bool empty = false;

    void popFront()
    {
        ++_index;
        if (_index - _offset >= _divisor * _alphabet.length)
        {
            _offset = _index;
            _divisor *= _alphabet.length;
            ++_wordLength;
        }
    }

    auto opIndex(size_t index)
    {
        // TODO: better algo for finding divisor
        index += _index;
        size_t offset = _offset;
        size_t divisor = _divisor;
        size_t numLetters = _alphabet.length;
        size_t wordLength = _wordLength;
        if (numLetters != 1)
        {
            size_t nextDivisor = divisor * numLetters;
            while (index - offset >= nextDivisor)
            {
                offset += nextDivisor;
                divisor = nextDivisor;
                nextDivisor *= numLetters;
                ++wordLength;
            }
        }
        else
        {
            wordLength += index;
        }
        return Word!Range(_alphabet.save, index - offset, divisor, wordLength);
    }

    private Range _alphabet;
    private size_t _index = 0;
    private size_t _divisor = 1;
    private size_t _offset = 0;
    private size_t _wordLength = 1;
}

private struct KWords(Range)
    if (isRandomAccessRange!Range &&
        hasLength!Range &&
        !isInfinite!Range)
{
    this(Range alphabet, size_t k)
    {
        assert(!alphabet.empty, "k-words over an empty alphabet not supported.");
        assert(k > 0, "0-length words not supported");

        size_t numLetters = alphabet.length;
        _alphabet = alphabet.save;
        _divisor = numLetters ^^ (k - 1);
        _length = _divisor * numLetters; // TODO: handle overflow
        _wordLength = k;
    }

    @property auto front() { return Word!Range(_alphabet.save, _index, _divisor, _wordLength); }
    @property bool empty() const { return _length == 0; }
    @property size_t length() const { return _length; }

    void popFront()
    {
        ++_index;
        --_length;
    }

    auto opIndex(size_t index)
    {
        return Word!Range(_alphabet.save, _index + index, _divisor, _wordLength);
    }

    // TODO: opSlice, save, back, popBack

    private Range _alphabet;
    private size_t _index = 0;
    private size_t _divisor = 1;
    private size_t _length = 0;
    private size_t _wordLength;
}

/**
   Returns an infinite range of all the words over $(D alphabet).
 */
auto words(Range)(Range alphabet)
{
    return Words!Range(alphabet);
}

unittest
{
    auto abc = words("abc"d);
    auto abc18 = ["a", "b", "c",
                  "aa", "ab", "ac",
                  "ba", "bb", "bc",
                  "ca", "cb", "cc",
                  "aaa", "aab", "aac",
                  "aba", "abb", "abc"];
    assert(equal!equal(abc.take(18), abc18));
    assert(equal(abc[0], "a"d));
    assert(equal(abc[1], "b"d));
    assert(equal(abc[3], "aa"d));
    assert(equal(abc[3+9], "aaa"d));
    assert(equal(abc[3+9+27], "aaaa"d));
    assert(equal(abc[3+9+27+81], "aaaaa"d));
    assert(equal(abc[3+9+27+81-1], "cccc"d));

    auto base2 = words([0, 1]);
    assert(equal(base2[0], [0]));
    assert(equal(base2[1], [1]));
    assert(equal(base2[2], [0, 0]));
    assert(equal(base2[3], [0, 1]));
    assert(equal(base2[4], [1, 0]));
    assert(equal(base2[5], [1, 1]));
    assert(equal(base2[2+4+8+16+32], [0, 0, 0, 0, 0, 0]));

    auto base1 = words("."d);
    assert(equal!equal(base1.take(5), [".", "..", "...", "....", "....."]));
    assert(equal(base1[0], "."));
    assert(equal(base1[1], ".."));
    assert(equal(base1[2], "..."));
    assert(equal(base1[3], "...."));
}

/**
   Returns a range of all $(D k)-length words over $(D alphabet).
 */
KWords!Range words(Range)(Range alphabet, size_t k)
{
    return KWords!Range(alphabet, k);
}

unittest
{
    auto abc_1 = ["a", "b", "c"];
    auto abc_2 = ["aa", "ab", "ac", "ba", "bb", "bc", "ca", "cb", "cc"];
    auto abc_3 = ["aaa", "aab", "aac", "aba", "abb", "abc", "aca", "acb", "acc",
                  "baa", "bab", "bac", "bba", "bbb", "bbc", "bca", "bcb", "bcc",
                  "caa", "cab", "cac", "cba", "cbb", "cbc", "cca", "ccb", "ccc"];

    assert(equal!equal(words("abc"d, 1), abc_1));
    assert(equal!equal(words("abc"d, 2), abc_2));
    assert(equal!equal(words("abc"d, 3), abc_3));
    assert(equal!equal(words(map!(a=>a)("abc"d), 3), abc_3));

    auto base2 = words("01"d, 32);
    assert(equal(base2[0x0000_0000U], "00000000000000000000000000000000"));
    assert(equal(base2[0x5555_5555U], "01010101010101010101010101010101"));
    assert(equal(base2[0xAAAA_AAAAU], "10101010101010101010101010101010"));
    assert(equal(base2[0xFFFF_FFFFU], "11111111111111111111111111111111"));

    auto base10 = words("0123456789"d, 10);
    assert(equal(base10[            0U], "0000000000"));
    assert(equal(base10[  999_999_999U], "0999999999"));
    assert(equal(base10[1_000_000_000U], "1000000000"));
    assert(equal(base10[1_234_567_890U], "1234567890"));
    assert(equal(base10[4_000_000_000U], "4000000000"));
    assert(equal(base10[uint.max], "4294967295")); // TODO: is this handled correctly?

    assert(equal!equal(words("."d, 1), ["."]));
    assert(equal!equal(words("."d, 10), [".........."]));
}

// TODO: belongs in std.order?
/**
Colex order
*/
bool colexOrder(alias pred = ((a, b) => a < b), Range)(Range lhs, Range rhs)
    if (isForwardRange!Range &&
        is(typeof(binaryFun!pred(lhs.front, rhs.front)) : bool))
{
    alias binaryFun!pred comp;
    alias binaryFun!((a, b) => !(comp(a, b) || comp(b, a))) equivalent;

    // TODO: probably faster to start at back for bidir ranges
    for (; !lhs.empty; lhs.popFront(), rhs.popFront())
    {
        assert(!rhs.empty, "Cannot determine colex order of ranges of different length.");
        if (comp(lhs.front, rhs.front))
        {
            Range lhsRem = lhs.save;
            Range rhsRem = rhs.save;
            lhsRem.popFront();
            rhsRem.popFront();
            if (equal!equivalent(lhsRem, rhsRem))
                return true;
        }
    }
    assert(rhs.empty, "Cannot determine colex order of ranges of different length.");
    return false;
}

unittest
{
    //123 < 124 < 134 < 234 < 125 < 135 < 235 < 145 < 245 < 345 <
    //126 < 136 < 236 < 146 < 246 < 346 < 156 < 256 < 356 < 456

    alias unaryFun!(a => isSorted!colexOrder(a)) isInColexOrder;

    int[][] colex4 = [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]];
    assert(isSorted!colexOrder(colex4));
    assert(permutations(colex4).map!(isInColexOrder)().count(true) == 1);
}

// http://theory.cs.uvic.ca/dis/programs.html (lots of algos)
// http://www.keithschwarz.com/binary-subsets/ (lex subset)
// http://www.math.dartmouth.edu/archive/m68f07/public_html/lectec.pdf
// http://www.informatik.uni-ulm.de/ni/Lehre/WS04/DMM/Software/partitions.pdf
// http://www.johndcook.com/computing_binomial_coefficients.html

// WANT:
// - ordering (lex, colex, grey, fastest)
// - set k-partitions
// - combinations
// - stirling numbers 1st + 2nd
// - fibonacci
// - integer compositions
// - algorithms for large/real values
// - k-permutations
// - integer compositions (ordered partitions)
// - set compositions
// - derangements
// - even permutations
// - range assertions
// - forward permutations
// - big O
// - const correctness
// - check over template constraints / add static asserts
// - FKM algorithm

// - random sampling
// - ranking/unranking
// - non-crossing partitions