module stdex.vectormap;

import core.exception : RangeError;
import std.array;
import std.typecons : Tuple;
import std.algorithm;

struct VectorMap(Key, Value)
{
public:
    alias Tuple!(Key, "key", Value, "value") Entry;

    this(Value[Key] init)
    {
        m_entries.length = init.length;
        size_t i = 0;
        foreach (key, value; init)
            m_entries[i++] = Entry(key, value);
        sort(m_entries);
    }

    this(Entry[] entries)
    {
        m_entries = entries.array;
        assert(entries.isSorted, "VectorMap provided with unsorted entries.");
    }

    Value opIndex(Key key)
    {
        Entry[] entries = m_entries;
        while (!entries.empty)
        {
            size_t n = entries.length;
            if (n == 1)
            {
                if (entries[0].key == key)
                    return entries[0].value;
                throw new RangeError();
            }
            size_t i = entries.length / 2;
            if (entries[i].key < key)
                entries = entries[i+1..$];
            else if (key < entries[i].key)
                entries = entries[0..i];
            else
                return entries[i].value;
        }
        throw new RangeError();
    }

private:
    Entry[] m_entries;
}

unittest
{
    import std.stdio;
    writeln("testing stdex.vectormap");

    int[string] daysInit = [
        "Sun" : 0,
        "Mon" : 1,
        "Tues" : 2,
        "Weds" : 3,
        "Thurs" : 4,
        "Fri" : 5,
        "Sat" : 6
    ];
    auto days = VectorMap!(string, int)(daysInit);
    foreach (key, value; daysInit)
        assert(days[key] == value, key);
}