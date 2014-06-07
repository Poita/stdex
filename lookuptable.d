module stdex.lookuptable;

import std.algorithm;
import std.array;
import std.typecons;

struct LookupTable(Key, Value, alias hash)
{
public:
    alias KeyValue = Tuple!(Key, "key", Value, "value");

    this(Value[Key] entries)
    {
        size_t n = entries.length;
        m_entries = new KeyValue[n];
        m_indices = new size_t[n + 1];
        
        size_t i = 0;
        foreach (key, value; entries)
            m_entries[i++] = KeyValue(key, value);
        sort!((a, b) => bucketIndex(a.key) < bucketIndex(b.key))(m_entries);

        i = 0;
        foreach (b; 0..m_indices.length)
        {
            m_indices[b] = i;
            while (i < m_entries.length && bucketIndex(m_entries[i].key) == b)
                ++i;
        }
    }

    inout(Value)* opBinaryRight(string op : "in")(Key key) inout
    {
        size_t b = bucketIndex(key);
        foreach (i; m_indices[b]..m_indices[b+1])
            if (m_entries[i].key == key)
                return &m_entries[i].value;
        return null;
    }

    ref inout(Value) opIndex(Key key) inout
    {
        auto p = key in this;
        assert(p !is null, "Key not in lookup table.");
        return *p;
    }

    @property size_t length() const
    {
        return m_entries.length;
    }

    size_t bucketIndex(Key key) const
    {
        return hash(key) % (m_indices.length - 1);
    }

private:
    KeyValue[] m_entries;
    size_t[] m_indices;
}

auto lookupTable(alias hash = (ref x) => typeid(Key).getHash(&x), Key, Value)(Value[Key] entries)
{
    return LookupTable!(Key, Value, hash)(entries);
}

unittest
{
    auto t = lookupTable([1:"one", 2:"two", 3:"three"]);
    assert(t[1] == "one");
    assert(t[2] == "two");
    assert(t[3] == "three");
    assert(0 !in t);
    assert(1 in t);
    assert(2 in t);
    assert(3 in t);
    assert(4 !in t);
    assert(t.length == 3);
}