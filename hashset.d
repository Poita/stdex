module stdex.hashset;

struct Node(T)
{
public:
    this(T value, Node!T* next, size_t hash)
    {
        m_value = value;
        m_next = next;
        m_hash = hash;
    }

    T m_value;
    Node* m_next;
    size_t m_hash;
}

struct HashSet(T,
               alias hash = defaultHash!T,
               alias pred = defaultPred!T)
{
public:
    this(size_t numBuckets)
    {
        m_buckets.length = numBuckets;
        m_buckets[] = null;
        m_pool = new Node!T[poolSize];
    }

    const(Node!T)* opIndex(T key)
    {
        assert(m_buckets.length != 0, "No buckets");
        size_t k = hash(key) & (m_buckets.length-1);
        Node!T* node = m_buckets[k];
        for ( ; node && (node.m_hash & (m_buckets.length-1)) == k; node = node.m_next)
            if (pred(node.m_value, key))
                return node;
        return null;
    }

    void insert(T key)
    {
        assert(m_buckets.length != 0, "No buckets");
        size_t maxK = m_buckets.length;
        size_t keyHash = hash(key);
        size_t k = hash(key) & (maxK-1);
        Node!T* node = &m_pool[0];
        m_pool = m_pool[1..$];
        if (m_pool.length == 0)
            m_pool = new Node!T[poolSize];
        *node = Node!T(key, m_buckets[k], keyHash);
        m_buckets[k] = node;
    }

private:
    enum poolSize = 1<<10;
    Node!T*[] m_buckets;
    Node!T[] m_pool;
}

size_t defaultHash(T)(T value)
{
    return typeid(T).getHash(cast(void*)&value);
}

bool defaultPred(T)(T x, T y)
{
    return x == y;
}

unittest
{
    import std.stdio;
    writeln("testing stdex.hashset");

    HashSet!int h = HashSet!int(2);
    assert(!h[0]);
    assert(!h[1]);
    assert(!h[2]);
    assert(!h[3]);
    h.insert(0);
    h.insert(1);
    h.insert(2);
    h.insert(3);
    assert(h[0]);
    assert(h[1]);
    assert(h[2]);
    assert(h[3]);
    assert(!h[4]);
    assert(!h[5]);
    assert(!h[6]);
    assert(!h[7]);
}