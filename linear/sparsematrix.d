module stdex.linear.sparsematrix;

import std.algorithm;
import std.array;
import std.range;
import std.typecons;
import stdex.range;
import stdex.util;
import stdex.vectormap;

struct SparseVector(Element)
{
public:
    this(Element[] values, size_t[] indices)
    {
        m_values = values;
        m_indices = indices;
    }

    Element opIndex(size_t index)
    {
        assert(0, "florp");
        return 0;
    }

    @property auto keyValues()
    {
        return zip(m_indices, m_values);
    }

private:
    Element[] m_values;
    size_t[] m_indices;
}

struct SparseMatrix(Element)
{
public:
    alias size_t[2] Key;

    this(size_t nrows, size_t ncolumns, Element[Key] init)
    {
        assert(nrows != 0, "Cannot have a matrix with 0 rows.");
        assert(ncolumns != 0, "Cannot have a matrix with 0 columns.");

        m_nrows = nrows;
        m_ncolumns = ncolumns;

        auto kvs = init.keyValueArray();
        m_values.length = kvs.length;
        m_columnIndex.length = kvs.length;
        m_rowOffset.length = nrows + 1;

        sort(kvs);
        auto rows = kvs.groupSlices!((a, b) => a[0][0] == b[0][0]);
        size_t i = 0;
        size_t offset = 0;
        foreach (row; rows)
        {
            assert(!row.empty, "Should not have an empty group slice.");
            size_t rowIndex = row.front[0][0];
            assert(rowIndex >= i, "Rows should be increasing.");
            m_rowOffset[i..rowIndex + 1] = offset;
            row.map!(e => e[1]).copy(m_values[offset..$]);
            row.map!(e => e[0][1]).copy(m_columnIndex[offset..$]);
            offset += row.length;
            i = rowIndex + 1;
        }
        m_rowOffset[i..$] = offset;
    }

    SparseVector!Element opIndex(size_t index)
    {
        assert(index < m_nrows, "Index out of bounds.");
        size_t i = m_rowOffset[index];
        size_t j = m_rowOffset[index + 1];
        return SparseVector!Element(m_values[i..j], m_columnIndex[i..j]);
    }

    Element opIndex(size_t row, size_t col)
    {
        assert(row < m_nrows, "Index out of bounds.");
        assert(col < m_ncolumns, "Index out of bounds.");
        size_t i = m_rowOffset[row];
        size_t j = m_rowOffset[row + 1];
        Element[] vals = m_values[i..j];
        size_t[] colIdx = m_columnIndex[i..j];
        auto idx = colIdx.countUntil(col);
        if (idx == -1)
            return cast(Element) 0;
        return vals[idx];
    }

    @property auto nonEmptyRows()
    {
        return iota(m_nrows).
               filter!(i => m_rowOffset[i] != m_rowOffset[i+1]).
               map!(i => tuple(i, this[i]));
    }

private:
    Element[] m_values;
    size_t[] m_columnIndex;
    size_t[] m_rowOffset;
    size_t m_nrows;
    size_t m_ncolumns;
}

auto sparseMatrix(Element)(size_t nrows, size_t ncolumns, Element[size_t[2]] init)
{
    return SparseMatrix!(Element)(nrows, ncolumns, init);
}

auto multiply(Element)(auto ref SparseMatrix!Element lhs, auto ref SparseMatrix!Element rhs)
{
    assert(lhs.m_ncolumns == rhs.m_nrows, "LHS columns must match RHS rows.");

    // ABij = Aik * Bkj

    Element[size_t[2]] result;
    foreach (i; 0..lhs.m_nrows)
    {
        foreach (j; 0..rhs.m_ncolumns)
        {
            Element e = 0;
            foreach (k; 0..lhs.m_ncolumns)
            {
                e += lhs[i, k] * rhs[k, j];
            }
            if (e != 0)
            {
                size_t[2] key = void;
                key[0] = i;
                key[1] = j;
                result[key] = e;
            }
        }
    }
    return result;
}

unittest
{
    import std.stdio;
    writeln("testing stdex.linear.sparsematrix.SparseMatrix");

    double[size_t[2]] aa = [
        [0, 0] : 1.0,
        [1, 1] : 1.0,
        [2, 2] : 1.0 /+,
        [3, 3] : 1.0,
        [1, 3] : 3.0,
        [2, 3] : 5.0+/
    ];

    auto a = sparseMatrix(3, 3, aa);
    auto b = sparseMatrix(3, 3, aa);
    //writeln(multiply(a, b));
}