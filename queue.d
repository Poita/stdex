module stdex.queue;

import std.array;
import stdex.heap;

struct Queue(T)
{
public:
	@property bool empty() const
	{
		return m_first == m_last;
	}

	@property size_t length() const
	{
		return m_length;
	}

	@property ref const(T) front() const
	{
		assert(m_length, "Queue is empty.");
		assert(m_first < m_buffer.length, "Internal error: m_first is out of range.");
		return m_buffer[m_first];
	}

	void push(T value)
	{
		if (m_length == m_buffer.length)
		{
			// Resize
			size_t n = m_length;
			m_buffer.length = n == 0 ? 8 : n * 2;
			size_t m = n - m_first;
			if (m)
			{
				m_buffer[$-m..$] = m_buffer[m_first..n];
				m_first = m_buffer.length - m;
			}
		}
		m_buffer[m_last++] = value;
		++m_length;
		if (m_last == m_buffer.length)
			m_last = 0;
	}

	void pop()
	{
		assert(m_length, "Queue is empty.");
		++m_first;
		if (m_first == m_buffer.length)
			m_first = 0;
		--m_length;
	}

private:
	T[] m_buffer = [];
	size_t m_first = 0;
	size_t m_last = 0;
	size_t m_length = 0;
}

unittest
{
	import std.stdio;
	writeln("testing stdex.queue.Queue");

	Queue!int q;
	foreach (i; 0..1000)
	{
		q.push(2*i);
		q.push(2*i+1);
		assert(q.front == i);
		q.pop();
		assert(q.length == i+1);
	}
}

struct PriorityQueue(T, alias compare = "a < b")
{
public:
	@property bool empty() const
	{
		return m_heap.empty;
	}

	@property size_t length() const
	{
		return m_heap.length;
	}

	@property ref const(T) front() const
	{
		return m_heap[0];
	}

	void push(T value)
	{
		m_heap ~= value;
		pushHeap!compare(m_heap);
	}

	void pop()
	{
		popHeap!compare(m_heap);
		m_heap = m_heap[0..$-1];
	}

private:
	T[] m_heap;
}

unittest
{
	import std.stdio;
	writeln("testing stdex.queue.PriorityQueue");

	PriorityQueue!int q;
	int[] a = [5, 4, 6, 3, 7, 2, 8, 1, 9];
	int[] b = [5, 5, 6, 6, 7, 7, 8, 8, 9];
	foreach (i; 0..9)
	{
		q.push(a[i]);
		assert(q.front == b[i]);
	}
	foreach (i; 0..9)
	{
		assert(q.front == 9 - i);
		q.pop();
	}
}