module stdex.heap;

import std.functional;

void pushHeap(alias compare = "a < b", T)(T[] heap)
{
	pushHeap!compare(heap, heap.length - 1, 0, heap[$ - 1]);
}

void popHeap(alias compare = "a < b", T)(T[] heap)
{
	if (heap.length > 1)
		popHeap!compare(heap[0..$-1], heap[$-1]);
}

void makeHeap(alias compare = "a < b", T)(T[] heap)
{
	if (heap.length < 2)
		return;

	size_t length = heap.length;
	size_t parent = (length - 2) / 2;

	while (true)
	{
		T value = heap[parent];
		adjustHeap!compare(heap, parent, value);
		if (parent == 0)
			return;
		--parent;
	}
}

void adjustHeap(alias compare = "a < b", T)(T[] heap, size_t hole, T value)
{
	alias binaryFun!compare comp;
	const size_t top = hole;
	size_t second = hole;
	size_t length = heap.length;
	while (second < (length - 1) / 2)
	{
  		second = 2 * (second + 1);
  		if (comp(heap[second], heap[second - 1]))
    	second--;
    	heap[hole] = heap[second];
    	hole = second;
	}
  	if ((length & 1) == 0 && second == (length - 2) / 2)
	{
	  	second = 2 * (second + 1);
	  	heap[hole] = heap[second - 1];
	  	hole = second - 1;
	}
	pushHeap!compare(heap, hole, top, value);
}

private void popHeap(alias compare = "a < b", T)(T[] heap, ref T result)
{
	T value = result;
	result = heap[0];
	adjustHeap!compare(heap, 0, value);
}

private void pushHeap(alias compare = "a < b", T)(T[] heap, size_t hole, size_t top, T value)
{
	alias binaryFun!compare comp;
	size_t parent = (hole - 1) / 2;
	while (hole > top && comp(heap[parent], value))
	{
		heap[hole] = heap[parent];
		hole = parent;
		parent = (hole - 1) / 2;
	}
	heap[hole] = value;
}