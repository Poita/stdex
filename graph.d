module stdex.graph;

import std.algorithm;
import std.functional;
import std.range;
import std.traits;
import stdex.queue;
import stdex.stack;

template VertexType(Graph)
{
	alias VertexType = ParameterTypeTuple!(Graph.adjacent)[0];
}

// Default implementation of edgeCost for unweighted graphs.
size_t edgeCost(Graph, Vertex)(Graph, Vertex, Vertex)
{
	return 1;
}

auto breadthFirstSearch(Graph, Vertex)(Graph graph, Vertex start)
{
	return forwardSearch!(Queue!Vertex)(graph, start);
}

auto depthFirstSearch(Graph, Vertex)(Graph graph, Vertex start)
{
	return forwardSearch!(Stack!Vertex)(graph, start);
}

auto dijkstraSearch(Graph, Vertex)(Graph graph, Vertex start)
{
	alias typeof(graph.edgeCost(Vertex.init, Vertex.init)) Cost;
	return aStarSearch!(v => cast(Cost)0)(graph, start);
}

auto aStarSearch(alias heuristic, Graph, Vertex)(Graph graph, Vertex start)
{
	return PrioritySearch!(Graph, heuristic)(graph, start);
}

auto forwardSearch(Queue, Graph, Vertex)(Graph graph, Vertex start)
{
	return ForwardSearch!(Queue, Graph)(graph, start);
}

struct ForwardSearch(Queue, Graph)
{
public:
	alias VertexType!Graph Vertex;

	this(Graph graph, Vertex start)
	{
		m_graph = graph;
		m_queue.push(start);
		m_visited[start] = true;
	}

	@property Vertex front()
	{
		return m_queue.front;
	}

	@property bool empty() const
	{
		return m_queue.empty;
	}

	void popFront()
	{
		assert(!m_queue.empty);
		Vertex u = m_queue.front;
		m_queue.pop();
		foreach (v; m_graph.adjacent(u))
		{
			if (v !in m_visited)
			{
				m_visited[v] = true;
				m_queue.push(v);
			}
		}
	}

private:
	Graph m_graph;
	bool[Vertex] m_visited;
	Queue m_queue;
}

struct PrioritySearch(Graph, alias heuristicFunction = (v) => 0)
{
	import std.container : RedBlackTree;
	import std.typecons : Tuple;

	alias unaryFun!heuristicFunction heuristic;

public:
	alias VertexType!Graph Vertex;
	alias typeof(Graph.init.edgeCost(Vertex.init, Vertex.init)) Cost;
	alias Tuple!(Cost, "cost", Vertex, "vertex") Node;

	this(Graph graph, Vertex start)
	{
		m_graph = graph;
		m_queue = new typeof(m_queue)();
		m_queue.insert(Node(heuristic(start), start));
		m_cost[start] = 0;
	}

	@property Vertex front()
	{
		return m_queue.front.vertex;
	}

	@property bool empty()
	{
		return m_queue.empty;
	}

	void popFront()
	{
		import std.stdio;
		assert(!m_queue.empty);
		Node nu = m_queue.front;
		Vertex u = nu.vertex;
		m_queue.removeFront();
		foreach (v; m_graph.adjacent(u))
		{
			auto edgeCost = m_graph.edgeCost(u, v);
			auto newCost = m_cost[u] + edgeCost;
			auto costv = v in m_cost;
			if (!costv || newCost < *costv)
			{
				auto hv = heuristic(v);
				if (costv)
					m_queue.removeKey(Node(*costv - hv, v));
				m_cost[v] = newCost;
				m_parent[v] = u;
				m_queue.insert(Node(newCost + hv, v));
			}
		}
	}

	@property auto path()
	{
		Vertex v = front;
		Vertex[] p = [v];
		while (v in m_parent)
		{
			v = m_parent[v];
			p ~= v;
		}
		reverse(p);
		return p;
	}

private:
	Graph m_graph;
	RedBlackTree!(Node) m_queue;
	Cost[Vertex] m_cost;
	Vertex[Vertex] m_parent;
}

unittest
{
	import std.stdio;
	writeln("testing stdex.graph");

	int[][] maze = [
		[1,1,1,1,1,1,1,1],
		[1,0,0,1,0,0,0,1],
		[1,0,0,0,0,1,0,1],
		[1,1,0,0,0,1,0,1],
		[1,1,1,1,1,0,1,1]
	];

	struct Graph
	{
	public:
		int[2][] adjacent(int[2] u)
		{
			immutable int[2][] deltas = [[-1,0], [1,0], [0,-1], [0,1]];
			int[2][] a;
			foreach (d; deltas)
			{
				int[2] v = [u[0]+d[0], u[1]+d[1]];
				if (v[0] < 0 || v[1] < 0) continue;
				if (v[0] >= m.length || v[1] >= m[0].length) continue;
				if (m[v[0]][v[1]]) continue;
				a ~= v;
			}
			return a;
		}

		int edgeCost(int[2] u, int[2] v)
		{
			return 1;
		}

		int[][] m;
	}

	Graph g = Graph(maze);
	int[2] u = [1, 1];
	import std.typetuple;
	foreach (search; TypeTuple!(breadthFirstSearch, depthFirstSearch, dijkstraSearch))
	{
		assert(search!(v => v == [3, 6])(g, u));
		assert(!search!(v => v == [4, 6])(g, u));
		assert(search!(v => v == [1, 6])(g, u));
		assert(search!(v => v == [1, 1])(g, u));
		assert(search!(v => v == [2, 3])(g, u));
		assert(!search!(v => v == [4, 5])(g, u));
	}
}

auto implicitGraph(Vertex, alias _adjacent, alias _edgeCost = (u, v) => cast(size_t)1)()
{
	struct Result
	{
		auto adjacent(Vertex v)
		{
			return _adjacent(v);
		}

		auto edgeCost(Vertex u, Vertex v)
		{
			return _edgeCost(u, v);
		}
	}
	return Result();
}

unittest
{
	import std.stdio, std.math, std.functional;
	writeln("testing stdex.graph");

	auto graph = implicitGraph!(int, u => [u + 1, u + 2]);
	assert(aStarSearch!(v => abs(v - 10) / 2)(graph, 1).find(10).path == [1, 3, 5, 7, 9, 10]);
}