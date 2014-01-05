module stdex.stack;

struct Stack(T)
{
public:
	@property bool empty() const
	{
		return m_top == 0;
	}

	@property size_t length() const
	{
		return m_top;
	}

	@property ref const(T) front() const
	{
		assert(m_top, "Stack is empty.");
		return m_buffer[m_top - 1];
	}

	void push(T value)
	{
		if (m_buffer.length == m_top)
			m_buffer.length = m_buffer.length == 0 ? 8 : m_buffer.length * 2;
		m_buffer[m_top++] = value;
	}

	void pop()
	{
		assert(m_top, "Stack is empty.");
		--m_top;
	}


private:
	T[] m_buffer;
	size_t m_top;
}