
#ifndef INCLUDE__DQ_HXX

#define INCLUDE__DQ_HXX

template<typename T> class dq_node
{
	public:
	T d_;
	dq_node<T> *prev_;
	dq_node<T> *next_;

	dq_node(T d, dq_node<T> *prev, dq_node<T> *next):
			d_(d), prev_(prev), next_(next)
	{
	}
};

template<typename T> class dq
{
	public:
	size_t sz_;
	dq_node<T> *first_;
	dq_node<T> *last_;

	dq(): sz_(0), first_(NULL), last_(NULL)
	{
	}

	dq(const dq<T> &orig): sz_(0), first_(NULL), last_(NULL)
	{
		for (dq_node<T> *iter = orig.first_; iter != NULL; iter = iter->next_)
		{
			push_back(iter->d_);
		}
	}

	dq<T> &operator=(const dq<T> &rhs)
	{
		if (this != (&rhs))
		{
			dq_node<T> *q;

			for (dq_node<T> *iter = first_; iter != NULL; iter = q)
			{
				q = iter->next_;

				delete iter;
			}

			sz_ = 0;
			first_ = NULL;
			last_ = NULL;

			for (dq_node<T> *iter = rhs.first_; iter != NULL; iter = iter->next_)
			{
				push_back(iter->d_);
			}
		}

		return *this;
	}

	void push_back(T d)
	{
		dq_node<T> *n = new dq_node<T>(d, last_, NULL);

		if (first_ == NULL)
		{
			first_ = last_ = n;
		}
		else
		{
			last_->next_ = n;
			last_ = n;
		}
	}

	T back()
	{
		return last_->d_;
	}

	void pop_back()
	{
		dq_node<T> *cur_last = last_;

		last_ = cur_last->prev_;

		if (last_ == NULL)
		{
			first_ = NULL;
		}
		else
		{
			last_->next_ = NULL;
		}

		delete cur_last;
	}

	void insert(dq_node<T> *n, T d)
	{
		dq_node<T> *p = n->prev_;
		dq_node<T> *c = new dq_node<T>(d, p, n);
		p->next_ = c;
		n->prev_ = c;
	}
};

#endif

