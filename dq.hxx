
#ifndef INCLUDE__DQ_HXX

#define INCLUDE__DQ_HXX

#include <cassert>

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

	dq_node(const dq_node<T> &orig):
			d_(orig.d_), prev_(orig.prev_), next_(orig.next_)
	{
#ifndef NDEBUG
		assert(false);
#endif
	}

	~dq_node()
	{
	}

	dq_node<T> &operator=(const dq_node<T> &rhs)
	{
#ifndef NDEBUG
		assert(false);
#endif

		if (this != (&rhs))
		{
			d_ = rhs.d_;
			prev_ = rhs.prev_;
			next_ = rhs.next_;
		}

		return *this;
	}
};

// Don't touch mah stuff.
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

	~dq()
	{
		dq_node<T> *q;

		for (dq_node<T> *iter = first_; iter != NULL; iter = q)
		{
			q = iter->next_;

			delete iter;
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
		
		++sz_;
	}

	T back()
	{
		return last_->d_;
	}

	void pop_back()
	{
#ifndef NDEBUG
		assert(last_ != NULL);
#endif

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

		--sz_;
	}

	// Make sure you're NOT using a pointer into a wrong dq.
	// 'Cause I ain't gonna catch that.
	void insert(dq_node<T> *n, T d)
	{
		dq_node<T> *p;

		if (n != NULL)
		{
			p = n->prev_;
		}
		else
		{
			p = last_;
		}

		dq_node<T> *c = new dq_node<T>(d, p, n);

		if (p != NULL)
		{
			p->next_ = c;
		}
		else
		{
			first_ = c;
		}

		if (n != NULL)
		{
			n->prev_ = c;
		}
		else
		{
			last_ = c;
		}

		++sz_;
	}
};

#endif

