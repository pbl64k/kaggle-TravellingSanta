
#undef NDEBUG

#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <deque>
#include <functional>
#include <list>
#include <numeric>
#include <random>
#include <unordered_set>
#include <utility>
#include <vector>

#include "kdt.hxx"
#include "dq.hxx"
#include "qxt.hxx"

using namespace std;

#undef DIAG_IMP_L
#define DIAG_PATHG_V
#define DIAG_PATHG_SUM
#define DIAG_OPTS_O
#define DIAG_OPTS_OMOD 1000
#undef DIAG_OPTS_GAIN
#undef DIAG_OPTS_WHOOPS
#define DIAG_OPTS_SUM
#define DIAG_OPTX_ITER
#define DIAG_OPTX_ITSC
#define DIAG_OPTX_INIT
#define DIAG_OPTX_SUM

#define PROB_SZ 150000
#define DIM 2
#define PATH_NUM 2
#define OPT_ITERS 300000
#define OPT_NN_K 64
#define GAIN_THRESH 10000
#define XOPT_ITERS 1
#define XOPT_TGT 8000000

typedef size_t vertex_id;
typedef unsigned long long edge;
typedef unsigned long long quad;

edge make_edge(vertex_id a, vertex_id b)
{
	return (static_cast<edge>(a) << 18) | static_cast<edge>(b);
}

quad make_quad(vertex_id a, vertex_id b, vertex_id c, vertex_id n)
{
	return static_cast<quad>(a) | (static_cast<quad>(b) << 18) | (static_cast<quad>(c) << 36) | (static_cast<quad>(n) << 54);
}

vector<point> cities(PROB_SZ);

class path
{
	public:
	double dist_;
	vertex_id cur_;
	dq<vertex_id> p_;
	vector<dq_node<vertex_id> *> vs_;
	vector<list<vertex_id> > es_;
	qxt<pair<vertex_id, vertex_id> > qxes_;

	path(): dist_(0.0), cur_(0), p_(), vs_(PROB_SZ), es_(PROB_SZ),
			qxes_()
	{
	}

	// qxt blows up on copy and assignment. just as planned.
	path(const path &orig): dist_(orig.dist_), cur_(orig.cur_),
			p_(orig.p_), vs_(orig.vs_), es_(orig.es_), qxes_()
	{
	}

	~path()
	{
	}

	path &operator=(const path &rhs)
	{
		if (this != (&rhs))
		{
			dist_ = rhs.dist_;
			cur_ = rhs.cur_;
			p_ = rhs.p_;
			vs_ = rhs.vs_;
			es_ = rhs.es_;
		}

		return *this;
	}

	void init(vertex_id city)
	{
		cur_ = city;
		p_.push_back(cur_);
		vs_[cur_] = p_.first_;
	}

	void add(const pair<coord, vertex_id> &pt)
	{
		es_[cur_].push_back(pt.second);
		es_[pt.second].push_back(cur_);
		qxes_.add_edge(make_pair(cities[cur_], cities[pt.second]), make_pair(cur_, pt.second));
		dist_ += pt.first;
		cur_ = pt.second;
		p_.push_back(cur_);
		vs_[cur_] = p_.last_;
	}
};

kdt<DIM, vertex_id> s_kdt;

path *pth_ptr;
vector<bool> *v_cities_ptr;
unordered_set<edge> *blacklist_ptr;

coord calc_dist_v(vertex_id a, vertex_id b)
{
	return calc_dist(cities[a], cities[b]);
}

bool ff(const vertex_id &x)
{
	return (! (*v_cities_ptr)[x]) && (blacklist_ptr->find(make_edge(pth_ptr->cur_, x)) == blacklist_ptr->end());
}

void gen_path(path *pthp, size_t vert_num, const kdt<DIM, vertex_id> &s_kdt, unordered_set<edge> &blacklist, size_t cand_num)
{
	vector<bool> v_cities(vert_num);

	random_device rng;
	uniform_int_distribution<vertex_id> rng_dst(0, vert_num - 1);

	pthp->init(rng_dst(rng));

	v_cities[pthp->cur_] = true;

	pth_ptr = pthp;
	v_cities_ptr = &v_cities;
	blacklist_ptr = &blacklist;

	for (size_t i = 1; i != vert_num; ++i)
	{
		deque<pair<coord, vertex_id> > cnds = s_kdt.find_knn(cities[pthp->cur_], cand_num, ff);

		vertex_id last_cur = pthp->cur_;

		size_t j = 0;
		coord w = 0;
		vector<coord> ws(cand_num);

		deque<pair<coord, vertex_id> >::const_iterator iter;

		for (iter = cnds.begin(); iter != cnds.end(); ++iter)
		{
			w += 1.0 / ((iter->first) * (iter->first));
			ws[j] = w;
			++j;
		}

		uniform_real_distribution<coord> rng_dst(0.0, 1.0);

		coord p = rng_dst(rng) * w;

		j = 0;

		for (iter = cnds.begin(); ws[j] < p; ++iter)
		{
			++j;
		}

		if (iter == cnds.end())
		{
			--iter;
		}

		pthp->add(*iter);

		v_cities[pthp->cur_] = true;
		blacklist.insert(make_edge(last_cur, pthp->cur_));
		blacklist.insert(make_edge(pthp->cur_, last_cur));

#ifdef DIAG_PATHG_V
		cout << (i + 1) << " dist: " << pthp->dist_ << endl;
#endif
	}
}

unordered_set<quad> q_exc;

void opt_path(const kdt<DIM, vertex_id> &s_kdt, path &pth, unordered_set<edge> &blacklist, size_t k, vertex_id seed)
{
	pth_ptr = &pth;

	(*v_cities_ptr)[pth.cur_] = true;
	deque<pair<coord, vertex_id> > cnds = s_kdt.find_knn(cities[pth.cur_], k, ff);
	(*v_cities_ptr)[pth.cur_] = false;

	coord bestgain = -100500;
	vertex_id bestnn1 = 0;
	vertex_id bestnn2 = 0;

	for (deque<pair<coord, vertex_id> >::const_iterator citer = cnds.begin(); citer != cnds.end(); ++citer)
	{
		vertex_id nn1 = citer->second;

		for (list<vertex_id>::const_iterator v_citer = pth.es_[nn1].begin(); v_citer != pth.es_[nn1].end(); ++v_citer)
		{
			vertex_id nn2 = (*v_citer);

			if ((blacklist.find(make_edge(pth.cur_, nn2)) != blacklist.end()) || (q_exc.find(make_quad(pth.cur_, nn1, nn2, seed)) != q_exc.end()))
			{
				continue;
			}

			coord gain = calc_dist_v(pth.cur_, pth.p_.last_->prev_->d_) + calc_dist_v(nn1, nn2)
					- calc_dist_v(pth.cur_, nn1) - calc_dist_v(pth.cur_, nn2);

			if (gain > bestgain)
			{
				bestgain = gain;
				bestnn1 = nn1;
				bestnn2 = nn2;
			}
		}
	}

#ifdef DIAG_OPTS_GAIN
	cout << "BG:" << bestgain << " ";
#endif

	if (bestgain < -20000)
	{
#ifdef DIAG_OPTS_WHOOPS
		cout << "Whoops!" << endl;
#endif
		return;
	}

	q_exc.insert(make_quad(pth.cur_, bestnn1, bestnn2, seed));

	pth.dist_ -= bestgain;

	coord cur = pth.p_.back();
	pth.p_.pop_back();
	coord new_cur = pth.p_.back();
	pth.cur_ = new_cur;

	blacklist.erase(blacklist.find(make_edge(cur, new_cur)));
	blacklist.erase(blacklist.find(make_edge(new_cur, cur)));
	blacklist.erase(blacklist.find(make_edge(bestnn1, bestnn2)));
	blacklist.erase(blacklist.find(make_edge(bestnn2, bestnn1)));
	blacklist.insert(make_edge(bestnn1, cur));
	blacklist.insert(make_edge(cur, bestnn1));
	blacklist.insert(make_edge(bestnn2, cur));
	blacklist.insert(make_edge(cur, bestnn2));

	pth.es_[cur].remove(new_cur);
	pth.es_[new_cur].remove(cur);
	pth.es_[bestnn1].remove(bestnn2);
	pth.es_[bestnn2].remove(bestnn1);
	pth.es_[bestnn1].push_back(cur);
	pth.es_[cur].push_back(bestnn1);
	pth.es_[bestnn2].push_back(cur);
	pth.es_[cur].push_back(bestnn2);

	dq_node<vertex_id> *ins = pth.vs_[bestnn1]->next_;

	if ((ins == NULL) || (ins->d_ != bestnn2))
	{
		ins = pth.vs_[bestnn2]->next_;
	}

	pth.p_.insert(ins, cur);

	pth.vs_[cur] = ins->prev_;

	/*
	// SANITY
	dq_node<vertex_id> *a = pth.p_.first_;
	dq_node<vertex_id> *b = pth.p_.first_->next_;

	for (int kkk = 0; kkk != PROB_SZ - 1; ++kkk)
	{
		if (blacklist.find(make_edge(a->d_, b->d_)) == blacklist.end())
		{
			cout << "th: " << kkk << " " << a->d_ << " " << b->d_ << endl;
			assert(false);
		}

		if (blacklist.find(make_edge(b->d_, a->d_)) == blacklist.end())
		{
			cout << "fr: " << kkk << " " << b->d_ << " " << a->d_ << endl;
			assert(false);
		}

		a = b;
		b = b->next_;
	}
	*/
}

void simul_opt(const kdt<DIM, vertex_id> &s_kdt, path ***paths, unordered_set<edge> &blacklist, size_t maxiters, size_t k)
{
	vector<bool> v_cities(PROB_SZ);
	path *bestpaths[PATH_NUM];
	coord bestscore = 0;

	for (int l = 0; l != PATH_NUM; ++l)
	{
		bestpaths[l] = new path();
		(*(bestpaths[l])) = (*((*paths)[l]));
	}

	v_cities_ptr = &v_cities;
	blacklist_ptr = &blacklist;

	for (size_t i = 0; i != PATH_NUM; ++i)
	{
		bestscore = max(bestscore, (*paths)[i]->dist_);
	}

	coord oldbest = bestscore;

	for (size_t i = 0; i != maxiters; ++i)
	{
		size_t n = 0;
		coord dn = (*paths)[n]->dist_;

		for (size_t nn = 1; nn != PATH_NUM; ++nn)
		{
			if ((*paths)[nn]->dist_ > dn)
			{
				n = nn;
				dn = (*paths)[n]->dist_;
			}
		}

		opt_path(s_kdt, *((*paths)[n]), blacklist, k, n);

		coord curscore = 0;

		for (size_t j = 0; j != PATH_NUM; ++j)
		{
			curscore = max(curscore, (*paths)[j]->dist_);
		}

		oldbest = bestscore;

		if (curscore <= (bestscore - GAIN_THRESH))
		{
			bestscore = curscore;

			for (int l = 0; l != PATH_NUM; ++l)
			{
				(*(bestpaths[l])) = (*((*paths)[l]));
			}
		}

#ifdef DIAG_OPTS_O
#ifdef DIAG_OPTS_OMOD
		if (((i + 1) % DIAG_OPTS_OMOD) == 0)
		{
#endif
		cout << (i + 1) << " " << (n + 1) << " best: " << bestscore << " gain: " << (oldbest - bestscore) << " cur: " << curscore
				<< " 1: " << (*paths)[0]->dist_ << " 2: " << (*paths)[1]->dist_ << endl;
#ifdef DIAG_OPTS_OMOD
		}
#endif
#endif
	}

	for (int l = 0; l != PATH_NUM; ++l)
	{
		(*((*paths)[l])) = (*(bestpaths[l]));
	}
}

void opt_x(const kdt<DIM, vertex_id> &s_kdt, path &pth, unordered_set<edge> &blacklist, vertex_id seed)
{
	random_device rng;
	uniform_int_distribution<vertex_id> rng_dst(0, PROB_SZ - 1);

	vertex_id start = rng_dst(rng);

#ifdef DIAG_OPTX_INIT
	cout << "Initial vertex: " << start << endl;
#endif

	/*
	for ()
	{
	}
	*/
}

void x_opt(const kdt<DIM, vertex_id> &s_kdt, path ***paths, unordered_set<edge> &blacklist, size_t maxiters, size_t tgt)
{
	size_t iter_num = 0;
	bool tgt_met = false;

	while ((iter_num < maxiters) && (! tgt_met))
	{
#ifdef DIAG_OPTX_ITER
		cout << "Iteration #" << (iter_num + 1) << endl;
#endif

		size_t n = 0;
		coord dn = (*paths)[n]->dist_;

		for (size_t nn = 1; nn != PATH_NUM; ++nn)
		{
			if ((*paths)[nn]->dist_ > dn)
			{
				n = nn;
				dn = (*paths)[n]->dist_;
			}
		}

		opt_x(s_kdt, *((*paths)[n]), blacklist, n);

		coord curscore = 0;

		for (size_t i = 0; i != PATH_NUM; ++i)
		{
			curscore = max(curscore, (*paths)[i]->dist_);
		}

#ifdef DIAG_OPTX_ITSC
		cout << "Score: " << curscore << endl;
#endif

		tgt_met = curscore <= tgt;

		++iter_num;
	}
}

int main()
{
	srand(time(NULL));

	for (size_t i = 0; i != PROB_SZ; ++i)
	{
		vertex_id id;
		double x, y;

		cin >> id >> x >> y;

#ifdef DIAG_IMP_L
		cout << "Importing line " << (i + 1) << " id: " << id << " (" << x << ", " << y << ")" << endl;
#endif

		point pt(DIM);

		pt[0] = x;
		pt[1] = y;

		cities[id] = pt;

		s_kdt.add(pt, id);
	}

	path **paths = new path *[PATH_NUM];

	unordered_set<edge> blacklist;

	for (int i = PATH_NUM; i >= 1; --i)
	{
		paths[i - 1] = new path();
		gen_path(paths[i - 1], PROB_SZ, s_kdt, blacklist, i);
	}

#ifdef DIAG_PATHG_SUM
	cout << "Pathgen complete: " << paths[0]->dist_ << " " << paths[1]->dist_ << endl;
#endif

/*
	simul_opt(s_kdt, &paths, blacklist, OPT_ITERS, OPT_NN_K);

#ifdef DIAG_OPTS_SUM
	cout << "Simulopt complete: " << paths[0]->dist_ << " " << paths[1]->dist_ << endl;
#endif
*/

	x_opt(s_kdt, &paths, blacklist, XOPT_ITERS, XOPT_TGT);

#ifdef DIAG_OPTX_SUM
	cout << "Xopt complete: " << paths[0]->dist_ << " " << paths[1]->dist_ << endl;
#endif

	fstream csv("./sol.csv", ios_base::out);

	csv << "path1,path2" << endl;

	dq_node<vertex_id> *p1v = paths[0]->p_.first_;
	dq_node<vertex_id> *p2v = paths[1]->p_.first_;

	for (size_t i = 0; i != PROB_SZ; ++i)
	{
		csv << (p1v->d_) << "," << (p2v->d_) << endl;

		p1v = p1v->next_;
		p2v = p2v->next_;
	}

	return 0;
}

