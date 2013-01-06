
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
#define DIAG_OPTX_CAND
#define DIAG_OPTX_SUM

#define PROB_SZ 150000
#define DIM 2
#define PATH_NUM 2
#define OPT_ITERS 300000
#define OPT_NN_K 64
#define GAIN_THRESH 10000
#define XOPT_ITERS 1000000
#define XOPT_TGT 7500000

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
	vector<vertex_id> vvs_;
	vector<list<vertex_id> > es_;
	qxt<pair<vertex_id, vertex_id> > qxes_;

	path(): dist_(0.0), cur_(0), p_(), vs_(PROB_SZ), vvs_(PROB_SZ), es_(PROB_SZ),
			qxes_()
	{
	}

	// qxt blows up on copy and assignment. just as planned.
	path(const path &orig): dist_(orig.dist_), cur_(orig.cur_),
			p_(orig.p_), vs_(orig.vs_), vvs_(orig.vvs_), es_(orig.es_), qxes_()
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
			vvs_ = rhs.vvs_;
			es_ = rhs.es_;
		}

		return *this;
	}

	void init(vertex_id city)
	{
		cur_ = city;
		p_.push_back(cur_);
		vs_[cur_] = p_.first_;
		vvs_[cur_] = p_.sz_ - 1;
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
		vvs_[cur_] = p_.sz_ - 1;
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

bool opt_x(const kdt<DIM, vertex_id> &s_kdt, path &pth, unordered_set<edge> &blacklist, vertex_id seed)
{
	random_device rng;
	uniform_int_distribution<vertex_id> rng_dst(0, PROB_SZ - 1);

	vertex_id start = rng_dst(rng);

#ifdef DIAG_OPTX_INIT
	cout << "Initial vertex: " << start << endl;
#endif

	for (vertex_id a0 = start; a0 != start + PROB_SZ; ++a0)
	{
		vertex_id a = a0 % PROB_SZ;

		for (list<vertex_id>::iterator iter = pth.es_[a].begin(); iter != pth.es_[a].end(); ++iter)
		{
			set<pair<vertex_id, vertex_id> > xs = pth.qxes_.find_x(make_pair(cities[a], cities[*iter]));

			for (set<pair<vertex_id, vertex_id> >::iterator s_iter = xs.begin(); s_iter != xs.end(); ++s_iter)
			{
				if ((s_iter->first == a) || (s_iter->first == (*iter)) || (s_iter->second == a) || (s_iter->second == (*iter)))
				{
					continue;
				}

				vertex_id v1, v2, v3, v4;

				if (pth.vvs_[a] < pth.vvs_[s_iter->first])
				{
					v1 = a;
					v2 = (*iter);
					v3 = s_iter->first;
					v4 = s_iter->second;
				}
				else
				{
					v1 = s_iter->first;
					v2 = s_iter->second;
					v3 = a;
					v4 = (*iter);
				}

				if (pth.vvs_[v1] > pth.vvs_[v2])
				{
					swap(v1, v2);
				}

				if (pth.vvs_[v3] > pth.vvs_[v4])
				{
					swap(v3, v4);
				}

				if ((blacklist.find(make_edge(v1, v3)) != blacklist.end()) ||
						(blacklist.find(make_edge(v2, v4)) != blacklist.end()))
				{
					continue;
				}

#ifdef DIAG_OPTX_CAND
				cout << "Candidate found: " << a << "-" << (*iter) << " intersects with " << s_iter->first << "-" << s_iter->second << endl;
				cout << "(" << v1 << "-" << v2 << ", " << v3 << "-" << v4 << ")" << endl;
				cout << "(" << pth.vvs_[v1] << "-" << pth.vvs_[v2] << ", " << pth.vvs_[v3] << "-" << pth.vvs_[v4] << ")" << endl;
				cout << "(" << v1 << "-" << v3 << ", " << v2 << "-" << v4 << ")" << endl;
				cout << "Exchanging edges..." << endl;
#endif

#ifndef NDEBUG
				/*
				if (pth.vvs_[v1] + 1 != pth.vvs_[v2])
				{
					cout << v1 << " ";
					for (list<vertex_id>::iterator it1 = pth.es_[v1].begin(); it1 != pth.es_[v1].end(); ++it1)
					{
						cout << (*it1) << " ";
					}
					cout << endl;
					cout << v2 << " ";
					for (list<vertex_id>::iterator it1 = pth.es_[v2].begin(); it1 != pth.es_[v2].end(); ++it1)
					{
						cout << (*it1) << " ";
					}
					cout << endl;
				}

				if (pth.vvs_[v3] + 1 != pth.vvs_[v4])
				{
					cout << v3 << " ";
					for (list<vertex_id>::iterator it1 = pth.es_[v3].begin(); it1 != pth.es_[v3].end(); ++it1)
					{
						cout << (*it1) << " ";
					}
					cout << endl;
					cout << v4 << " ";
					for (list<vertex_id>::iterator it1 = pth.es_[v4].begin(); it1 != pth.es_[v4].end(); ++it1)
					{
						cout << (*it1) << " ";
					}
					cout << endl;
				}
				*/

				assert(pth.vvs_[v1] + 1 == pth.vvs_[v2]);
				assert(pth.vvs_[v3] + 1 == pth.vvs_[v4]);
				assert(pth.vvs_[v2] < pth.vvs_[v3]);
#endif

				pth.dist_ += calc_dist(cities[v1], cities[v3]) + calc_dist(cities[v2], cities[v4]) -
						calc_dist(cities[v1], cities[v2]) - calc_dist(cities[v3], cities[v4]);

				pth.qxes_.del_edge(make_pair(cities[v1], cities[v2]));
				pth.qxes_.del_edge(make_pair(cities[v3], cities[v4]));
				pth.qxes_.add_edge(make_pair(cities[v1], cities[v3]), make_pair(v1, v3));
				pth.qxes_.add_edge(make_pair(cities[v2], cities[v4]), make_pair(v2, v4));

				assert(blacklist.erase(make_edge(v1, v2)) == 1);
				assert(blacklist.erase(make_edge(v2, v1)) == 1);
				assert(blacklist.erase(make_edge(v3, v4)) == 1);
				assert(blacklist.erase(make_edge(v4, v3)) == 1);
				blacklist.insert(make_edge(v1, v3));
				blacklist.insert(make_edge(v3, v1));
				blacklist.insert(make_edge(v2, v4));
				blacklist.insert(make_edge(v4, v2));
			
				pth.es_[v1].remove(v2);
				pth.es_[v2].remove(v1);
				pth.es_[v3].remove(v4);
				pth.es_[v4].remove(v3);
#ifndef NDEBUG
				/*
				if (pth.es_[v1].size() > 1 || pth.es_[v2].size() > 1 || pth.es_[v3].size() > 1 || pth.es_[v4].size() > 1)
				{
					cout << v1 << " ";
					for (list<vertex_id>::iterator it1 = pth.es_[v1].begin(); it1 != pth.es_[v1].end(); ++it1)
					{
						cout << (*it1) << " ";
					}
					cout << endl;
					cout << v2 << " ";
					for (list<vertex_id>::iterator it1 = pth.es_[v2].begin(); it1 != pth.es_[v2].end(); ++it1)
					{
						cout << (*it1) << " ";
					}
					cout << endl;
					cout << v3 << " ";
					for (list<vertex_id>::iterator it1 = pth.es_[v3].begin(); it1 != pth.es_[v3].end(); ++it1)
					{
						cout << (*it1) << " ";
					}
					cout << endl;
					cout << v4 << " ";
					for (list<vertex_id>::iterator it1 = pth.es_[v4].begin(); it1 != pth.es_[v4].end(); ++it1)
					{
						cout << (*it1) << " ";
					}
					cout << endl;
				}
				*/

				assert(pth.es_[v1].size() < 2);
				assert(pth.es_[v2].size() < 2);
				assert(pth.es_[v3].size() < 2);
				assert(pth.es_[v4].size() < 2);
#endif
				pth.es_[v1].push_back(v3);
				pth.es_[v3].push_back(v1);
				pth.es_[v2].push_back(v4);
				pth.es_[v4].push_back(v2);

				dq_node<vertex_id> *n1 = pth.vs_[v2];
				dq_node<vertex_id> *n2 = pth.vs_[v3];

				while (n1 != n2)
				{
					//cout << "XCHG " << n1 << " " << n2 << " N" << n1->d_ << " N" << n2->d_ << " P" << pth.vvs_[n1->d_] << " P" << pth.vvs_[n2->d_] << endl;
					swap(pth.vvs_[n1->d_], pth.vvs_[n2->d_]);
					swap(pth.vs_[n1->d_], pth.vs_[n2->d_]);
					swap(n1->d_, n2->d_);

					if (n1->next_ == n2)
					{
						break;
					}

					n1 = n1->next_;
					n2 = n2->prev_;
				}

				return true;
			}
		}
	}

	return false;
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

		if (! opt_x(s_kdt, *((*paths)[n]), blacklist, n))
		{
			return;

			if (! opt_x(s_kdt, *((*paths)[1 - n]), blacklist, 1 - n))
			{
				return;
			}
		}

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

