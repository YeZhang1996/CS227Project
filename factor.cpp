#include "factor.h"
#include <assert.h>
#include <string.h>
#include <iomanip>
#include <algorithm>

using namespace std;

static int scopesize(const factor::scope &s) {
	int ret = 1;
	for(factor::scope::const_iterator i=s.begin();i!=s.end();++i)
		ret *= i->second;
	return ret;
}

factor::factor(const scope &sc, double d)
		: s(sc), f(scopesize(s),d), stride(s.size()) {
	sz = 1;
	int ii=0;
	for(scope::const_iterator i=s.begin();i!=s.end();++i) {
		stride[ii++] = sz;
		sz *= i->second;
	}
}

double factor::operator()(const index &i) const {
	return f[i.i];
}

double &factor::operator()(const index &i) {
	return f[i.i];
}

int factor::assign2index(const factor::assign &i) const {
	int ii=0;
	int ind = 0;
	for(scope::const_iterator si=s.begin();si!=s.end();++si,++ii) {
		assign::const_iterator ai = i.find(si->first);
		if (ai!=i.end()) ind += ai->second*stride[ii];
	}
	return ind;
}

double factor::operator()(const factor::assign &i) const {
	return f[assign2index(i)];
}

double &factor::operator()(const factor::assign &i) {
	return f[assign2index(i)];
}

void factor::print(ostream &os) const {
	for(scope::const_iterator si=s.begin();si!=s.end();si++)
		os << "X" << setw(3) << left << si->first << ' ';
	os << "   value" << endl;
	for(scope::const_iterator si=s.begin();si!=s.end();si++)
		os << "-----";
	os << "----------------" << endl;
	for(index i(*this);i.i!=sz;++i) {
		for (int j=0;j<i.a.size();++j)
			os << ' ' << setw(3) << right << i.a[j] << ' ';
		os << "    " << *i << endl;
	}
}

static inline bool cmpfirst(const pair<int,int> &p1, const pair<int,int> &p2) {
	return p1.first < p2.first;
}

factor factor::reduce(const assign &a) const {
	scope news;
	set_difference(s.begin(),s.end(),a.begin(),a.end(),
			inserter(news,news.begin()), cmpfirst);
	factor ret(news);
	int newi = 0;
	for(index i(*this,a);i.i!=sz;i.incexcept(a),++newi)
		ret.f[newi] = *i;
	return ret;
}

factor factor::marginalize(const scope &tosumout) const {
	scope news;
	set_difference(s.begin(),s.end(),tosumout.begin(),tosumout.end(),
			inserter(news,news.begin()), cmpfirst);
	factor ret(news);
	for(index i(*this),j(ret);i.i!=sz;i.incalso(j))
		ret.f[j.i] += *i;
	return ret;
}

// C++98 does not support templates parameterized by operators:
// (well, it does, but the amount of extra verbiage necessary is silly in 
//  this case.)
#define GENERALOP(OP) \
	scope news; \
	set_union(s.begin(),s.end(),f2.s.begin(),f2.s.end(), \
			inserter(news,news.begin()),cmpfirst); \
	factor ret(news); \
	for(index i(ret),j1(*this),j2(f2);i.i!=ret.sz;i.incalso(j1,j2)) \
		ret.f[i.i] = f[j1.i] OP f2.f[j2.i]; \
	return ret;

factor factor::operator*(const factor &f2) const { GENERALOP(*); }
factor factor::operator/(const factor &f2) const { GENERALOP(/); }
factor factor::operator+(const factor &f2) const { GENERALOP(+); }
factor factor::operator-(const factor &f2) const { GENERALOP(-); }

factor &factor::operator*=(double d) {
	for(int i=0;i<sz;i++) f[i] *= d;
	return *this;
}
factor &factor::operator/=(double d) {
	for(int i=0;i<sz;i++) f[i] /= d;
	return *this;
}
factor &factor::operator+=(double d) {
	for(int i=0;i<sz;i++) f[i] += d;
	return *this;
}
factor &factor::operator-=(double d) {
	for(int i=0;i<sz;i++) f[i] -= d;
	return *this;
}

//-----------------------------------------------------
factor::index::index(const factor &fact) : f(fact), a(fact.s.size(),0) {
     i = 0;
}

factor::index::index(const factor &fact, const assign &a)
			: f(fact), a(fact.s.size()) {
	i = 0;
	int ii=0;
	assign::const_iterator ai = a.begin();
	for(scope::const_iterator si=f.s.begin();si!=f.s.end();++si,++ii) {
		while (ai != a.end() && ai->first < si->first)
			++ai;
		if (ai != a.end() && ai->first == si->first) {
			i += ai->second*f.stride[ii];
			this->a[ii] = ai->second;
		} else {
			this->a[ii] = 0;
		}
	}
}


factor::index &factor::index::operator++() {
	int ii=0;
	for(scope::const_iterator si=f.s.begin();si!=f.s.end();++si,++ii)
		if (incelem(ii,si->second-1)) return *this;
	i = f.sz; // to indicate we've finished
	return *this;
}

const factor::index factor::index::operator++(int) {
	index ret(*this);
	++(*this);
	return ret;
}

double factor::index::operator*() const  {
	return f.f[i];
}

bool factor::index::incelem(int ii, int maxs) {
	if (a[ii] == maxs) {
		a[ii] = 0;
		i -= maxs*f.stride[ii];
		return false;
	} else {
		a[ii]++;
		i += f.stride[ii];
		return true;
	}
}

void factor::index::incexcept(const assign &except) {
	int ii=0;
	assign::const_iterator ei = except.begin();
	for(scope::const_iterator si=f.s.begin();si!=f.s.end();++si,++ii) {
		while (ei != except.end() && ei->first < si->first)
			++ei;
		if (ei != except.end() && ei->first == si->first) continue;
		if (incelem(ii,si->second-1)) return;
	}
	i = f.sz; // to indicate we've finished
}

void factor::index::incalso(index &i2) {
	int ii1=0,ii2=0;
	scope::const_iterator s2i = i2.f.s.begin();
	for(scope::const_iterator s1i = f.s.begin();s1i != f.s.end();++s1i,++ii1) {
		while(s2i != i2.f.s.end() && s2i->first < s1i->first)
			++s2i, ++ii2;
		if (s2i != i2.f.s.end() && s2i->first == s1i->first)
			i2.incelem(ii2,s1i->second-1);
		if (incelem(ii1,s1i->second-1)) return;
	}
	i = f.sz;
	i2.i = i2.f.sz;
}

void factor::index::incalso(index &i2, index &i3) {
	int ii1=0,ii2=0,ii3=0;
	scope::const_iterator s2i = i2.f.s.begin();
	scope::const_iterator s3i = i3.f.s.begin();
	for(scope::const_iterator s1i = f.s.begin();s1i != f.s.end();++s1i,++ii1) {
		while(s2i != i2.f.s.end() && s2i->first < s1i->first)
			++s2i, ++ii2;
		while(s3i != i3.f.s.end() && s3i->first < s1i->first)
			++s3i, ++ii3;
		if (s2i != i2.f.s.end() && s2i->first == s1i->first)
			i2.incelem(ii2,s1i->second-1);
		if (s3i != i3.f.s.end() && s3i->first == s1i->first)
			i3.incelem(ii3,s1i->second-1);
		if (incelem(ii1,s1i->second-1)) return;
	}
	i = f.sz;
	i2.i = i2.f.sz;
	i3.i = i3.f.sz;
}
