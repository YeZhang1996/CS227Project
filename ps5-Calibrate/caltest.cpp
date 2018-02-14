#include "bn.h"
#include <cstdarg>
#include <algorithm>

using namespace std;

// handy way to set an assignment quickly
// (variable arguments at the end are all integers: the assignments "in order")
factor::assign A(const factor::scope s, ...) {
     factor::assign ret;
     va_list as;
     va_start(as,s);
     for(factor::scope::const_iterator si=s.begin();si!=s.end();si++) {
          int v = va_arg(as,int);
          ret[si->first] = v;
     }
     va_end(as);
     return ret;
}

static inline bool cmpfirst(pair<int,int> p1, pair<int,int> p2) {
     return p1.first < p2.first;
}

// scope union:
factor::scope operator+(const factor::scope &s1, const factor::scope &s2) {
	factor::scope ret;
	set_union(s1.begin(),s1.end(),s2.begin(),s2.end(), 
               inserter(ret,ret.begin()),cmpfirst); 
	return ret;
}


// Make the student example from the assignment (not quite the same as
//  the one in the book)
void makeex(bn &n, bn::cliquetree &ct, map<int,int> &alpha) {
	// variable mapping C=0, D=1, I=2, T=3, G=4, S=5, L=6, J=7
	enum {C=0,D=1,I=2,T=3,G=4,S=5,L=6,J=7};

	factor::scope cvar; cvar[C] = 2;
	factor::scope dvar; dvar[D] = 2;
	factor::scope ivar; ivar[I] = 2;
	factor::scope tvar; tvar[T] = 2;
	factor::scope gvar; gvar[G] = 3;
	factor::scope svar; svar[S] = 2;
	factor::scope lvar; lvar[L] = 2;
	factor::scope jvar; jvar[J] = 2;

	factor::scope c = cvar;
	factor::scope d = cvar+dvar;
	factor::scope i = ivar;
	factor::scope t = ivar+tvar;
	factor::scope g = dvar+tvar+gvar;
	factor::scope s = tvar+svar;
	factor::scope l = gvar+lvar;
	factor::scope j = svar+lvar+jvar;

	factor pc(c),pd(d),pi(i),pt(t),pg(g),ps(s),pl(l),pj(j);
	pc(A(c,0)) = 0.5; pc(A(c,1)) = 0.5;

	pd(A(d,0,0)) = 0.4; pd(A(d,0,1)) = 0.6;
	pd(A(d,1,0)) = 0.8; pd(A(d,1,1)) = 0.2;

	pi(A(i,0)) = 0.6; pi(A(i,1)) = 0.4;

	pt(A(t,0,0)) = 0.9; pt(A(t,0,1)) = 0.1;
	pt(A(t,1,0)) = 0.4; pt(A(t,1,1)) = 0.6;

	pg(A(g,0,0,0)) = 0.3; pg(A(g,0,0,1)) = 0.4; pg(A(g,0,0,2)) = 0.3;
	pg(A(g,1,0,0)) = .05; pg(A(g,1,0,1)) = .25; pg(A(g,1,0,2)) = 0.7;
	pg(A(g,0,1,0)) = 0.9; pg(A(g,0,1,1)) = .08; pg(A(g,0,1,2)) = .02;
	pg(A(g,1,1,0)) = 0.5; pg(A(g,1,1,1)) = 0.3; pg(A(g,1,1,2)) = 0.2;

	ps(A(s,0,0)) = 0.95; ps(A(s,0,1)) = 0.05;
	ps(A(s,1,0)) = 0.2;  ps(A(s,1,1)) = 0.8;

	pl(A(l,0,0)) = 0.1; pl(A(l,0,1)) = 0.9;
	pl(A(l,1,0)) = 0.4; pl(A(l,1,1)) = 0.6;
	pl(A(l,2,0)) = 0.99; pl(A(l,2,1)) = 0.01;

	pj(A(j,0,0,0)) = 0.9; pj(A(j,0,0,1)) = 0.1;
	pj(A(j,1,0,0)) = 0.4; pj(A(j,1,0,1)) = 0.6;
	pj(A(j,0,1,0)) = 0.3; pj(A(j,0,1,1)) = 0.7;
	pj(A(j,1,1,0)) = 0.1; pj(A(j,1,1,1)) = 0.9;

	n.addvar(C,pc);
	n.addvar(D,pd);
	n.addvar(I,pi);
	n.addvar(T,pt);
	n.addvar(G,pg);
	n.addvar(S,ps);
	n.addvar(L,pl);
	n.addvar(J,pj);

	// cliquetree:
	ct.nodes.push_back(cvar+dvar);    //0
	ct.nodes.push_back(gvar+dvar+tvar);  //1
	ct.nodes.push_back(ivar+tvar);    //2
	ct.nodes.push_back(gvar+tvar+svar);  //3
	ct.nodes.push_back(svar+lvar+gvar);  //4
	ct.nodes.push_back(jvar+lvar+svar);  //5

	ct.adj.insert(make_pair(0,1)); ct.adj.insert(make_pair(1,0));
	ct.adj.insert(make_pair(1,3)); ct.adj.insert(make_pair(3,1));
	ct.adj.insert(make_pair(2,3)); ct.adj.insert(make_pair(3,2));
	ct.adj.insert(make_pair(3,4)); ct.adj.insert(make_pair(4,3));
	ct.adj.insert(make_pair(4,5)); ct.adj.insert(make_pair(5,4));

	alpha[C] = 0;
	alpha[D] = 0;
	alpha[I] = 2;
	alpha[T] = 2;
	alpha[G] = 1;
	alpha[S] = 3;
	alpha[L] = 4;
	alpha[J] = 5;

}

void printmarginals(const factor &f) {
	factor::scope s = f.getscope();
	cout << "corresponding marginals:" << endl;
	for(factor::scope::iterator j=s.begin();j!=s.end();++j) {
		factor::scope s2 = s;
		s2.erase(j->first);
		(f.marginalize(s2)).print(cout);
	}
	cout << endl;
}

int main(int argc, char **argv) {
	bn n;
	bn::cliquetree ct;
	//alpha assign factors to cliques. <i,j> i is factor index in fmap, j is the clique index in nodes
	std::map<int,int> alpha;
	makeex(n,ct,alpha);

	bn::ctpotentials pot = n.calibrate(ct,alpha);

	// print all marginals from each factor
	for(vector<factor>::iterator i=pot.beta.begin();i!=pot.beta.end();++i) {
		cout << "joint from clique node:" << endl;
		i->print(cout);
		printmarginals(*i);
	}

	for(map<pair<int,int>,factor>::iterator i=pot.mu.begin();i!=pot.mu.end();++i) {
		if (i->first.first > i->first.second) continue;
		cout << "sepset joint:" << endl;
		(i->second).print(cout);
		printmarginals(i->second);
	}
}
