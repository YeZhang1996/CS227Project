#include "bn.h"
#include <exception>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

class badfactorexception : public exception {
public:
	badfactorexception(const string &reason) : w(reason) { }
	virtual ~badfactorexception() throw() { }
	virtual const char *what() const throw() { return w.c_str(); }
private:
	string w;
};
		
void bn::addvar(int v, const factor &f) {
	factor::scope s = f.getscope();
	if (s.find(v) == s.end()) {
		stringstream ss;
		ss << v;
		string reason;
		ss >> reason;
		reason = string("adding variable ")+reason+" with factor that does not include it";
		throw new badfactorexception(reason);
	}
	map<int,factor>::iterator i = factors.lower_bound(v);
	if (i!=factors.end() && i->first==v) {
		stringstream ss;
		ss << v;
		string reason;
		ss >> reason;
		reason = string("adding variable ")+reason+" to BN that already includes that variable";
		throw new badfactorexception(reason);
	}
	factors.insert(i,make_pair(v,f));
}

static inline bool cmpfirst(pair<int,int> p1, pair<int,int> p2) {
     return p1.first < p2.first;
}

bn::ctpotentials bn::calibrate(const bn::cliquetree &ct,
			const map<int,int> &alpha) const {
/*---------------------------Initialize beta and mu---------------------------*/

					bn::ctpotentials initpot;
					//initialize beta
					int i = 0;
					for(vector<factor::scope>::const_iterator ni = ct.nodes.begin(); ni != ct.nodes.end(); ni++){
						initpot.beta.push_back(factor(*ni,1));
						for(map<int,int>::const_iterator ai = alpha.begin(); ai != alpha.end(); ai++){
							if (ai->second == i){
								initpot.beta[i] = initpot.beta[i]*factors.at(ai->first);
							}
						}
						i++;
					}

					//init mu
					int j = 0;
					for(multimap<int,int>::const_iterator ei = ct.adj.begin(); ei != ct.adj.end(); ei++){
						i = ei->first;
						j = ei->second;
						if(i < j){
							factor::scope sij;
							set_intersection(ct.nodes[i].begin(),ct.nodes[i].end(),ct.nodes[j].begin(),ct.nodes[j].end(),
											  inserter(sij,sij.begin()),cmpfirst);
							initpot.mu.insert(make_pair(make_pair(i,j),factor(sij,1.0)));
						}
					}

/*----------------------------Initialization done------------------------------*/


/*----------------------------Upward belif propagation-------------------------*/

                //pick the last clique as root
                int root = ct.nodes.size()-1;
                //this map stores the number of neighbors for cliques. <Ci, # of neighbors of Ci>
                map<int,int> neighbors;
                for(i = 0; i <= root; i++){
                    neighbors[i] = ct.adj.count(i);
                }
                multimap<int,int> unmarked_edges = ct.adj;
                vector< pair<int,int> > reverse_order;
                for(map<int,int>::iterator ni = neighbors.begin(); ni != neighbors.end(); ni++){
                    if(ni->second == 1){
                        i = ni->first;
                        j = unmarked_edges.find(i)->second;
                        pair<int,int> i_j = make_pair(i,j);
                        if(i>j){
                            i_j.first = j;
                            i_j.second = i;
                        }
                        //mark i->j and j->i, meanwhile push j->i into reverse order
                        unmarked_edges.erase(i);
                        for(multimap<int,int>::iterator ei = unmarked_edges.begin();
                            ei != unmarked_edges.end(); ei++){
                            if(ei->second == i){
                                reverse_order.push_back(*ei);
                                unmarked_edges.erase(ei);
                            }
                        }
                        //update beta and mu
                        factor::scope tosumout; 
                        factor::scope sij = (initpot.mu.find(i_j))->second.getscope();
                        factor::scope ci = initpot.beta[i].getscope();
                        set_difference(ci.begin(), ci.end(), sij.begin(), sij.end(),
                                        inserter(tosumout,tosumout.begin()), cmpfirst);
                        factor sigmaij(sij);
                        sigmaij = initpot.beta[i].marginalize(tosumout);
                        initpot.beta[j] = initpot.beta[j]*sigmaij/initpot.mu.find(i_j)->second;
                        initpot.mu.find(i_j)->second = sigmaij;
                        //update finished
                        
                        neighbors.erase(ni);
                        neighbors.find(j)->second --;
                    }
                }

/*----------------------------Upward done-------------------------------------*/

/*----------------------------Downward belief propagation---------------------*/

                for(vector<pair<int,int> >::reverse_iterator ri = reverse_order.rbegin(); ri != reverse_order.rend(); ri++){
                    i = ri->first;
                    j = ri->second;
                    pair<int,int> i_j = make_pair(i,j);
                    if(i>j){
                        i_j.first = j;
                        i_j.second = i;
                    }
                    //update beta and mu
                    factor::scope tosumout; 
                    factor::scope sij = (initpot.mu.find(i_j))->second.getscope();
                    factor::scope ci = initpot.beta[i].getscope();
                    set_difference(ci.begin(), ci.end(), sij.begin(), sij.end(),
                                    inserter(tosumout,tosumout.begin()), cmpfirst);
                    factor sigmaij(sij);
                    sigmaij = initpot.beta[i].marginalize(tosumout);
                    initpot.beta[j] = initpot.beta[j]*sigmaij/initpot.mu.find(i_j)->second;
                    initpot.mu.find(i_j)->second = sigmaij;
                }
                    
/*----------------------------Downward done-------------------------------------*/
            //}
	return initpot;				
}

