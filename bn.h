#ifndef BN_H
#define BN_H

#include "factor.h"
#include<map>
#include<set>
#include<vector>

class bn { // class representing a Bayesian network
public:
	// f needs to be a factor that includes v in its scope
	// (and represents the conditional distribution of v)
	void addvar(int v, const factor &f);

	// produces the marginal distribution over the variables
	// in remaining, as a factor
	// employs variable elimination with the elimination order
	// given by order:
	//   all variables in the parameter order are eliminated first
	//         in the order they appear in the vector order
	//         (if a variable in remaining appears in order, it
	//          is ignored -- as if it did not appear in order)
	//   then, if any remaining variables need to be eliminated
	//     they may be eliminated in any order
	factor varelim(const std::set<int> &remaining,
		const std::vector<int> &order = std::vector<int>()) const;

	// a clique tree
	class cliquetree {
	public:
		std::vector<factor::scope> nodes; // nodes are consecutively
			// numbered from 0.  nodes[i] is the set of variables
			// assigned to clique tree node i
		std::multimap<int,int> adj; // adjacency of cliques
			// for ease both i->j and j->i are included
	};

	// clique tree potentials
	class ctpotentials {
	public:
		std::vector<factor> beta; // beliefs for each node
		std::map<std::pair<int,int>, factor> mu; // beliefs for each sepset
		// (note that mu[<1,2>] and mu[<2,1>] should be the same
		//  but we don't want to store both -- only store mu[<i,j>]
		//  if i<j)
	};

	// takes in a cliquetree and returns the calibrated potentials
	// from it by running belief propagation
	// (may assume ct is family preserving and obeys the running
	//  intersection property wrt this bn)
	// alpha is the mapping from BN variable to the ct node which
	//   preserves the variable's CPD
	ctpotentials calibrate(const cliquetree &ct,
			const std::map<int,int> &alpha) const;
	
private:
	// invariant:
	// if factors maps i to f, then f has i in its scope
	//    and if f is marginalized over i, it becomes a factor of all 1s
	//    (that is, it is a conditional distribution over i)
	// Note, the members of the scope of f that are not i are the
	//   parents of i in the BN
	typedef std::map<int,factor> fmap;
	fmap factors;
};
	
#endif
