#ifndef FACTOR_H
#define FACTOR_H

#include <iostream>
#include <map>
#include <cstdarg>
#include <vector>

class factor {
public:
	typedef std::map<int,int> scope;
	typedef std::map<int,int> assign;

	// build blank factors (all values=d) over scope 
	factor(const scope &sc, double d=0.0);


	double operator()(const assign &i) const;
	double &operator()(const assign &i);


	void print(std::ostream &os) const;

	factor reduce(const assign &a) const;
	factor marginalize(const scope &tosumout) const;
	factor operator*(const factor &f2) const;
	factor operator/(const factor &f2) const;
	factor operator+(const factor &f2) const;
	factor operator-(const factor &f2) const;

	factor &operator*=(double d);
	factor &operator/=(double d);
	factor &operator+=(double d);
	factor &operator-=(double d);

	scope getscope() const { return s; }
private:
	struct index {
		index(const factor &f);
		index(const factor &f, const assign &a);

		index &operator++();
		const index operator++(int);
		void incexcept(const assign &except);
		void incalso(index &i2);
		void incalso(index &i2, index &i3);
		double operator*() const;

		int i;
		const factor &f;
		std::vector<int> a;

		bool incelem(int ii, int maxs);
	};

	// index must have been created with this factor
	double operator()(const index &i) const;
	double &operator()(const index &i);

	int assign2index(const assign &a) const;

	scope s;
	int sz;
	std::vector<double> f;
	std::vector<int> stride;
};

#endif
