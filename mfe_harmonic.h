#pragma once
#include "not_symmetrical_sparse_matrix.h"
#include <vector>
#include <functional>
#include <fstream>	

struct FiniteElem {
	int left;
	int right;
};

class HarmonicTask {
	using Matrix = not_symmetrical_sparse_matrix;
	using LocalMatrix = std::vector<std::vector<double>>;
	using func = std::function<double(const double&)>;
	using func2 = std::function<double(const double&, const double&)>;
	using func5 = std::function<double(const double, const double, const double, const double, const double)>;
	using func3d1fun = std::function<double(const double&, const double&, const double&, const func2&)>;

public:
	void init();
	void setParams();
	void resetTask();

	void setFirstBoundaryConditions();
	void solve();
	void explore();

	void saveResult(int);


	void calculateError(double & sinErr, double & cosErr);

	void saveExploreResult(const double time, const double sinErr, const double cosErr);
	void saveExploreResult(const int iter, const double time, const double sinErr, const double cosErr);

private:
	void calculateGlobalMatrixAndRightPart();

	void calculateLocalMatrix(uint32_t elemNum);
	void calculateLocalRightPart(uint32_t elemNum);

	void addLocalRigtPartToGlobal(uint32_t elemNum);

private:
	// result output stream
	std::ofstream fout;
	std::ofstream exploreout;
	std::vector<double> qExact;  // for calculate norm of error

	// space grid:
	std::vector<int> subareas;
	std::vector<FiniteElem> elems;
	std::vector<double> nodes;	// vector of Nodes

	// solutions:
	std::vector<double> q;
	//std::vector<double> temp;

	// local matrix of mass and rigidity, and local vector:
	LocalMatrix localMatrix;
	std::vector<double> fLocal;

	// global matrix and vector of right part:
	std::vector<double> f;
	Matrix globalMatrix;

	// vectors for solve:
	std::vector<double> x0; // q0
	std::vector<double> x_min; // q0
	std::vector<double> r;
	std::vector<double> r_0;
	std::vector<double>	z;
	std::vector<double>	p;
	std::vector<double> temp;
	std::vector<double> temp2;

	std::vector<double> ud;
	std::vector<double> l;
	std::vector<double>	u;

	std::vector<double> v;
	std::vector<double> y; 
	std::vector<double> h; 
	std::vector<double> s;
	std::vector<double> t;

	// parameters of equation in different subareas:
	int amountSubareas;

	func5 fFuncSin;
	func5 fFuncCos;
	func uSinExact;
	func uCosExact;

	std::vector<double> lambda;
	std::vector<double> sigma;
	std::vector<double> epsilon;

	double w;	// frequency of outer source

	// first boundary conditions
	double leftUsin;
	double leftUcos;
	double rightUsin;
	double rightUcos;
};