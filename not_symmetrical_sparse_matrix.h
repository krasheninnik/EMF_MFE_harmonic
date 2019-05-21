#include "pch.h"
#include <vector>

#define CALCULATE

class not_symmetrical_sparse_matrix {
	using LocalMatrix = std::vector<std::vector<double>>;
public:
	void init(const int elemNum, const int max_iter, const double eps);

	void load();
	void copy(not_symmetrical_sparse_matrix source);
	void addLocalMatrix(int elemNum, LocalMatrix M);
	int  get_dim();
	int  get_amount_elems();
	int  get_method();

	void setFirstBoundaryConditionsLeft();
	void setFirstBoundaryConditionsRight();
	void resetMatrix();

	void directSolveLU(std::vector<double>& l, std::vector<double>& ud,
		std::vector<double>& u, std::vector<double>& temp, std::vector<double>& q, std::vector<double>& f);

	void LU_partial_decomposite(std::vector<double>& l,
		std::vector<double>& ud, std::vector<double>& u);

	void diag_decomposite(std::vector<double>& diag);

	void forward_gauss(std::vector<double>& matrix, std::vector<double>& x,
		std::vector<double>& f);
	void forward_gauss(std::vector<double>& matrix, std::vector<double>& diag,
		std::vector<double>& x, std::vector<double>& f);

	void backward_gauss(std::vector<double>& matrix, std::vector<double>& x,
		std::vector<double>& f);
	void backward_gauss(std::vector<double>& matrix, std::vector<double>& diag,
		std::vector<double>& x, std::vector<double>& f);

	void diag_gauss(std::vector<double>& diag, std::vector<double>& x, std::vector<double>& f);

	void multiplicate_ufactor(std::vector<double>& u, std::vector<double>& diag,
		std::vector<double>& b, std::vector<double>& result);

	void LOS(std::vector<double>& x0, std::vector<double>& f,
		std::vector<double>& r,	std::vector<double>& z, 
		std::vector<double>& p, std::vector<double>& temp);

	void LOS_LU(std::vector<double>& l, std::vector<double>& ud, std::vector<double>& u,
		std::vector<double>& x0, std::vector<double>& f,
		std::vector<double>& r, std::vector<double>& z,
		std::vector<double>& p, std::vector<double>& temp,
		std::vector<double>& temp2);

	void LOS_DIAG(std::vector<double>& diag, std::vector<double>& x, std::vector<double>& f,
		std::vector<double>& r, std::vector<double>& z,
		std::vector<double>& p, std::vector<double>& temp,
		std::vector<double>& temp2);

	int MCG(std::vector<double>& x0, std::vector<double>& x_min,
		std::vector<double>& f, std::vector<double>& r, std::vector<double>& z,
		std::vector<double>& temp, std::vector<double> temp2);

	void MCG_LU(std::vector<double>& l, std::vector<double>& ud, std::vector<double>& u,
		std::vector<double>& x, std::vector<double>& x_min, std::vector<double>& f,
		std::vector<double>& r, std::vector<double>& z,
		std::vector<double>& temp, std::vector<double>& temp2);

	void MCG_DIAG(std::vector<double>& diag, std::vector<double>& x, std::vector<double>& x_min,
		std::vector<double>& f, std::vector<double>& r, std::vector<double>& z,
		std::vector<double>& temp, std::vector<double>& temp2);


	int BCGSTAB(std::vector<double> &l, std::vector<double> & ud, std::vector<double> & u,
		std::vector<double> & x, std::vector<double> & f, std::vector<double> & r, std::vector<double> &r_0,
		std::vector<double> & z, std::vector<double> & temp, std::vector<double> & temp2, std::vector<double> &v, std::vector<double> &p,
		std::vector<double> & y, std::vector<double> &h, std::vector<double> &s, std::vector<double> &t);

	void multiplicate_with_vector(std::vector<double>& b, std::vector<double>& result);
	void multiplicate_trans_with_vector(std::vector<double>& b, std::vector<double>& result);
	void multiplicate_diag_with_vector(std::vector<double>& b, std::vector<double>& result);

private:
	std::vector<int> ig;		// index elems - start of row
	std::vector<int> jg;		// column index of elem
	std::vector<double> di;		// diag elems
	std::vector<double> ggu;	// elem's of upper trinagle
	std::vector<double> ggl;	// elem's of lower trinagle

	int n;						// dimension
	int max_iter;				// max amount of iter
	double eps;					// min relative discrepancy
	const double step_eps = 1e-25;
	int method;
};