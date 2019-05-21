#include "pch.h"
#include "not_symmetrical_sparse_matrix.h"
#include <stdio.h>
#include <iostream>

#pragma region vector_operations
inline void vector_add(std::vector<double>& a, std::vector<double>& b, std::vector<double>& result);
inline void vector_subtract(std::vector<double>& a, std::vector<double>& b, std::vector<double>& result);
inline void vector_assigment(std::vector<double>& a, std::vector<double>& b);
inline void vector_mult_on_const(std::vector<double>& a, double constant, std::vector<double>& result);
inline double scalar_product(std::vector<double>& a, std::vector<double>& b);
inline double calc_norm(std::vector<double>& vector);
#pragma endregion

#pragma region  someUtilits
int not_symmetrical_sparse_matrix::get_dim() { return n; }
int not_symmetrical_sparse_matrix::get_amount_elems() { return ggl.size(); }
int not_symmetrical_sparse_matrix::get_method() { return method; }
void not_symmetrical_sparse_matrix::copy(not_symmetrical_sparse_matrix source) {
	ig = source.ig;				// index elems - start of row
	jg = source.jg;				// column index of elem
	di = source.di;				// diag elems
	ggu = source.ggu;			// elem's of upper trinagle
	ggl = source.ggl;			// elem's of lower trinagle

	n = source.n;						// dimension
	max_iter = source.max_iter;			// max amount of iter
	eps = source.eps;					// min relative discrepancy
	method = source.method;
}
void not_symmetrical_sparse_matrix::load() {
	FILE *in;
	fopen_s(&in, "kuslau.txt", "r");
	fscanf_s(in, "%d %d %le", &n, &max_iter, &eps);
	fclose(in);

	fopen_s(&in, "method.txt", "r");
	fscanf_s(in, "%d", &method);
	fclose(in);

	di.resize(n);
	ig.resize(n + 1);

	fopen_s(&in, "ig.txt", "r");
	for (int i = 0; i < n + 1; i++) fscanf_s(in, "%d", &ig[i]);
	fclose(in);

	fopen_s(&in, "di.txt", "r");
	for (int i = 0; i < n; i++) fscanf_s(in, "%le", &di[i]);
	fclose(in);

	jg.resize(ig[n]);
	ggu.resize(ig[n]);
	ggl.resize(ig[n]);

	fopen_s(&in, "jg.txt", "r");
	for (int i = 0; i < jg.size(); i++) {
		fscanf_s(in, "%d", &jg[i]);
	}
	fclose(in);

	fopen_s(&in, "ggu.txt", "r");
	for (int i = 0; i < ggu.size(); i++) fscanf_s(in, "%le", &ggu[i]);
	fclose(in);

	fopen_s(&in, "ggl.txt", "r");
	for (int i = 0; i < ggl.size(); i++) fscanf_s(in, "%le", &ggl[i]);
	fclose(in);

	if (ig[0] == 1) {	// transform FORTRAN format --> C Format.
		for (int i = 0; i < ig.size(); i++) ig[i]--;
		for (int i = 0; i < jg.size(); i++) jg[i]--;
	}
}
#pragma endregion

void not_symmetrical_sparse_matrix::resetMatrix() {
	for (auto &el : ggl) el = 0;
	for (auto &el : ggu) el = 0;
	for (auto &el : di) el = 0;
}

void not_symmetrical_sparse_matrix::init(const int elemNum, const int _max_iter, const double _eps) {
	// 6 elems of first finite elem
	// 5 elements for each next finite elem
	int amountElemsInMatrix = elemNum * 5 + 1; // 6 + (num - 1)*5 = num*5 + 1:

	// 4 dim for first finite elem
	// 2 dim for each next finite elem
	n = elemNum * 2 + 2;	// 4 + (num - 1)*2 = num*2 + 2:
	max_iter = _max_iter;
	eps = _eps;

	// memory allocation
	di = std::vector<double>(n);
	ig = std::vector<int>(n + 1);
	jg = std::vector<int>(amountElemsInMatrix);
	ggu = std::vector<double>(amountElemsInMatrix);
	ggl = std::vector<double>(amountElemsInMatrix);

	// forming ig array
	ig[0] = 0;
	ig[1] = 0;
	ig[2] = 1;
	for (int pos = 3; pos < n; pos += 2) {
		ig[pos] = ig[pos - 1] + 2;
		ig[pos + 1] = ig[pos] + 3;
	}

	// formin jg array
	jg[0] = 0;
	for (int pos = 1, i = 0, column = 0; pos < amountElemsInMatrix; pos++) {
		jg[pos] = column;				// columns:
		jg[++pos] = column + 1;			// j j+1
		jg[++pos] = column;				// j j+1 j+2
		jg[++pos] = ++column;
		jg[++pos] = ++column;
	}
}

void not_symmetrical_sparse_matrix::addLocalMatrix(int elemNum, LocalMatrix M) {
	// consider diag elems
	const int diagPos = elemNum * 2;
	di[diagPos] += M[0][0];
	di[diagPos + 1] += M[1][1];
	di[diagPos + 2] += M[2][2];
	di[diagPos + 3] += M[3][3];

	const int elemPos = elemNum * 5;
	ggl[elemPos] += M[1][0];
	ggu[elemPos] += M[0][1];

	ggl[elemPos + 1] += M[2][0];
	ggu[elemPos + 1] += M[0][2];

	ggl[elemPos + 2] += M[2][1];
	ggu[elemPos + 2] += M[1][2];

	ggl[elemPos + 3] += M[3][0];
	ggu[elemPos + 3] += M[0][3];

	ggl[elemPos + 4] += M[3][1];
	ggu[elemPos + 4] += M[1][3];

	ggl[elemPos + 5] += M[3][2];
	ggu[elemPos + 5] += M[2][3];
}

#pragma region multiplications

void not_symmetrical_sparse_matrix::multiplicate_with_vector(std::vector<double>& b,
	std::vector<double>& result) {
	for (int i = 0; i < result.size(); i++) result[i] = di[i] * b[i];	// init result vect

	for (int i = 0; i < ig.size() - 1; i++) {								// multiplicate:
		for (int j = ig[i]; j < ig[i + 1]; j++) {
			result[i] += ggl[j] * b[jg[j]];
			result[jg[j]] += ggu[j] * b[i];
		}
	}
}

void not_symmetrical_sparse_matrix::multiplicate_trans_with_vector(std::vector<double>& b,
	std::vector<double>& result) {
	for (int i = 0; i < result.size(); i++) result[i] = di[i] * b[i];		// init result vect

	for (int i = 0; i < ig.size() - 1; i++) {								// multiplicate:
		for (int j = ig[i]; j < ig[i + 1]; j++) {
			result[i] += ggu[j] * b[jg[j]];
			result[jg[j]] += ggl[j] * b[i];
		}
	}
}

void not_symmetrical_sparse_matrix::multiplicate_diag_with_vector(
			std::vector<double>& b, std::vector<double>& result) {
	for (int i = 0; i < di.size(); i++) result[i] = di[i] * b[i];
}

void  not_symmetrical_sparse_matrix::
		multiplicate_ufactor(std::vector<double>& u, std::vector<double>& diag,
		std::vector<double>& b, std::vector<double>& result) {

	for(int i = 0; i < result.size(); i++) result[i] = diag[i] * b[i];	// diag elems

	for (int i = 0; i < ig.size() - 1; i++) {		// multiplicate:		
		for (int j = ig[i]; j < ig[i + 1]; j++) {	// U elems
			result[jg[j]] += u[j] * b[i];
		}
	}
}

#pragma endregion

#pragma region solve_methods

void not_symmetrical_sparse_matrix::directSolveLU(std::vector<double>& l, std::vector<double>& ud,
	std::vector<double>& u, std::vector<double>& temp, std::vector<double>& q, std::vector<double>& f) {
	// note:
	// this method work only if LU partial decompose equalent to exact LU decompose!
	LU_partial_decomposite(l, ud, u);
	forward_gauss(l, temp, f);
	backward_gauss(u, ud, q, temp);
}


void not_symmetrical_sparse_matrix::forward_gauss(std::vector<double>& matrix,
		std::vector<double>& x,	std::vector<double>& f) {
	int i, j;
	double sum;

	//x[i] = f[i] - sum[0;i-1](Lki * Xi);

	for (i = 0; i < n; i++) {
		sum = 0;
		for (j = ig[i]; j < ig[i + 1]; j++) {
			sum += matrix[j] * x[jg[j]];	// 
		}
		x[i] = f[i] - sum;
	}
}

void not_symmetrical_sparse_matrix::forward_gauss(std::vector<double>& matrix,
	std::vector<double>& diag, std::vector<double>& x, std::vector<double>& f) {
	int i, j;
	double sum;

	//x[i] = f[i] - sum[0;i-1](Lki * Xi);

	for (i = 0; i < n; i++) {
		sum = 0;
		for (j = ig[i]; j < ig[i + 1]; j++) {
			sum += matrix[j] * x[jg[j]];	// 
		}
		x[i] = (f[i] - sum) / diag[i];
	}
}

void not_symmetrical_sparse_matrix::backward_gauss(std::vector<double>& matrix,
	std::vector<double>& x, std::vector<double>& f) {

	// method will change vector f
	for (int i = n - 1; i >= 0; i--) {
		x[i] = f[i];

		for (int j = ig[i]; j < ig[i + 1]; j++) {
			f[jg[j]] -= matrix[j] * x[i];
		}
	}
}

void not_symmetrical_sparse_matrix::backward_gauss(std::vector<double>& matrix,
	std::vector<double>& diag,	std::vector<double>& x, std::vector<double>& f) {
	// method will change vector f
	for (int i = n - 1; i >= 0; i--) {
		x[i] = f[i] / diag[i];

		for (int j = ig[i]; j < ig[i + 1]; j++) {
			f[jg[j]] -= matrix[j] * x[i];
		}
	}
}

void not_symmetrical_sparse_matrix::diag_gauss(std::vector<double>& diag, std::vector<double>& x,
	std::vector<double>& f) {
	for (int i = 0; i < diag.size(); i++) x[i] = f[i] / diag[i];
}

void not_symmetrical_sparse_matrix::diag_decomposite(std::vector<double>& diag) {
	for (int i = 0; i < diag.size(); i++) {
		if (diag[i] >= 0)
			diag[i] = sqrt(di[i]);
		else {
			printf("error, diag elem < 0");
			return;
		}
	}
}

void not_symmetrical_sparse_matrix::LU_partial_decomposite(std::vector<double>& l,
			std::vector<double>& ud, std::vector<double>& u) {

	int column, row;
	int to_column, to_row;
	int first_column, first_row;
	int amount_column, amount_row;

	double diag_sum, l_sum, u_sum;
	int i, j, k, t;

	for (i = 0; i < n; i++) {
		diag_sum = 0;
	
		first_column = ig[i];
		//to_column = 0;

		for (j = first_column, amount_column = 0; j < ig[i + 1]; j++, amount_column++) {
			l_sum = 0;
			u_sum = 0;

			column = jg[j];
			amount_row = ig[column + 1] - ig[column];
			first_row = ig[column];

			// find elems to multiplicate (row == column)
			to_column = to_row = 0;

			// while don't yet reviewed all the columns and rows and
			//	current row < max column
			while	(to_column < amount_column && to_row < amount_row &&
						jg[first_row + to_row] < jg[j]) {	// j = first_column + amount_column
				
				if (jg[first_row + to_row] < jg[first_column + to_column]) to_row++;
				if (jg[first_row + to_row] > jg[first_column + to_column]) to_column++;
				if (jg[first_row + to_row] == jg[first_column + to_column]) {
					l_sum += l[first_column + to_column] * u[first_row + to_row];			// accumulate sum for calc L elems
					u_sum += u[first_column + to_column] * l[first_row + to_row];			// accumulate sum for calc U elems

					to_column++;
					to_row++;
				}
			}

			//printf("calc l elem: column: %d\n", column);
			l[j] = (ggl[j] - l_sum) / ud[column];		// L elem
			u[j] =  ggu[j] - u_sum;						// U elem
			
			diag_sum += l[j] * u[j];
		}

		// Chech diag elem:

		if (abs(ud[i] - di[i] + diag_sum) < ud[i] * 1e-15) {
			printf_s("error in LDU decomposition. zero on diag.\n");
			return;
		}
		else {
			ud[i] = di[i] - diag_sum;						// Udiag elem
		}
	}
}

#include <time.h>
void not_symmetrical_sparse_matrix::LOS(std::vector<double>& x,
	std::vector<double>& f, std::vector<double>& r, std::vector<double>& z,
	std::vector<double>& p, std::vector<double>& temp) {
	

	double a, b, scalar_pp, scalar_rr, discrepancy, norm_f;
	int iter = 0;

	// r0 = f - Ax0:
	multiplicate_with_vector(x, temp);		// Ax0
	vector_subtract(f, temp, r);			// r0 = f - Ax0

	z = r;									// z0 = r0
	multiplicate_with_vector(z, p);			// p0 = Az0	

	norm_f = calc_norm(f);
	scalar_rr = scalar_product(r, r);

	double CalcTimeStart = clock() / (double)CLOCKS_PER_SEC;
	// iterations:
	do {
		if (iter > max_iter) break;

		//scalar_rr = scalar_product(r, r);
		scalar_pp = scalar_product(p, p);		// scalar_pp = (pk-1,pk-1)
		a = scalar_product(p, r) / scalar_pp; // ak = (pk-1, rk-1)/(pk-1,pk-1)
		
		vector_mult_on_const(z, a, temp);		// temp = a*zk-1
		vector_add(x, temp, x);				// x0 = x0 + a*zk-1

		// step out:
		if (calc_norm(temp) / calc_norm(x) < step_eps) {
#ifdef CALCULATE
			printf("\nstep out\n");
#endif
			break;
		}

		vector_mult_on_const(p, a, temp);		// temp = a*pk-1
		vector_subtract(r, temp, r);			//  rk = rk-1 - a*pk-1

		multiplicate_with_vector(r, temp);		// temp = Ark
		b = -1.0 * scalar_product(p, temp) / scalar_pp;	// bk = -(pk-1, Ark)/(pk-1,pk-1)

		vector_mult_on_const(z, b, temp);		// temp = b*zk-1
		vector_add(r, temp, z);				// zk  = zk-1  + b*zk-1

		// pk = Ark + Bkpk-1 <==> 1. pk = Bkpk-1; 2. pk += Ark
		vector_mult_on_const(p, b, p);
		multiplicate_with_vector(r, temp);		// temp = Ark
		vector_add(p, temp, p);


		// check discrepancy
		scalar_rr -= a * a * scalar_pp;
		if (scalar_rr < 0) scalar_rr = scalar_product(r, r);

		discrepancy = sqrt(scalar_rr) / norm_f;
	
		iter++;
#ifdef _DEBUG
		printf("\riter: %d\t %.15e", iter, discrepancy);

		



#endif
	} while (discrepancy > eps);
	double CalcTimeStop = clock() / (double)CLOCKS_PER_SEC;
	//printf("\nEx Time: %f sec\n", (CalcTimeStop - CalcTimeStart));
}

void not_symmetrical_sparse_matrix::
	LOS_LU(std::vector<double>& l, std::vector<double>& ud, std::vector<double>& u,
	std::vector<double>& x, std::vector<double>& f,
	std::vector<double>& r, std::vector<double>& z,
	std::vector<double>& p, std::vector<double>& temp,
	std::vector<double>& temp2) {

	int iter = 0;
	double scalar_pp, scalar_rr, discrepancy, a, b, norm_f;

	// partial LU decomposition	
	LU_partial_decomposite(l, ud, u);

	//r = ..
	multiplicate_with_vector(x, temp);
	vector_subtract(f, temp, temp2);	
	forward_gauss(l, r, temp2);

	// z = ..
	vector_assigment(temp, r);
	backward_gauss(u, ud, z, temp);

	// p = ...
	multiplicate_with_vector(z, temp);
	forward_gauss(l, p, temp);

	norm_f = calc_norm(f);

	scalar_rr = scalar_product(r, r);
	do {
		if (iter > max_iter) break;
		
		scalar_pp = scalar_product(p, p);		// scalar_pp = (pk-1,pk-1)
		
		a = scalar_product(p, r) / scalar_pp;

		// x = ...
		vector_mult_on_const(z, a, temp);
		vector_add(x, temp, x);

		// step out:
		if (calc_norm(temp) / calc_norm(x) < step_eps) {
#ifdef CALCULATE
			printf("\nstep out\n");
#endif
			break;
		}

		//r = ...
		vector_mult_on_const(p, a, temp);
		vector_subtract(r, temp, r);

		// b = ...
		vector_assigment(temp, r);
		backward_gauss(u, ud, temp2, temp);
		multiplicate_with_vector(temp2, temp);
		forward_gauss(l, temp2, temp);

		b = scalar_product(p, temp2) / scalar_pp;

		// p = ...
		vector_mult_on_const(p, b, temp);
		vector_add(temp2, temp, p);

		// z = ...
		vector_assigment(temp, r);
		backward_gauss(u, ud, temp2, temp);
		vector_mult_on_const(z, b, temp);
		vector_add(temp, temp2, z);

		// check discrepancy
		scalar_rr -= a * a * scalar_pp;
		//scalar_rr = scalar_product(r, r);

		if (scalar_rr < 0) scalar_rr = scalar_product(r, r);
		discrepancy = sqrt(scalar_rr) / norm_f;
		//discrepancy = scalar_rr / norm_f;

		//for (int i = 0; i < x0.size(); i++) printf("%.15lf\n", x0[i]);

		iter++;
#ifdef CALCULATE
		printf("\riter: %d\t discr: %.15e", iter, discrepancy);
#endif
	} while (discrepancy > eps);
}


void not_symmetrical_sparse_matrix::
	LOS_DIAG(std::vector<double>& diag, std::vector<double>& x, std::vector<double>& f,
	std::vector<double>& r, std::vector<double>& z,
	std::vector<double>& p, std::vector<double>& temp,
	std::vector<double>& temp2) {

	int iter = 0;
	double scalar_pp, scalar_rr, discrepancy, a, b, norm_f;

	// L = U = sqrt(di)
	diag_decomposite(diag);

	//r = ..
	multiplicate_with_vector(x, temp);
	vector_subtract(f, temp, temp2);	

	diag_gauss(diag, r, temp2);

	// z = ..
	vector_assigment(temp, r);
	diag_gauss(diag, z, temp);

	// p = ...
	multiplicate_with_vector(z, temp);
	diag_gauss(diag, p, temp);

	norm_f = calc_norm(f);
	scalar_rr = scalar_product(r, r);
	
	double CalcTimeStart = clock() / (double)CLOCKS_PER_SEC;
	do {
		if (iter > max_iter) break;
		
		scalar_pp = scalar_product(p, p);		// scalar_pp = (pk-1,pk-1)
		a = scalar_product(p, r) / scalar_pp;

		// x = ...
		vector_mult_on_const(z, a, temp);
		vector_add(x, temp, x);

		// step out:
		if (calc_norm(temp) / calc_norm(x) < step_eps) {
#ifdef CALCULATE
			printf("\nstep out\n");
#endif
			break;
		}


		//r = ...
		vector_mult_on_const(p, a, temp);
		vector_subtract(r, temp, r);

		// b = ...
		vector_assigment(temp, r);
		diag_gauss(diag, temp2, temp);
		multiplicate_with_vector(temp2, temp);
		diag_gauss(diag, temp2, temp);

		b = scalar_product(p, temp2) / scalar_pp;

		// p = ...
		vector_mult_on_const(p, b, temp);
		vector_add(temp2, temp, p);

		// z = ...
		vector_assigment(temp, r);
		diag_gauss(diag, temp2, temp);
		vector_mult_on_const(z, b, temp);
		vector_add(temp, temp2, z);

		// check discrepancy
		scalar_rr -= a * a * scalar_pp;

		if (scalar_rr < 0) printf("!"), scalar_rr = scalar_product(r, r);
		discrepancy = sqrt(scalar_rr) / norm_f;

		//discrepancy = scalar_rr / norm_f;


		iter++;
#ifdef CALCULATE
		printf("\riter: %d\t discr: %.15e", iter, discrepancy);
#endif
	} while (discrepancy > eps);

	double CalcTimeStop = clock() / (double)CLOCKS_PER_SEC;
	printf("\nEx Time: %f sec\n", (CalcTimeStop - CalcTimeStart));
}


int not_symmetrical_sparse_matrix::MCG(std::vector<double>& x, std::vector<double>& x_min,
	std::vector<double>& f, std::vector<double>& r, std::vector<double>& z,
	std::vector<double>& temp, std::vector<double> temp2) {

	double a, b, scalar_rr, scalar_rr_prev, norm_f, discrepancy, min_discrepancy;
	int iter = 0;
												// system symmetrization:
	//vector_assigment(temp, f);				// temp = f	
	//multiplicate_trans_with_vector(temp, f);	// f = At * f

	multiplicate_with_vector(x, temp);			// r0 = At * (f -  A * x0)
	vector_subtract(f, temp, temp);
	multiplicate_trans_with_vector(temp, r);

	vector_assigment(z, r);						// z0 = r0

	min_discrepancy = std::numeric_limits<double>::infinity();

	norm_f = calc_norm(f);
	scalar_rr_prev = scalar_product(r, r);		// (rk-1, rk-1)
	do {
		if (iter > max_iter) {
			//vector_assigment(x, x_min);
			break;
		}

		multiplicate_with_vector(z, temp);			// temp2 = At * A * z
		multiplicate_trans_with_vector(temp, temp2);

		a = scalar_rr_prev / scalar_product(temp2, z);	// a = ...

		vector_mult_on_const(z, a, temp);			// xk = xk-1 + ak * zk-1
		vector_add(x, temp, x);

		// step out:
		if (calc_norm(temp) / calc_norm(x) < step_eps) {
			//std::cout << "step OUT!\n"; 
			break;
		}

		vector_mult_on_const(temp2, a, temp);		// rk = rk-1 - At * A * zk-1
		vector_subtract(r, temp, r);


		scalar_rr = scalar_product(r, r);			// b = (rk,rk)/(rk-1, rk-1)
		b = scalar_rr / scalar_rr_prev;
		scalar_rr_prev = scalar_rr;

		vector_mult_on_const(z, b, temp);			// zk = rk + b * zk-1
		vector_add(r, temp, z);

		discrepancy = calc_norm(r) / norm_f;		// discrepancy
		iter++;
//#ifdef CALCULATE
//		printf("\riter: %d\t discr: %.6e", iter, discrepancy);
//#endif
		if (discrepancy < min_discrepancy) {
			min_discrepancy = discrepancy;
			vector_assigment(x_min, x);
		}

#ifdef _DEBUG
		multiplicate_with_vector(x, temp);
		vector_subtract(temp, f, temp);
		discrepancy = calc_norm(temp);

		
		std::cout << "rr: " << discrepancy << " Discrepancy " << calc_norm(temp) << std::endl;

#endif


	} while (discrepancy > eps);

	vector_assigment(x, x_min);
	return iter;
}

void not_symmetrical_sparse_matrix::
	MCG_LU(std::vector<double>& l, std::vector<double>& ud, std::vector<double>& u,
	std::vector<double>& x, std::vector<double>& x_min, std::vector<double>& f,
	std::vector<double>& r, std::vector<double>& z,
	std::vector<double>& temp,std::vector<double>& temp2)
{
	double discrepancy, min_discrepancy, norm_f, scalar_rr, scalar_rr_prev;
	double a, b;
	int iter = 0;

	// partial LU decomposition	
	LU_partial_decomposite(l, ud, u);

	// By = g
	// B = U-t * At * L-t * L-1 * A * U-1
	// y = U * x
	// g = U-t * At * L-t * L-1 * f

	// r = U-t * At * L-t * L-1 * ( f - Ax )
	multiplicate_with_vector(x, temp);
	vector_subtract(f, temp, temp);  
	forward_gauss(l, temp2, temp);
	backward_gauss(l, temp, temp2);
	multiplicate_trans_with_vector(temp, temp2);
	forward_gauss(u, ud, r, temp2);		// r = ...

	vector_assigment(z, r);				//  z = r

	vector_assigment(temp, x);
	multiplicate_ufactor(u, ud, temp, x);			// x = Ux

	min_discrepancy = std::numeric_limits<double>::infinity();
	norm_f = calc_norm(f);
	scalar_rr_prev = scalar_product(r, r);

	do {
		if (iter > max_iter) {
			vector_assigment(x, x_min);
#ifdef CALCULATE
			printf("\riter > max_iter = %d", max_iter);
#endif
			break;
		}

		// ak = ...
		vector_assigment(temp, z);
		backward_gauss(u, ud, temp2, temp);
		multiplicate_with_vector(temp2, temp);
		forward_gauss(l, temp2, temp);
		backward_gauss(l, temp, temp2);
		multiplicate_trans_with_vector(temp, temp2);
		forward_gauss(u, ud, temp, temp2);

		a = scalar_rr_prev / scalar_product(temp, z);

		// rk = ...
		vector_mult_on_const(temp, a, temp); 
		vector_subtract(r, temp, r);

		//x'k = x'k-1 - a * z
		vector_mult_on_const(z, a, temp);
		vector_add(x, temp, x);

		// step out:
		if (calc_norm(temp) / calc_norm(x) < step_eps) {
#ifdef CALCULATE
			printf("\nstep out\n");
#endif
			break;
		}

		// b = ...
		scalar_rr = scalar_product(r, r);
		b = scalar_rr / scalar_rr_prev;
		scalar_rr_prev = scalar_rr;

		// z = ...
		vector_mult_on_const(z, b, temp);
		vector_add(r, temp, z);

		discrepancy = calc_norm(r) / norm_f;		// discrepancy
		iter++;
#ifdef CALCULATE
		printf("\riter: %d\t discr: %.6e", iter, discrepancy);
#endif
		if (discrepancy < min_discrepancy) {
			min_discrepancy = discrepancy;
			vector_assigment(x_min, x);
		}
	} while (discrepancy > eps);

	///// x = U-1 * x'
	vector_assigment(temp, x);
	backward_gauss(u, ud, x, temp);
}



void not_symmetrical_sparse_matrix::
	MCG_DIAG(std::vector<double>& diag, std::vector<double>& x, std::vector<double>& x_min, std::vector<double>& f,
	std::vector<double>& r, std::vector<double>& z, 
	std::vector<double>& temp, std::vector<double>& temp2) {

	double discrepancy, min_discrepancy, norm_f, scalar_rr, scalar_rr_prev;
	double a, b;
	int iter = 0;

	// L = U = sqrt(di)
	diag_decomposite(diag);

	// r = U-t * At * L-t * L-1 * ( f - Ax )
	multiplicate_with_vector(x, temp);
	vector_subtract(f, temp, temp);
	diag_gauss(diag, temp2, temp);
	diag_gauss(diag, temp, temp2);
	multiplicate_trans_with_vector(temp, temp2);
	diag_gauss(diag, r, temp2);		// r = ...

	vector_assigment(z, r);				//  z = r

	vector_assigment(temp, x);
	//multiplicate_ufactor(u, ud, temp, x);			// x = Ux
	multiplicate_diag_with_vector(temp, x);

	min_discrepancy = std::numeric_limits<double>::infinity();
	norm_f = calc_norm(f);
	scalar_rr_prev = scalar_product(r, r);

	do {
		if (iter > max_iter) {
			vector_assigment(x, x_min);
#ifdef CALCULATE
			printf("\niter > max_iter = %d\n", max_iter);
#endif
			break;
		}

		// ak = ...
		vector_assigment(temp, z);
		diag_gauss(diag, temp2, temp);
		multiplicate_with_vector(temp2, temp);
		diag_gauss(diag, temp2, temp);
		diag_gauss(diag, temp, temp2);
		multiplicate_trans_with_vector(temp, temp2);
		diag_gauss(diag, temp, temp2);

		a = scalar_rr_prev / scalar_product(temp, z);

		// rk = ...
		vector_mult_on_const(temp, a, temp);
		vector_subtract(r, temp, r);

		//x'k = x'k-1 - a * z
		vector_mult_on_const(z, a, temp);
		vector_add(x, temp, x);

		// step out:
		if (calc_norm(temp) / calc_norm(x) < step_eps) {
#ifdef CALCULATE
			printf("\nstep out\n");
#endif
			break;
		}
		// b = ...
		scalar_rr = scalar_product(r, r);
		b = scalar_rr / scalar_rr_prev;
		scalar_rr_prev = scalar_rr;

		// z = ...
		vector_mult_on_const(z, b, temp);
		vector_add(r, temp, z);

		discrepancy = calc_norm(r) / norm_f;		// discrepancy
		iter++;
#ifdef CALCULATE
		printf("\riter: %d\t discr: %.6e", iter, discrepancy);
#endif
		if (discrepancy < min_discrepancy) {
			min_discrepancy = discrepancy;
			vector_assigment(x_min, x);
		}
	} while (discrepancy > eps);

	// x = U-1 * x'
	vector_assigment(temp, x);
	diag_gauss(diag, x, temp);
}


int not_symmetrical_sparse_matrix::BCGSTAB(std::vector<double> &l, std::vector<double> & ud, std::vector<double> & u,
	std::vector<double> & x, std::vector<double> & f, std::vector<double> & r, std::vector<double> &r_0, 
	std::vector<double> & z, std::vector<double> & temp, std::vector<double> & temp2, std::vector<double> &v, std::vector<double> &p,
	std::vector<double> & y, std::vector<double> &h, std::vector<double> &s, std::vector<double> &t) {

	//eps = 1e-15;

	double discrepancy;
	int iter = 0;

	// partial LU decomposition	
	LU_partial_decomposite(l, ud, u);

	// r0 = Ax - b
	multiplicate_with_vector(x, temp);
	vector_subtract(f, temp, r_0);

	vector_assigment(r, r_0);
	// z0 = r0
	z = r_0;

	double ro = 1, roi;
	double a = 1;
	double w = 1;
	double b = 0;

	std::vector<double> xPrev = x;

	// iteration process:
	do {
		iter++;

		// pi = 
		roi = scalar_product(r_0, r);

		// b = 
		b = (roi / ro) * (a / w);
		ro = roi;

		// p = r(i-1) + b(p(i-1) - w(i-1)*v(i-1))
		vector_mult_on_const(v, w, temp);
		vector_subtract(p, temp, p);
		vector_mult_on_const(p, b, p);
		vector_add(p, r, p);

		// y = U^-1  * L^-1 * p <--> LUy = p
		forward_gauss(l, temp, p);
		backward_gauss(u, ud, y, temp);

		// v = Ay
		multiplicate_with_vector(y, v);

		// 
		a = roi / scalar_product(r_0, v);

		// 
		vector_mult_on_const(y, a, temp);
		vector_add(x, temp, h);

		/* // check if "h' is accurate enough then xi = h and QUIT
		multiplicate_with_vector(h, temp);
		vector_subtract(temp, f, temp);
		discrepancy = calc_norm(temp);
		

		if (discrepancy < eps) {
			vector_assigment(x, h);
			break;
		}
		*/

		// 
		vector_mult_on_const(v, a, temp);
		vector_subtract(r, temp, s);

		//
		forward_gauss(l, temp, s);
		backward_gauss(u, ud, z, temp);

		// 
		multiplicate_with_vector(z, t);

		// 
		forward_gauss(l, temp, t);
		forward_gauss(l, temp2, s);

		w = scalar_product(temp, temp2) / scalar_product(temp, temp);

		// 
		vector_mult_on_const(z, w, temp);
		vector_add(h, temp, x);

		// check step exit:
		vector_subtract(xPrev, x, temp);
		double stepValue = calc_norm(temp) / calc_norm(x);

		if (stepValue < step_eps) {
			break;
		}


		vector_assigment(xPrev, x); 

		//
		vector_mult_on_const(t, w, temp);
		vector_subtract(s, temp, r);


		//std::cout << "dscr: " << discrepancy << "r:" << calc_norm(r) << std::endl;

		discrepancy = calc_norm(r);
		if (discrepancy < eps) {
			break;
		}

	} while (iter < max_iter); // exit by max iter

	return iter;
}

#pragma endregion

void not_symmetrical_sparse_matrix::setFirstBoundaryConditionsLeft() {
	di[0] = 1;	// set first sin elem
	di[1] = 1;	// set first cos elem

	// zeroing another off-diagonal elems:
	ggl[0] = 0;
	for(int i = 0; i <= 4; i++) ggu[i] = 0;
}

void not_symmetrical_sparse_matrix::setFirstBoundaryConditionsRight() {
	int last = di.size() - 1;
	di[last - 1] = 1;	// set last sin elem
	di[last] = 1;		// set last cos elem


	// zeroing another off-diagonal elems:
	last = ggl.size() - 1;
	ggu[last] = 0;
	for(int i = 0 ; i <= 4; i++) ggl[last--] = 0;
}