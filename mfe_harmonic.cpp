#include "pch.h"
#include "mfe_harmonic.h"
#include <fstream>
#include "assert.h"
#include <time.h>
#include <iostream>
#include <chrono>

void HarmonicTask::setParams() {
	//lambda = std::vector<double>(amountSubareas);
	//sigma = std::vector<double>(amountSubareas);

	w = 2;

	lambda[0] = 3;
	sigma[0] = 10;
	epsilon[0] = 6;

	uSinExact = [](const double& x) {return x; };
	uCosExact = [](const double& x) {return 3*x; };

	//  fSin = -div(lmd * grad(Usin) - w * sgm * Ucos - w*w*eps*Usin
	fFuncSin = [](const double x, const double lmd, const double sgm, const double eps, const double w)
		{return -lmd * 0 - w * sgm * 3*x - w * w * eps * 1*x; };
	//  fCos = -div(lmd * grad(Ucos) + w * sgm * Usin - w*w*eps*Ucos
	fFuncCos = [](const double x, const double lmd, const double sgm, const double eps, const double w) 
		{return -lmd * 0 + w * sgm * 1*x - w * w * eps * 3*x; };

}

void HarmonicTask::init() {
#pragma region outputResultFile
	fout.open(R"(output/result.txt)");
	exploreout.open(R"(output/explore.txt)");
	exploreout.precision(15);
	//fout.precision(17);
#pragma endregion

#pragma region InitSpaceGrid
{
	nodes = std::vector<double>();
	elems = std::vector<FiniteElem>();
	subareas = std::vector<int>();

	std::fstream fin(R"(input/grid_space.txt)");

	int numSpaceGridDividion = 1;
	int numOfAreas = 1;

	fin >> numSpaceGridDividion >> numOfAreas;

	int numOfElems;
	double xStart, step, coef;
	fin >> xStart >> numOfElems >> step >> coef;

	int k = pow(2, numSpaceGridDividion - 1);

	// calculate grid parameters for unevenness:
	numOfElems *= k;
	coef = pow(coef, 1.0 / k);
	// calculate first step
	double stepsCoef = 0;
	for (int i = 0; i < k; i++) stepsCoef += pow(coef, i);
	step /= stepsCoef;

	double x = xStart;	
	nodes.push_back(x);		// add x0 in nodes
	for (int i = 0; i < numOfElems; i++) {
		x += step;
		nodes.push_back(x);

		step *= coef;					// change step
	}

	// fill elems array:
	int j = 0;
	for (; j < nodes.size() - 1; j++) {
		elems.push_back(FiniteElem{ j, j + 1});
	}

	assert(j == nodes.size() - 1);

	// fill subareas:
	int  numOfFiniteElems = 0;
	int sum = 0;

	fin >> amountSubareas;
	for (int i = 0; i < amountSubareas; i++) {
		fin >> numOfFiniteElems;
		numOfFiniteElems *= k; // consideration of grid dividion
		sum += numOfFiniteElems;
		for(int j = 0; j < numOfFiniteElems; j++) subareas.push_back(0);
	}	
	
	assert(sum == elems.size());

	// init vector of params of equals in subareas:
	lambda = std::vector<double>(amountSubareas);
	sigma = std::vector<double>(amountSubareas);
	epsilon = std::vector<double>(amountSubareas);

	fin.close();
}
#pragma endregion

#pragma region MatrixInit
{
	// memory allocation for matrix
	std::fstream fin(R"(input/solveParams.txt)");
	double eps = 1e-10;
	int maxIter = 3000;
	fin >> eps >> maxIter;

	globalMatrix.init(elems.size(), maxIter, eps);
	fin.close();

	// memory allocation for solver: (LOS with LU)
	const int matrixDim = globalMatrix.get_dim();
	const int amountElems = globalMatrix.get_amount_elems();
	x_min = std::vector<double>(matrixDim);
	r = std::vector<double>(matrixDim);
	r_0 = std::vector<double>(matrixDim);

	z = std::vector<double>(matrixDim);
	p = std::vector<double>(matrixDim);
	temp = std::vector<double>(matrixDim);
	temp2 = std::vector<double>(matrixDim);

	v = std::vector<double>(matrixDim);
	y = std::vector<double>(matrixDim);
	h = std::vector<double>(matrixDim);
	s = std::vector<double>(matrixDim);
	t = std::vector<double>(matrixDim);

	ud = std::vector<double>(matrixDim);
	l = std::vector<double>(amountElems);
	u = std::vector<double>(amountElems);
}
#pragma endregion

#pragma region MemoryAllocation
	const int matrixDim = globalMatrix.get_dim();
	// for result:
	q = std::vector<double>(matrixDim);

	// local matrix of mass and rigidity, and local vector;
	const int pointsToFiniteElem = 4;
	localMatrix = LocalMatrix(pointsToFiniteElem);
	for (auto& el : localMatrix) el = std::vector<double>(pointsToFiniteElem);
	
	// rigth part:
	fLocal = std::vector<double>(pointsToFiniteElem);
	f = std::vector<double>(matrixDim);

#pragma endregion
}

void HarmonicTask::resetTask() {
	globalMatrix.resetMatrix();	// reset matrix
	for (auto& el : f) el = 0; // reset right part
	calculateGlobalMatrixAndRightPart();

	// zeroing help vectors:
	std::fill(q.begin(), q.end(), 0);
	std::fill(x_min.begin(), x_min.end(), 0);

	std::fill(r.begin(), r.end(), 0);
	std::fill(r_0.begin(), r_0.end(), 0);
	std::fill(z.begin(), z.end(), 0);
	std::fill(p.begin(), p.end(), 0);
	std::fill(temp.begin(), temp.end(), 0);
	std::fill(temp2.begin(), temp2.end(), 0);
	std::fill(v.begin(), v.end(), 0);
	std::fill(y.begin(), y.end(), 0);
	std::fill(h.begin(), h.end(), 0);
	std::fill(s.begin(), s.end(), 0);
	std::fill(t.begin(), t.end(), 0);

	std::fill(ud.begin(), ud.end(), 0);
	std::fill(u.begin(), u.end(), 0);
	std::fill(l.begin(), l.end(), 0);
}

void HarmonicTask::calculateGlobalMatrixAndRightPart() {
	// calculate global matrix and vector:
	for (int elemNum = 0; elemNum < elems.size(); elemNum++) {
		// calculate local parts:
		calculateLocalMatrix(elemNum);	
		calculateLocalRightPart(elemNum);
	
		// add local parts to global
		globalMatrix.addLocalMatrix(elemNum, localMatrix);
		addLocalRigtPartToGlobal(elemNum);
	}

	// set boundary conditions
	setFirstBoundaryConditions();
}


void HarmonicTask::solve() {

	//int not_symmetrical_sparse_matrix::MCG(std::vector<double>& x, std::vector<double>& x_min,
	//	std::vector<double>& f, std::vector<double>& r, std::vector<double>& z,
	//	std::vector<double>& temp, std::vector<double> temp2) {

	w = 1;

	lambda[0] = 1;
	sigma[0] = 1;
	epsilon[0] = 1;



	enum methods {MCG, LU, BCG, MCG_LU};

	methods method = MCG;


	int iter = 0;

	resetTask();

	switch(method) {
		case methods::MCG: iter = globalMatrix.MCG(q, x_min, f, r, z, temp, temp2); break;
		case methods::LU:	globalMatrix.directSolveLU(l, ud, u, temp, q, f);  break;
		case methods::BCG: iter = globalMatrix.BCGSTAB(l, ud, u, q, f, r, r_0, z, temp, temp2, v, p, y, h, s, t);  break;
	}

	saveResult(iter);

	//void not_symmetrical_sparse_matrix::LOS(std::vector<double>& x,
	//	std::vector<double>& f, std::vector<double>& r, std::vector<double>& z,
	//	std::vector<double>& p, std::vector<double>& temp) {

	//resetTask();
	//globalMatrix.LOS(q, f, r, z, p, temp);

	//saveResult(iter);

}

void HarmonicTask::explore() {
	// method for explore function on differet solvers.

	// set difference param to explore:
	w = 1;

	lambda[0] = 1;
	sigma[0] = 1;
	epsilon[0] = 1;
	//


	std::vector<double> vW = { 1e-4, 1e-3, 1e-2, 1e-1, 10, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9 };
	std::vector<double> vLambda = { 1e2, 1e3, 1e4, 1e5, 8 * 1e5 };
	std::vector<double> vSigma = { 0, 10, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8 };
	std::vector<double> vEpsilon = { 8.81*1e-12,  1e-12, 1e-11, 1e-10 };


	// table title:
	exploreout << "Parameter\tMSGiter\tMSGtime\tMSGsinErr\tMSGcosErr\t";
	//exploreout << "LOSiter\tLOStime\tLOSsinErr\tLOScosErr\t";
	exploreout << "LUtime\tLUsinErr\tLUcosErr\t";
	exploreout << "BSGSTABiter\tBSGSTABtime\tBSGSTABsinErr\tBSGSTABcosErr\n";

	double cosErr, sinErr;
	int i = 0, iter = 0;

	const auto & vExploredParameter = vEpsilon;
	auto &param_place = epsilon[0];
	const auto sizeExplore = vExploredParameter.size();

	using timer = std::chrono::high_resolution_clock;
	timer::time_point CalcTimeStart, CalcTimeStop;
	double time;

	// exploration:
	for(auto value: vExploredParameter) {
		exploreout << value << "\t";
		param_place = value;
		//exploreout << w << "\t";

		// explore MSG
		// reset matrix and right part
		std::cout << "\riteration " << i << " of " << sizeExplore << " MSG is calculating";

		resetTask();
		CalcTimeStart = timer::now();
		iter = globalMatrix.MCG(q, x_min, f, r, z,temp, temp2);
		CalcTimeStop = timer::now();

		time = std::chrono::duration_cast<std::chrono::microseconds>(CalcTimeStop - CalcTimeStart).count() / 1000000.0;
		calculateError(sinErr, cosErr);
		saveExploreResult(iter, time, sinErr, cosErr);
			   

		// explore direcxt LU
		// reset matrix and right part
		std::cout << "\riteration " << i << " of " << sizeExplore << " LU is calculating";

		resetTask();
		CalcTimeStart = timer::now();
		globalMatrix.directSolveLU(l, ud, u, temp, q, f);
		CalcTimeStop = timer::now();

		time = std::chrono::duration_cast<std::chrono::microseconds>(CalcTimeStop - CalcTimeStart).count() / 1000000.0;
		calculateError(sinErr, cosErr);
		saveExploreResult(time, sinErr, cosErr);

		// explore BCGSTAB with LU
		// reset matrix and right part
		std::cout << "\riteration " << i << " of " << sizeExplore << "\t\tBCGSTAB is calculating";

		resetTask();
		CalcTimeStart = timer::now();
		iter = globalMatrix.BCGSTAB(l, ud, u, q, f, r, r_0, z, temp, temp2, v, p, y, h, s, t);
		CalcTimeStop = timer::now();

		time = std::chrono::duration_cast<std::chrono::microseconds>(CalcTimeStop - CalcTimeStart).count() / 1000000.0;
		calculateError(sinErr, cosErr);
		saveExploreResult(iter, time, sinErr, cosErr);

		exploreout << std::endl;	
		i++;
	}
}

void HarmonicTask::calculateLocalMatrix(uint32_t elemNum) {
	// approximate div(lambda grad) part:								
	const auto& elem = elems[elemNum];	
	const double step = nodes[elem.right] - nodes[elem.left];

	const double coefLambda = lambda[subareas[elemNum]] / (step); 
	const double coefSigma = (step / 6) * sigma[subareas[elemNum]];
	const double coefEpsilon = (step / 6) * epsilon[subareas[elemNum]];

	// block: i = 1, j = 1
	localMatrix[0][0] = coefLambda - w*w * 2 * coefEpsilon;
	localMatrix[1][0] = w * 2 * coefSigma;
	localMatrix[0][1] = -localMatrix[1][0];
	localMatrix[1][1] = localMatrix[0][0];

	// block: i = 1, j = 2
	localMatrix[0][2] = -1 * coefLambda - w * w * 1 * coefEpsilon;
	localMatrix[1][2] = w * 1 * coefSigma;
	localMatrix[0][3] = -localMatrix[1][2];
	localMatrix[1][3] = localMatrix[0][2];

	// block: i = 2, j = 1
	localMatrix[2][0] = -1 * coefLambda - w * w * 1 * coefEpsilon;
	localMatrix[3][0] = w * 1 * coefSigma;
	localMatrix[2][1] = -localMatrix[3][0];
	localMatrix[3][1] = localMatrix[2][0];

	// block: i = 2, j = 2
	localMatrix[2][2] = coefLambda - w * w * 2 * coefEpsilon;
	localMatrix[3][2] = w * 2 * coefSigma;
	localMatrix[2][3] = -localMatrix[3][2];
	localMatrix[3][3] = localMatrix[2][2];
}

void HarmonicTask::calculateLocalRightPart(uint32_t elemNum) {
	const auto& elem = elems[elemNum];
	const auto& leftNode = nodes[elem.left];
	const auto& rightNode = nodes[elem.right];

	const double lmd = lambda[subareas[elemNum]];
	const double sgm = sigma[subareas[elemNum]];
	const double eps = epsilon[subareas[elemNum]];

	const double f0Sin = fFuncSin( leftNode, lmd, sgm, eps, w);
	const double f1Sin = fFuncSin(rightNode, lmd, sgm, eps, w);
	const double f0Cos = fFuncCos (leftNode, lmd, sgm, eps, w);
	const double f1Cos = fFuncCos(rightNode, lmd, sgm, eps, w);

	const auto k =  (rightNode - leftNode) / 6;

	fLocal[0] = k * (2 * f0Sin + f1Sin);
	fLocal[1] = k * (2 * f0Cos + f1Cos);
	fLocal[2] = k * (f0Sin + 2 * f1Sin);
	fLocal[3] = k * (f0Cos + 2 * f1Cos);
}

void HarmonicTask::addLocalRigtPartToGlobal(uint32_t elemNum) {
	const int elemsInLocal = 4;
	int place = elemNum * 2;

	for (int i = 0; i < elemsInLocal; i++, place++) f[place] += fLocal[i];
}

void HarmonicTask::setFirstBoundaryConditions() {
	// set first bounday conditions in left side:
	globalMatrix.setFirstBoundaryConditionsLeft();
	f[0] = uSinExact(nodes.front());
	f[1] = uCosExact(nodes.front());

	// set first bounday conditions in right side:
	globalMatrix.setFirstBoundaryConditionsRight();
	f[f.size() - 2] = uSinExact(nodes.back()); /*[nodes.size() - 1]*/
	f[f.size() - 1] = uCosExact(nodes.back());
}

void vectorSubtraction(std::vector<double>& result, const std::vector<double>& a){
	for (int i = 0; i < result.size(); i++) result[i] -= a[i];
}

double calcNorm(const std::vector<double> &x) {
	double norm = 0;
	for (int i = 0; i < x.size(); i++) {
		norm += x[i] * x[i];
	}
	norm = sqrt(norm);
	return norm;
}

/*
//#include <iostream>
bool HarmonicTask::SimpleIterationDiscrepOut() {
	// || A(qi) * qi - b(qi) || / || b(qi) || < eps => out:
	
	temp = globalMatrix.multiplicate_with_vector(q, temp);
	vectorSubtraction(temp, f);

	double resultNorm = calcNorm(temp) / calcNorm(f);
	//std::cout << resultNorm << std::endl;

	return resultNorm < epsDiscrep;
	
	return true;
}
*/


// +
void HarmonicTask::saveResult(int iter) {
	double sumCos = 0;
	double sumSin = 0;

	double uSinAccurate = 0;
	double uCosAccurate = 0;
	double uSinError = 0;
	double uCosError = 0;
	double coordinate = 0;

	fout.scientific;
	fout.precision(3);

	fout << "x\texactSin\tnumSin\t errorSin\texactCos\tnumCos\terrorCos\n";
	for (int i = 0; i < nodes.size(); i++) {
		// output node: x - coordinate
		coordinate = nodes[i];
		fout << coordinate << '\t';

		// output accurate sin result
		uSinAccurate = uSinExact(coordinate);
		fout << uSinAccurate << '\t';

		// output numerical sin result
		fout << q[2 * i] << '\t';

		// outpu sin absolute error
		uSinError = uSinAccurate - q[2 * i];
		fout << uSinError << '\t';

		// output accurate sin result
		uCosAccurate = uCosExact(coordinate);
		fout << uCosAccurate << '\t';

		// output numerical sin result
		fout << q[2 * i + 1] << '\t';

		// outpu Cos absolute error
		uCosError = uCosAccurate - q[2 * i + 1];
		fout << uCosError << "\t\n";

		// accumulate error:
		sumSin += uSinError * uSinError;
		sumCos += uCosError * uCosError;

	}
	sumSin = sqrt(sumSin);
	sumCos = sqrt(sumCos);

	fout << std::endl << "iters:" << iter;
	fout << std::endl << "sin norm error = " << sumSin;
	fout << std::endl << "cos norm error = " << sumCos << std::endl;
}

void HarmonicTask::calculateError(double & sinErr, double & cosErr) {
	double sumCos = 0;
	double sumSin = 0;

	double uSinAccurate = 0;
	double uCosAccurate = 0;
	double uSinError = 0;
	double uCosError = 0;
	double coordinate = 0;

	for (int i = 0; i < nodes.size(); i++) {
		// output node: x - coordinate
		coordinate = nodes[i];

		// sin absolute error
		uSinError = uSinExact(coordinate) - q[2 * i];

		// cos absolute error
		uCosError = uCosExact(coordinate) - q[2 * i + 1];

		// accumulate error:
		sumSin += uSinError * uSinError;
		sumCos += uCosError * uCosError;
	}

	// return errors
	sinErr = sqrt(sumSin);
	cosErr = sqrt(sumCos);
}

void HarmonicTask::saveExploreResult(const double time, const double sinErr, const double cosErr) {
	exploreout << time << '\t' << sinErr << '\t' << cosErr << '\t';
}

void HarmonicTask::saveExploreResult(const int iter, const double time, const double sinErr, const double cosErr) {
	exploreout  << iter << '\t' << time << '\t' << sinErr << '\t' << cosErr << '\t';
}

