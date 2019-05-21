#include "pch.h"
#include <iostream>
#include "mfe_harmonic.h"

// solve 1D harmonic task with linear basis function
int main()
{
	HarmonicTask T;
	T.init();
	T.setParams();
	//T.solve();
	T.explore();
	std::cout << "\nend work";
}


