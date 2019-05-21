#include "pch.h"
#include <vector>

void save_vector(std::vector<double> &vect, std::string filename, std::string frm) {
	FILE* out;
	fopen_s(&out, filename.c_str(), "w");
	for (int i = 0; i < vect.size(); i++) fprintf(out, frm.c_str(), vect[i]);
	fclose(out);
}

void load_vector(int size, std::string filename, std::vector<double> &vect) {
	FILE *in;
	fopen_s(&in, filename.c_str(), "r");
	for (int i = 0; i < size; i++) fscanf_s(in, "%le ", &vect[i]);
	fclose(in);
}

void vector_subtract(std::vector<double>& a, std::vector<double>& b, std::vector<double>& result) {
	for (int i = 0; i < result.size(); i++) result[i] = a[i] - b[i];
}

void vector_add(std::vector<double>& a, std::vector<double>& b, std::vector<double>& result) {
	for (int i = 0; i < result.size(); i++) result[i] = a[i] + b[i];
}

double scalar_product(std::vector<double>& a, std::vector<double>& b) {
	double sum = 0;
	for (int i = 0; i < a.size(); i++)	sum += a[i] * b[i];
	return sum;
}

void vector_mult_on_const(std::vector<double>& a, double constant, std::vector<double>& result) {
	for (int i = 0; i < result.size(); i++) result[i] = a[i] * constant;
}

void vector_assigment(std::vector<double>& a, std::vector<double>& b) {
	for (int i = 0; i < a.size(); i++) a[i] = b[i];
}

double calc_norm(std::vector<double>& vector) {
	double sum = 0;
	for (int i = 0; i < vector.size(); i++) sum += pow(vector[i], 2);
	return sqrt(sum);
}