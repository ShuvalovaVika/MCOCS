#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <complex>
#include <string>
#include <algorithm>

typedef std::complex<double> complex_temp;
typedef std::vector<complex_temp> ComplexVect;

class My_DFT {
public:
	ComplexVect x, y;
	ComplexVect coef;

	My_DFT(std::string fileName);
	
	complex_temp inline omega(int k, int l, int sign);
	complex_temp inline exponent(int N, int j, int k, int sign);
	ComplexVect DiscreteFourierTransform(int sign);
	ComplexVect InverseDiscreteFourierTransform(int sign);
	ComplexVect FastFourierTransform(int sign);
	ComplexVect InverseFastFourierTransform(int sign);
	ComplexVect Convolution(ComplexVect vector_x1, ComplexVect vector_x2);
	ComplexVect ConvolutionFFT(ComplexVect vector_x1, ComplexVect vector_x2);
	ComplexVect ConvolutionFFT(ComplexVect vector_x1, ComplexVect vector_x2, int result_size, int max_size);

	void ReadFile(std::string filename);
	bool WriteFile(ComplexVect vector_r, std::string filename);
};

bool Comparison(std::string file1, std::string file2);
void InaccuracyComparison(std::string file1, std::string file2, std::string Inaccuracy);