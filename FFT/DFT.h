#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <complex>
#include <string>
#include <algorithm>
#include <chrono>
#include <ctime>

typedef std::complex<double> complex;
typedef std::vector<complex> ComplexVect;

class My_DFT {
public:
	ComplexVect x_input;
	ComplexVect x_IDFT, y_DFT;
	ComplexVect x_IFFT, y_FFT;
	size_t N;

	My_DFT(std::string fileName);
	My_DFT(ComplexVect input);

	complex inline omega(int k, int l, int sign);
	complex inline exponent(int N, int j, int k, int sign);
	ComplexVect DiscreteFourierTransform(int sign);
	ComplexVect InverseDiscreteFourierTransform(int sign);
	ComplexVect FastFourierTransform(int sign);
	ComplexVect InverseFastFourierTransform(int sign);
	ComplexVect Convolution(ComplexVect vector_x1, ComplexVect vector_x2);
	ComplexVect ConvolutionFFT(ComplexVect vector_x1, ComplexVect vector_x2);
	ComplexVect ConvolutionFFT(ComplexVect vector_x1, ComplexVect vector_x2, int result_size, int max_size);
};

bool WriteFile(ComplexVect vector_r, std::string filename);
bool WriteFile(std::string ss, std::string filename);
ComplexVect ReadFile(std::string filename);
bool Comparison(std::string file1, std::string file2, std::string fileName);
void InaccuracyComparison(std::string file1, std::string file2, std::string Inaccuracy);