#include "DFT.h"

My_DFT::My_DFT(std::string fileName)
{
	x_input = ReadFile(fileName);

	N = x_input.size();

	x_IDFT.resize(N, 0.0);
	y_DFT.resize(N, 0.0);
	x_IFFT.resize(N, 0.0);
	y_FFT.resize(N, 0.0);
}

My_DFT::My_DFT(ComplexVect input)
{
	x_input = input;

	N = x_input.size();

	x_IDFT.resize(N, 0.0);
	y_DFT.resize(N, 0.0);
	x_IFFT.resize(N, 0.0);
	y_FFT.resize(N, 0.0);
}

inline complex My_DFT::omega(int k, int l, int sign)
{
	return exp(complex(0, (double)sign * 2.0 * M_PI * (double)l / (double)(pow(2., k))));
}

inline complex My_DFT::exponent(int N, int j, int k, int sign)
{
	return exp(complex(0, (double)sign * 2.0 * M_PI * (double)j * (double)k / ((double)N)));
}

ComplexVect My_DFT::DiscreteFourierTransform(int sign) // DPF formula 1
{
	for (int k = 0; k < N; k++)
	{
		for (int j = 0; j < N; j++)
		{
			y_DFT[k] = y_DFT[k] + x_input[j] * exponent(N, j, k, sign);
		}

		y_DFT[k] = y_DFT[k] / sqrt(N);
	}

	return y_DFT;
}

ComplexVect My_DFT::InverseDiscreteFourierTransform(int sign)
{
	for (int j = 0; j < N; j++)
	{
		for (int k = 0; k < N; k++)
		{
			x_IDFT[j] = x_IDFT[j] + y_DFT[k] * exponent(N, j, k, sign);
		}
		x_IDFT[j] = x_IDFT[j] / sqrt(N);
	}
	
	return x_IDFT;
}

ComplexVect My_DFT::FastFourierTransform(int sign)
{
	int n = log2(N);

	for (int k = 1; k < n + 1; k++)
	{
		for (int j = 0; j < (1 << (k - 1)); j++)
		{
			for (int l = 0; l < (1 << (n - k)); l++)
			{
				int even = j * (1 << (n + 1 - k)) + l;
				int odd = even + (1 << (n - k));
				y_FFT[even] = x_input[even] + x_input[odd];
				y_FFT[odd] = (x_input[even] - x_input[odd]) * omega(n + 1 - k, l, sign);
			}
		}
	}

	ComplexVect z(y_FFT);

	for (int k = 1; k < N - 1; k++)
	{
		int temp = 1;
		int inverse = 0;
		for (int j = 0; j < n; j++)
		{
			if (temp & k)
				inverse += (1 << (n - 1 - j));
			temp = temp << 1;
		}
		z[k] = y_FFT[inverse];
	}

	y_FFT = z;

	for (size_t i = 0; i < N; i++)
		y_FFT[i] /= (double)sqrt(N);
	/*for (size_t i = 0; i < N; i++)
		y_FFT[i] *= 2;*/

	return y_FFT;
}

ComplexVect My_DFT::InverseFastFourierTransform(int sign)
{
	int n = log2(N);

	for (int k = 1; k < n + 1; k++)
	{
		for (int j = 0; j < (1 << (k - 1)); j++)
		{
			for (int l = 0; l < (1 << (n - k)); l++)
			{
				int even = j * (1 << (n + 1 - k)) + l;
				int odd = even + (1 << (n - k));
				x_IFFT[even] = y_FFT[even] + y_FFT[odd];
				x_IFFT[odd] = (y_FFT[even] - y_FFT[odd]) * omega(n + 1 - k, l, sign);
			}
		}
	}

	ComplexVect z(x_IFFT);

	for (int k = 1; k < N - 1; k++)
	{
		int temp = 1;
		int inverse = 0;
		for (int j = 0; j < n; j++)
		{
			if (temp & k)
				inverse += (1 << (n - 1 - j));
			temp = temp << 1;
		}
		z[k] = x_IFFT[inverse];
	}

	x_IFFT = z;

	for (size_t i = 0; i < N; i++)
		x_IFFT[i] /= (double)sqrt(N);// x[i] = x[i] / sqrt
	/*for (size_t i = 0; i < N; i++)
		x_IFFT[i] *= 2;// x[i] = x[i] / sqrt*/

	return x_IFFT;
}

ComplexVect My_DFT::Convolution(ComplexVect x, ComplexVect y)
{
	int size_x1 = x.size();
	int size_x2 = y.size();
	int max_size = std::max(size_x1, size_x2);

	ComplexVect var1, var2;
	ComplexVect result(2 * max_size);

	var1 = x;
	var2 = y;
	var1.resize(2 * max_size);
	var2.resize(2 * max_size);

	for (int i = 0; i < 2 * max_size; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			if (i - j >= 0)
				result[i] = result[i] + var1[j] * var2[i - j];
		}
	}

	return result;
}

ComplexVect My_DFT::ConvolutionFFT(ComplexVect x, ComplexVect y)
{
	int size_x1 = x.size();
	int size_x2 = y.size();
	int max_size = std::max(size_x1, size_x2);

	ComplexVect result(2 * max_size);

	ComplexVect var1 = x;
	ComplexVect var2 = y;
	var1.resize(2 * max_size);
	var2.resize(2 * max_size);

	int count = 1;
	max_size = 2 * max_size;

	while (true)
	{
		if (max_size < pow(2., count))
		{
			max_size = pow(2., count);
			var1.resize(max_size);
			var2.resize(max_size);
			result.resize(max_size);
			break;
		}

		count++;
	}

	var1 = My_DFT(var1).FastFourierTransform(-1);
	var2 = My_DFT(var2).FastFourierTransform(-1);

	for (unsigned int k = 0; k < max_size; k++)
	{
		result[k] = var1[k] * var2[k] * sqrt(max_size);
	}

	result = My_DFT(result).FastFourierTransform(1);

	return result;
}

ComplexVect My_DFT::ConvolutionFFT(ComplexVect x, ComplexVect y, int result_size, int max_size)
{
	ComplexVect vector_y;
	for (int i = 0; i < 2 * max_size; i++)
	{
		y.push_back(sqrt(2 * max_size) * x[i] * y[i]);
	}

	auto vector_result = My_DFT(y).FastFourierTransform(-1);

	for (int i = 0; i < 2 * max_size - result_size; i++)
	{
		vector_result.pop_back();
	}

	return vector_result;
}

ComplexVect ReadFile(std::string filename)
{
	std::ifstream in(filename);
	ComplexVect tmp;
	if (in.is_open())
	{
		double var2, var3;

		if (!in.eof())
			in >> var2 >> var3;

		while (!in.eof())
		{
			tmp.push_back(complex(var2, var3));
			in >> var2 >> var3;
		}

		in.close();
	}
	else
	{
		throw "File couldn't be opened!";
	}
	return tmp;
}

bool WriteFile(ComplexVect vector_r, std::string filename)
{
	std::ofstream in(filename);

	if (!in) 
	{
		throw "File couldn't be opened!";
		return false;
	}

	for (unsigned int i = 0; i < vector_r.size(); i++) 
	{
		in << "    " << vector_r[i].real() << "    " << vector_r[i].imag() << std::endl;
	}

	in.close();
	return true;
}

bool Comparison(std::string file1, std::string file2)
{
	std::ifstream in(file1);
	ComplexVect first = ReadFile(file1);
	ComplexVect second = ReadFile(file2);
	size_t size = first.size();
	for (size_t i = 0; i < size; i++)
	{
		if (!(std::fabs(first[i].real() - second[i].real()) < std::numeric_limits<double>::epsilon()))
			return false;
		if (!(std::fabs(first[i].imag() - second[i].imag()) < std::numeric_limits<double>::epsilon()))
			return false;
	}
	return true;
}

void InaccuracyComparison(std::string file1, std::string file2, std::string Inaccuracy)
{
	std::ifstream in(file1);
	ComplexVect first = ReadFile(file1);
	ComplexVect second = ReadFile(file2);
	size_t size = first.size();
	ComplexVect result(second);
	for (size_t i = 0; i < size; i++)
	{
		result[i] -= first[i];
	}

	auto max = max_element(result.begin(), result.end(), [](complex v1, complex v2) {return abs(v1) < abs(v2); });

	for (size_t i = 0; i < size; i++)
	{
		std::cout << abs(result[i]) << std::endl;
	}

	std::cout << "Max inaccuracy: ";
	//std::cout.setf(std::ios::fixed);
	std::cout << abs(*max) << std::endl;

	std::ofstream out(Inaccuracy);

	if (!out)
	{
		throw "File couldn't be opened!";
	}

	for (unsigned int i = 0; i < result.size(); i++)
	{
		out << "    " << result[i].real() << "    " << result[i].imag() << std::endl;
	}

	out.close();
}
