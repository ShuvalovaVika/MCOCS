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


/*inline complex My_DFT::omega(int k, int l, int sign)
{
	return exp(complex(0, (double)sign * 2.0 * M_PI * (double)l / (double)(powf(2., k))));
}*/

inline complex My_DFT::omega(int l, int k, int sign)
{
	complex i = complex(0, 1);
	complex omega = exp((-2 * (double)sign * M_PI * l / pow(2, k)) * i);
	return omega;
};


inline complex My_DFT::exponent(int N, int j, int k, int sign)
{
	return exp(complex(0, (double)sign * 2 * M_PI * (double)j * (double)k / ((double)N)));
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

ComplexVect My_DFT::InverseDiscreteFourierTransform(int sign) //ОДПФ 2 формула
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


int invert(int i, int n)
{
	int u = 0;
	int q;
	for (int k = 0; k < n; ++k)
	{
		q = i % 2;
		if (q == 1) u = u + pow(2, n - k - 1);
		i = i / 2;
	}
	return u;
}


ComplexVect My_DFT::FastFourierTransform(int sign)
{
	int n = log2(N);
	
	ComplexVect tmp(x_input);

	for (int k = 1; k < n + 1; k++)
	{
		for (int j = 0; j < (1 << (k - 1)); j++)
		{
			for (int l = 0; l < (1 << (n - k)); l++)
			{
				int even = j * (1 << (n + 1 - k)) + l;
				int odd = even + (1 << (n - k));
				y_FFT[even] = tmp[even] + tmp[odd];
				y_FFT[odd] = (tmp[even] - tmp[odd]) * omega(l, n + 1 - k, sign);
			}
		}

		tmp = y_FFT;
	}

	/*ComplexVect z(tmp);
	z[0] = tmp[0];

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
		z[k] = tmp[inverse];
	}
	tmp[N - 1] =z[N - 1];

	//z = tmp;
	//tmp = z;

	for (size_t i = 0; i < N; i++)
		tmp[i] /= (double)sqrt(N);
	return tmp;
}*/

	int u;
	complex d;
	for (int i = 0; i < N; ++i) //перестановка
	{
		u = invert(i, n);
		if (u >= i)
		{
			d = y_FFT[i];
			y_FFT[i] = y_FFT[u];
			y_FFT[u] = d;
		}
	}
	
	for (int i = 0; i < N; ++i) 
	{
		y_FFT[i] = y_FFT[i] / sqrt(N);
	}
	return y_FFT;
}


ComplexVect My_DFT::InverseFastFourierTransform(int sign)
{
	int n = log2(N);
	ComplexVect tmp(y_FFT);

	for (int k = 1; k < n + 1; k++)
	{
		for (int j = 0; j < (1 << (k - 1)); j++)
		{
			for (int l = 0; l < (1 << (n - k)); l++)
			{
				int even = j * (1 << (n + 1 - k)) + l;
				int odd = even + (1 << (n - k));
				x_IFFT[even] = (tmp[even] + tmp[odd]);
				x_IFFT[odd] = (tmp[even] - tmp[odd]) * omega(l, n - k + 1, sign);
			}
		}
		tmp = x_IFFT;
	}

	int u;
	complex d;
	for (int i = 0; i < N; ++i) //перестановка
	{
		u = invert(i, n);
		if (u >= i)
		{
			d = tmp[i];
			tmp[i] = tmp[u];
			tmp[u] = d;
		}
	}

	for (int i = 0; i < N; ++i) 
	{
		tmp[i] = tmp[i] / sqrt(N);
	}
	return tmp;
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

	var1 = My_DFT(var1).FastFourierTransform(1);
	var2 = My_DFT(var2).FastFourierTransform(1);

	for (unsigned int k = 0; k < max_size; k++)
	{
		result[k] = var1[k] * var2[k] * sqrt(max_size);
	}

	result = My_DFT(result).FastFourierTransform(-1);

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

		if (!in.eof()) //in.eof() - проверка конца файла
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
	std::numeric_limits< double > tmp;
	in.precision(std::streamsize(tmp.max_digits10));
	for (unsigned int i = 0; i < vector_r.size(); i++)
	{
		in << "    " << vector_r[i].real() << "    " << vector_r[i].imag() << std::endl;
	}

	in.close();
	return true;
}

bool WriteFile(std::string ss, std::string filename)
{
	std::ofstream in(filename);

	if (!in)
	{
		throw "File couldn't be opened!";
		return false;
	}

	in.precision(std::numeric_limits<double>::max_digits10);

	in << ss;

	in.close();
	return true;
}

bool Comparison(std::string file1, std::string file2, std::string fileName)
{
	/*
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
	return true;*/
	ComplexVect first = ReadFile(file1);
	ComplexVect second = ReadFile(file2);
	ComplexVect err(first.size(), 0);
	size_t size = first.size();
	for (size_t i = 0; i < size; i++)
	{
		err[i] = (first[i] - second[i]) * (first[i] - second[i]);
	}
	WriteFile(err, fileName);
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

	std::cout.precision(std::numeric_limits<double>::max_digits10);

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