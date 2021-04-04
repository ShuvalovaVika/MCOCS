#include "DFT.h"

My_DFT::My_DFT(std::string fileName)
{
	ReadFile(fileName);
	y.resize(x.size(), 0.0);
}

inline complex_temp My_DFT::omega(int k, int l, int sign)
{
	return exp(complex_temp(0, (double)sign * 2.0 * M_PI * (double)l / (double)(pow(2., k))));
}

inline complex_temp My_DFT::exponent(int N, int j, int k, int sign)
{
	return exp(complex_temp(0, (double)sign * 2 * M_PI * (double)j * (double)k / ((double)N)));
}

ComplexVect My_DFT::DiscreteFourierTransform(int sign) // DPF formula 1
{
	size_t N = x.size();
	for (int k = 0; k < N; k++)
	{
		for (int j = 0; j < N; j++)
		{
			y[k] = y[k] + x[j] * exponent(N, j, k, sign);
		}

		y[k] = y[k] / sqrt(N);
	}

	return y;
}

ComplexVect My_DFT::InverseDiscreteFourierTransform(int sign) //ОДПФ 2 формула
{
	size_t N = y.size();
	x.clear();
	for (int i = 0; i < N; i++)
	{
		x.push_back(0);
	}
	for (int j = 0; j < N; j++)
	{
		for (int k = 0; k < N; k++)
		{
			x[j] = x[j] + y[k] * exponent(N, j, k, sign);
		}
		x[j] = x[j] / sqrt(N);
	}
	
	return x;
}

ComplexVect My_DFT::FastFourierTransform(int sign)
{
	int N = x.size();
	int n = log2(N);

	for (int k = 1; k < n + 1; k++)
	{
		for (int j = 0; j < (1 << (k - 1)); j++)
		{
			for (int l = 0; l < (1 << (n - k)); l++)
			{
				int even = j * (1 << (n + 1 - k)) + l;
				int odd = even + (1 << (n - k));
				y[even] = x[even] + x[odd];
				y[odd] = (x[even] - x[odd]) * omega(n + 1 - k, l, sign);
			}
		}
	}

	ComplexVect z(y);

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
		z[k] = y[inverse];
	}

	y = z;

	for (size_t i = 0; i < N; i++)
		y[i] /= (double)sqrt(N);

	return y;
}

ComplexVect My_DFT::InverseFastFourierTransform(int sign)
{
	int N = y.size();
	int n = log2(N);

	for (int k = 1; k < n + 1; k++)
	{
		for (int j = 0; j < (1 << (k - 1)); j++)
		{
			for (int l = 0; l < (1 << (n - k)); l++)
			{
				int even = j * (1 << (n + 1 - k)) + l;
				int odd = even + (1 << (n - k));
				x[even] = y[even] + y[odd];
				x[odd] = (y[even] - y[odd]) * omega(n + 1 - k, l, sign);
			}
		}
	}

	ComplexVect z(x);

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
		z[k] = x[inverse];
	}

	x = z;

	for (size_t i = 0; i < N; i++)
		x[i] /= (double)sqrt(N);// x[i] = x[i] / sqrt

	return x;
}

ComplexVect My_DFT::Convolution(ComplexVect vector_x1, ComplexVect vector_x2)
{
	int size_x1 = vector_x1.size();
	int size_x2 = vector_x2.size();
	int max_size = std::max(size_x1, size_x2);

	ComplexVect var1, var2;
	ComplexVect vector_y(2 * max_size);

	var1 = vector_x1;
	var2 = vector_x2;
	var1.resize(2 * max_size);
	var2.resize(2 * max_size);

	for (int i = 0; i < 2 * max_size; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			if (i - j >= 0)
				vector_y[i] = vector_y[i] + var1[j] * var2[i - j];
		}
	}

	return vector_y;
}

ComplexVect My_DFT::ConvolutionFFT(ComplexVect vector_x1, ComplexVect vector_x2)
{
	int size_x1 = vector_x1.size();
	int size_x2 = vector_x2.size();
	int max_size = max(size_x1, size_x2);

	ComplexVect var1, var2;

	ComplexVect vector_y(2 * max_size);

	var1 = vector_x1;
	var2 = vector_x2;
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
			vector_y.resize(max_size);
			break;
		}

		count++;
	}

	var1 = FastFourierTransform(var1, -1);
	var2 = FastFourierTransform(var2, -1);

	for (unsigned int k = 0; k < max_size; k++)
	{
		vector_y[k] = var1[k] * var2[k] * sqrt(max_size);
	}

	vector_y = FastFourierTransform(vector_y, 1);

	return vector_y;
}

ComplexVect My_DFT::ConvolutionFFT(ComplexVect vector_x1, ComplexVect vector_x2, int result_size, int max_size)
{
	ComplexVect vector_y;
	for (int i = 0; i < 2 * max_size; i++)
	{
		vector_y.push_back(sqrt(2 * max_size) * vector_x1[i] * vector_x2[i]);
	}

	auto vector_result = FastFourierTransform(vector_y, -1);

	for (int i = 0; i < 2 * max_size - result_size; i++)
	{
		vector_result.pop_back();
	}

	return vector_result;
}

void My_DFT::ReadFile(std::string filename)
{
	std::ifstream in(filename);
	if (in.is_open())
	{
		double var2, var3;

		if (!in.eof()) //in.eof() - проверка конца файла
			in >> var2 >> var3;

		while (!in.eof())
		{
			x.push_back(complex_temp(var2, var3));
			in >> var2 >> var3;
		}

		y.resize(x.size());

		in.close();
	}
	else
	{
		throw "File couldn't be opened!";
	}
}

bool My_DFT::WriteFile(ComplexVect vector_r, std::string filename)
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
	ComplexVect first, second;
	if (in.is_open())
	{
		double var2, var3;

		if (!in.eof()) //in.eof() - проверка конца файла
			in >> var2 >> var3;

		while (!in.eof())
		{
			first.push_back(complex_temp(var2, var3));
			in >> var2 >> var3;
		}

		in.close();
	}
	else
	{
		throw "File couldn't be opened!";
	}
	in.open(file2);
	if (in.is_open())
	{
		double var2, var3;

		if (!in.eof()) //in.eof() - проверка конца файла
			in >> var2 >> var3;

		while (!in.eof())
		{
			second.push_back(complex_temp(var2, var3));
			in >> var2 >> var3;
		}

		in.close();
	}
	else
	{
		throw "File couldn't be opened!";
	}
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
	ComplexVect first, second;
	if (in.is_open())
	{
		double var2, var3;

		if (!in.eof()) //in.eof() - проверка конца файла
			in >> var2 >> var3;

		while (!in.eof())
		{
			first.push_back(complex_temp(var2, var3));
			in >> var2 >> var3;
		}

		in.close();
	}
	else
	{
		throw "File couldn't be opened!";
	}
	in.open(file2);
	if (in.is_open())
	{
		double var2, var3;

		if (!in.eof()) //in.eof() - проверка конца файла
			in >> var2 >> var3;

		while (!in.eof())
		{
			second.push_back(complex_temp(var2, var3));
			in >> var2 >> var3;
		}

		in.close();
	}
	else
	{
		throw "File couldn't be opened!";
	}
	size_t size = first.size();
	ComplexVect result(second);
	for (size_t i = 0; i < size; i++)
	{
		result[i] -= first[i];
	}

	auto max = max_element(result.begin(), result.end(), [](complex_temp v1, complex_temp v2) {return abs(v1) < abs(v2); });

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
