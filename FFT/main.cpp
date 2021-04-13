#include "DFT.h"

int main()
{
	std::cout << "Please, enter input file name" << std::endl;
	std::string fileName;
	std::cin >> fileName;
	std::cout << "Calculating DFT and IDFT. Task 1" << std::endl;
	My_DFT* DFT;
	try
	{
		DFT = new My_DFT(fileName);
		WriteFile(DFT->DiscreteFourierTransform(-1), "DFT.txt");
		WriteFile(DFT->InverseDiscreteFourierTransform(1), "IDFT.txt");
		std::cout << "Calculating DFT and IDFT is over. You can see the results in DFT.txt, IDFT.txt" << std::endl;
	}
	catch (const char* ex)
	{
		std::cout << ex;
	}

	std::cout << "Calculating FFT and IFFT. Task 2" << std::endl;
	try
	{
		DFT = new My_DFT(fileName);
		WriteFile(DFT->FastFourierTransform(1), "FFT.txt");
		WriteFile(DFT->InverseFastFourierTransform(-1), "IFFT.txt");
		std::cout << "Calculating FFT and IFFT is over. You can see the results in FFT.txt, IFFT.txt" << std::endl;
	}
	catch (const char* ex)
	{
		std::cout << ex;
	}

	try
	{
		std::cout << "Comparison input vector and IDFT. Task 3" << std::endl;
		if (Comparison(fileName, "IDFT.txt", "IDFT_err.txt"))
		{
			std::cout << "Comparison input vector and IDFT is successfully" << std::endl;
		}
		else
		{
			std::cout << "Comparison input vector and IDFT is unsuccessfully" << std::endl;
		}
		std::cout << "Comparison input vector and IFFT. Task 3" << std::endl;
		if (Comparison(fileName, "IFFT.txt", "IFFT_err.txt"))
		{
			std::cout << "Comparison input vector and IFFT is successfully" << std::endl;
		}
		else
		{
			std::cout << "Comparison input vector and IFFT is unsuccessfully" << std::endl;
		}
		std::cout << "Comparison DFT and FFT. Task 3" << std::endl;
		InaccuracyComparison("DFT.txt", "FFT.txt", "Inaccuracy.txt");
		std::cout << "Comparison DFT and FFT is over. You can see the results in Inaccuracy.txt" << std::endl;
	}
	catch (const char* ex)
	{
		std::cout << ex;
	}

	
	std::cout << "Calculating DFT and FFT performance. Task 4" << std::endl;
	try
	{
		std::cout << "Performance test on diferent vectors" << std::endl;
		ComplexVect tmp = ReadFile("test_1.txt");
		std::chrono::time_point <std::chrono::steady_clock> start, end;
		std::stringstream ss;
		ss.precision(std::numeric_limits<double>::max_digits10);
		for (int i = 1; i < 15; i++)
		{
			DFT = new My_DFT(ComplexVect(tmp.begin(), tmp.begin() + (1 << i)));
			start = std::chrono::steady_clock::now();
			DFT->DiscreteFourierTransform(-1);
			end = std::chrono::steady_clock::now();

			ss << "N = " << (1 << i) << " Time: " <<
				std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;
		}
		WriteFile(ss.str(), "DFT_perf_test.txt");
		ss.clear();
		for (int i = 1; i < 16; i++)
		{
			DFT = new My_DFT(ComplexVect(tmp.begin(), tmp.begin() + (1 << i)));
			start = std::chrono::steady_clock::now();
			DFT->FastFourierTransform(-1);
			end = std::chrono::steady_clock::now();

			ss << "N = " << (1 << i) << " Time: " <<
				std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;
		}
		WriteFile(ss.str(), "FFT_perf_test.txt");
		std::cout << "Calculating performance test is over. You can see the results in DFT_perf_test.txt, FFT_perf_test.txt.txt" << std::endl;
	}
	catch (const char* ex)
	{
		std::cout << ex;
	}

	std::cout << "Calculating Convolution. Task 5" << std::endl;
	try
	{
		DFT = new My_DFT(fileName);
		ComplexVect x, y;
		x = ReadFile("Convolution_x.txt");
		y = ReadFile("Convolution_y.txt");
		WriteFile(DFT->Convolution(x, y), "Convolution.txt");
		std::cout << "Calculating Convolution is over. You can see the results in Convolution.txt" << std::endl;
	}
	catch (const char* ex)
	{
		std::cout << ex;
	}

	std::cout << "Calculating ConvolutionFFT. Task 6" << std::endl;
	try
	{
		DFT = new My_DFT(fileName);
		ComplexVect x, y;
		x = ReadFile("Convolution_x.txt");
		y = ReadFile("Convolution_y.txt");
		WriteFile(DFT->ConvolutionFFT(x, y), "ConvolutionFFT.txt");
		std::cout << "Calculating ConvolutionFFT is over. You can see the results in ConvolutionFFT.txt" << std::endl;
	}
	catch (const char* ex)
	{
		std::cout << ex;
	}

	std::cout << "Calculating convolution and convolutionFFT performance. Task 8" << std::endl;
	try
	{
		ComplexVect x, y;
		std::chrono::time_point <std::chrono::steady_clock> start1, end1, start2, end2;

		long long m_fft[10], m_ft[10];
		std::vector<long long> ms_fft(10), ms_dft(10);

		x = ReadFile("Convolution_x.txt");
		y = ReadFile("Convolution_y.txt");

		std::cout << "One vector with fixed length:" << std::endl;
		std::cout << "ConvFFT\tConv" << std::endl;
		for (int j = 0; j < 2; j++)
		{
			ComplexVect vector_x1(x.begin(), x.begin() + pow(2, j)), vector_x2(y.begin(), y.begin() + pow(2, 3));

			int result_size = vector_x1.size() + vector_x2.size() - 1;
			int max_size = vector_x1.size() > vector_x2.size() ? vector_x1.size() : vector_x2.size();
			for (int i = vector_x1.size(); i < 2 * max_size; i++)
			{
				vector_x1.push_back(0);
			}
			for (int i = vector_x2.size(); i < 2 * max_size; i++)
			{
				vector_x2.push_back(0);
			}
		  // auto var1 = My_DFT(vector_x1).FastFourierTransform(1);
			//auto var2 = My_DFT(vector_x2).FastFourierTransform(1);
			auto var1 = My_DFT(vector_x1).FastFourierTransform(-1);
			auto var2 = My_DFT(vector_x2).FastFourierTransform(-1);
			start1 = std::chrono::steady_clock::now();
			My_DFT(vector_x1).ConvolutionFFT(var1, var2, result_size, max_size);
			end1 = std::chrono::steady_clock::now();
			m_fft[j] = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - start1).count();

			start2 = std::chrono::steady_clock::now();
			My_DFT(vector_x1).Convolution(vector_x1, vector_x2);
			end2 = std::chrono::steady_clock::now();
			m_ft[j] = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - start2).count();
			std::cout << m_fft[j] << "\t" << m_ft[j] << std::endl;
		}

		std::cout << "Both vectors change the lengths:" << std::endl;
		std::cout << "ConvFFT\tConv" << std::endl;
		for (int j = 0; j < 4; j++)
		{
			ComplexVect vector_x1(x.begin(), x.begin() + pow(2, j)), vector_x2(y.begin(), y.begin() + pow(2, j));

			int result_size = vector_x1.size() + vector_x2.size() - 1;
			int max_size = vector_x1.size() > vector_x2.size() ? vector_x1.size() : vector_x2.size();
			for (int i = vector_x1.size(); i < 2 * max_size; i++)
			{
				vector_x1.push_back(0);
			}
			for (int i = vector_x2.size(); i < 2 * max_size; i++)
			{
				vector_x2.push_back(0);
			}

			//auto var1 = My_DFT(vector_x1).FastFourierTransform(1);
			//auto var2 = My_DFT(vector_x2).FastFourierTransform(1);
			auto var1 = My_DFT(vector_x1).FastFourierTransform(-1);
			auto var2 = My_DFT(vector_x2).FastFourierTransform(-1);

			start1 = std::chrono::steady_clock::now();
			My_DFT(vector_x1).ConvolutionFFT(var1, var2, result_size, max_size);
			end1 = std::chrono::steady_clock::now();
			m_fft[j] = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - start1).count();

			start2 = std::chrono::steady_clock::now();
			My_DFT(vector_x1).Convolution(vector_x1, vector_x2);;
			end2 = std::chrono::steady_clock::now();
			m_ft[j] = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - start2).count();
			std::cout << m_fft[j] << "\t" << m_ft[j] << std::endl;
		}
	}
	catch (const char* ex)
	{
		std::cout << ex;
	}
}

