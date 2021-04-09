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
	} catch (const char* ex)
	{
		std::cout << ex;
	}

	std::cout << "Calculating FFT and IFFT. Task 2" << std::endl;
	try
	{
		DFT = new My_DFT(fileName);
		WriteFile(DFT->FastFourierTransform(-1), "FFT.txt");
		WriteFile(DFT->InverseFastFourierTransform(1), "IFFT.txt");
		std::cout << "Calculating FFT and IFFT is over. You can see the results in FFT.txt, IFFT.txt" << std::endl;
	}
	catch (const char* ex)
	{
		std::cout << ex;
	}
	
	try
	{
		std::cout << "Comparison input vector and IDFT. Task 3" << std::endl;
		if (Comparison(fileName, "IDFT.txt"))
		{
			std::cout << "Comparison input vector and IDFT is successfully" << std::endl;
		} else 
		{
			std::cout << "Comparison input vector and IDFT is unsuccessfully" << std::endl;
		}
		std::cout << "Comparison input vector and IFFT. Task 3" << std::endl;
		if (Comparison(fileName, "IFFT.txt"))
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
		std::cout << "Calculating FFT and IFFT is over. You can see the results in FFT.txt, IFFT.txt" << std::endl;
	}
	catch (const char* ex)
	{
		std::cout << ex;
	}
}