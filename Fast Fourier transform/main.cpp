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
		DFT->WriteFile(DFT->DiscreteFourierTransform(-1), "DFT.txt");
		DFT->WriteFile(DFT->InverseDiscreteFourierTransform(1), "IDFT.txt");
		std::cout << "Calculating DFT and IDFT is over. You can see the results in DFT.txt, IDFT.txt" << std::endl;
	} catch (const char* ex)
	{
		std::cout << ex;
	}
	std::cout << "Calculating FFT and IFFT. Task 2" << std::endl;
	try
	{
		DFT = new My_DFT(fileName);
		DFT->WriteFile(DFT->FastFourierTransform(-1), "FFT.txt");
		DFT->WriteFile(DFT->InverseFastFourierTransform(1), "IFFT.txt");
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
}