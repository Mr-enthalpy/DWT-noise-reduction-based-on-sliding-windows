// Description: DWT

#include "DWT.h"

#include <cmath>
#include <string>
#include <iomanip>
#include <iostream>

template<typename T, typename THRMAKE>
Decompose<T, THRMAKE>::Decompose(std::vector<T> WAVE, Wavelet wlet, std::function<T(T, T)> func) :functhr(func), wavelet(wlet)
{
	if (!WAVE.size()) return;
	deConvolution(WAVE);
}

template<typename T, typename THRMAKE>
std::vector<T> Decompose<T, THRMAKE>::mra() const
{
	return approxCoeffs;
}

template<typename T, typename THRMAKE>
void Decompose<T, THRMAKE>::deConvolution(std::vector<T> WAVE)
{
	if (!WAVE.size()) return;
	size_t idx = 0;
	detailCoeffs.reserve((WAVE.size() + wavelet.size()) / 2 + ADD_CIR);
	approxCoeffs.reserve((WAVE.size() + wavelet.size()) / 2 + ADD_CIR);
	Rwave.reserve(WAVE.size() + ADD_CIR);
	detailCoeffs.resize((WAVE.size() + wavelet.size()) / 2, 0);
	approxCoeffs.resize((WAVE.size() + wavelet.size()) / 2, 0);
	Rwave.resize(WAVE.size(), 0);
	for (long long i = 0; i < detailCoeffs.size(); ++i)
	{
		for (long long j = 0; j < wavelet.size(); ++j)
		{
			idx = wave_edge(WAVE.size(), 2 * i - j);
			detailCoeffs[i] += WAVE[idx] * wavelet.HI_D[j];
			approxCoeffs[i] += WAVE[idx] * wavelet.LO_D[j];
		}
	}
	reflag = false;
}

template<typename T, typename THRMAKE>
std::vector<T> Decompose<T, THRMAKE>::RDWT(T thr)
{
	return RDWT(approxCoeffs, thr);
}

template<typename T, typename THRMAKE>
std::vector<T> Decompose<T, THRMAKE>::RDWT(std::vector<T> approxCoeffs, T thr)
{
	size_t idx = 0;
	if (reflag) return Rwave;
	for (size_t i = 0; i < Rwave.size(); ++i)
	{
		for (size_t j = 0; j < wavelet.size(); ++j)
		{
			idx = i - j + wavelet.size() - 1;
			if (idx % 2 == 0 && idx / 2 < detailCoeffs.size())
				Rwave[i] += functhr(detailCoeffs[idx / 2], thr) * wavelet.HI_R[j] + approxCoeffs[idx / 2] * wavelet.LO_R[j];
		}
	}
	reflag = true;
	return Rwave;
}

template<typename T, typename THRMAKE>
std::vector<T> Decompose<T, THRMAKE>::csignalband()
{
	if (!reflag) throw std::runtime_error("No RDWT");
	return Rwave;
}

template<typename T, typename THRMAKE>
void Decompose<T, THRMAKE>::resethr(std::function<T(T, T)> newfun)
{
	functhr = newfun;
}

template<typename T, typename THRMAKE>
void Decompose<T, THRMAKE>::resetwavelet(Wavelet wlet)
{
	wavelet = wlet;
}

template<typename T, typename THRMAKE>
DWT<T, THRMAKE>::DWT(std::vector<T>& WAVE, int J, Wavelet wlet, std::function<T(T, T)> func) :dwt(J), functhr(func), j(J), wavelet(wlet)
{
	if (!WAVE.size()) return;
	dwt[0].resetwavelet(wavelet);
	dwt[0].deConvolution(WAVE);
	dwt[0].resethr(functhr);
	for (long long i = 1; i < j; ++i)
	{
		dwt[i].resetwavelet(wavelet);
		dwt[i].deConvolution(dwt[i - 1].mra());
		dwt[i].resethr(functhr);
	}
	S = THRMAKE(*this);
	RDWT();
}

template<typename T, typename THRMAKE>
void DWT<T, THRMAKE>::reset(std::vector<T>& WAVE, int J, Wavelet wlet, std::function<T(T, T)> newfunc)
{
	if (!WAVE.size()) return;
	functhr = newfunc;
	wavelet = wlet;
	dwt.resize(j = J);
	dwt[0].resetwavelet(wavelet);
	dwt[0].deConvolution(WAVE);
	dwt[0].resethr(functhr);
	for (size_t i = 1; i < j; ++i)
	{
		dwt[i].resetwavelet(wavelet);
		dwt[i].deConvolution(dwt[i - 1].mra());
		dwt[i].resethr(functhr);
	}
	S = THRMAKE(*this);
	RDWT();
}

template<typename T, typename THRMAKE>
void DWT<T, THRMAKE>::RDWT()
{
	dwt[j - 1].RDWT(S.thrshold(j - 1));
	for (size_t i = j - 2; i; --i)
		dwt[i].RDWT(dwt[i + 1].csignalband(), S.thrshold(i));
	dwt[0].RDWT(dwt[1].csignalband(), S.thrshold(0));
}

template<typename T, typename THRMAKE>
void DWT<T, THRMAKE>::print(std::string prefilename)
{
	std::string filename;
	for (int i = 0; i < j; i++)
	{
		std::cout << "direct DWT thr:" << ' ' << std::setprecision(30) << S.thrshold(i) << '\n';
		filename = prefilename + "d" + std::to_string(i) + ".txt";
		std::ofstream outfile(filename);
		for (size_t j = 0; j < dwt[i].detailCoeffs.size(); ++j)
			outfile << dwt[i].detailCoeffs[j] << '\n';
		filename = prefilename + "a" + std::to_string(i) + ".txt";
		outfile.close();
		outfile.open(filename);
		for (size_t j = 0; j < dwt[i].approxCoeffs.size(); ++j)
			outfile << dwt[i].approxCoeffs[j] << '\n';
		outfile.close();
		filename = prefilename + "r" + std::to_string(i) + ".txt";
		outfile.open(filename);
		for (size_t j = 0; j < dwt[i].Rwave.size(); ++j)
			outfile << dwt[i].Rwave[j] << '\n';
		outfile.close();
	}
}

template<typename T, typename THRMAKE>
std::vector<T> DWT<T, THRMAKE>::csignal()
{
	return dwt[0].csignalband();
}

template<typename T, typename THRMAKE>
void DWT<T, THRMAKE>::relevel(int J)
{
	if (J <= j)
		dwt.resize(j = J);
	else
	{
		dwt.resize(J);
		for (size_t i = j; i < J; ++i)
			dwt[i].deConvolution(dwt[i - 1].mra());
		j = J;
		RDWT();
	}
}

template<typename T, typename THRMAKE>
size_t DWT<T, THRMAKE>::level() const
{
	return j;
}