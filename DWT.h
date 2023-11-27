#pragma once
#include <set>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>

#include "vectocirc.hpp"
#include "wavelet_filters.h"

static constexpr double Q4 = 0.6744897501960817;

template<typename T, typename THRMAKE>
class DWT;

template<typename T, typename THRMAKE>
class SDecompose;

template<typename T, typename THRMAKE, size_t stweep>
class SDWT;

template<typename T, typename THRMAKE>
class Decompose
{
private:
	friend class DWT<T, THRMAKE>;
	friend class SDecompose<T, THRMAKE>;
	friend THRMAKE;
	Wavelet wavelet;
	bool reflag = false;
	std::function<T(T, T)> functhr;
	std::vector<T> detailCoeffs, approxCoeffs, Rwave;
	size_t wave_edge(long long range, long long num)
	{
		if (num >= 0 && num < range)
			return num;
		if (num < 0)
			return  -num;
		return 2 * range - num - 2;
	}
public:
	Decompose() {}
	Decompose(std::vector<T> WAVE, Wavelet wlet, std::function<T(T, T)> func);
	std::vector<T> mra() const;
	void deConvolution(std::vector<T> WAVE);
	std::vector<T> RDWT(T thr);
	std::vector<T> RDWT(std::vector<T> approxCoeffs, T thr);
	std::vector<T> csignalband();
	void resethr(std::function<T(T, T)> newfun);
	void resetwavelet(Wavelet newavelet);
};

template<typename T, typename THRMAKE>
class DWT
{
private:
	THRMAKE S;
	size_t j = 3;
	Wavelet wavelet;
	std::function<T(T, T)> functhr;
	std::vector<Decompose<T, THRMAKE>> dwt;
	friend THRMAKE;
	template<typename T1, typename THRMAKE1, size_t stweep>
	friend class SDWT;
public:
	DWT() {};
	DWT(std::vector<T>& WAVE, int J, Wavelet wlet, std::function<T(T, T)> func);
	size_t level() const;
	void relevel(int J);
	void reset(std::vector<T>& WAVE, int J, Wavelet wlet, std::function<T(T, T)> newfunc);
	std::vector<T> csignal();
	void RDWT();
	void print(std::string prefilename = "s");
};

template<typename T, typename THRMAKE>
class SDecompose
{
private:
	Wavelet wavelet;
	std::function<T()> thrshold;
	std::function<T(T, T)> functhr;
	std::function<void(T)> insert, earse;
	size_t rescs = 0, mrastyle = 0, rers = 0;
	VectoCirc<T> detailCoeffs;
	TailVectoCirc<T> approxCoeffs, Rwave;
	size_t wave_edge(long long range, long long num)
	{
		if (num >= 0 && num < range)
			return num;
		if (num < 0)
			return  -num;
		return 2 * range - num - 2;
	}
	template<typename T1, typename THRMAKE1, size_t stweep>
	friend class SDWT;
public:
	SDecompose() {};
	SDecompose(Decompose<T, THRMAKE>& dwt, THRMAKE& s, size_t J, size_t reclac);
	SDecompose(Decompose<T, THRMAKE>&& dwt, THRMAKE& s, size_t J, size_t reclac);
	SDecompose& operator =(SDecompose&& other) noexcept;
	TailVectoCirc<T>& mra();
	void SweepConvoSub(TailVectoCirc<T>& WAVE);
	TailVectoCirc<T>& RSDWT();
	TailVectoCirc<T>& RSDWT(TailVectoCirc<T>& approxCoeffs);
	TailVectoCirc<T>& csignalband();
};

template<typename T, typename THRMAKE, size_t stweep>
class SDWT
{
private:
	THRMAKE S;
	Wavelet wavelet;
	size_t j, windowscale = 0;
	std::function<T(T, T)> functhr;
	std::vector<SDecompose<T, THRMAKE>> sdwt;
public:
	SDWT() :j(0) {}
	SDWT(SDWT& other) noexcept;
	SDWT(SDWT&& other) noexcept;
	SDWT& operator =(SDWT& other) noexcept;
	SDWT& operator =(SDWT&& other) noexcept;
	SDWT(std::vector<T>& WAVE, size_t J, Wavelet wlet, std::function<T(T, T)> func);
	int level() const;
	VectoCirc<T>& csignal();
	void slide(TailVectoCirc<T>& newave);
	void print(std::string prefilename = "o");
};

template<typename T>
class Dynamedian
{
private:
	struct cmp1
	{
		bool operator()(const T& a, const T& b) const
		{
			return abs(a) > abs(b);
		}
	};
	struct cmp2
	{
		bool operator()(const T& a, const T& b) const
		{
			return abs(a) < abs(b);
		}
	};
	std::multiset<T, cmp1> minelements;
	std::multiset<T, cmp2> maxelements;
	void balance()
	{
		if (maxelements.size() > minelements.size() + 1)
		{
			minelements.insert(*maxelements.begin());
			maxelements.erase(maxelements.begin());
		}
		else if (minelements.size() > maxelements.size())
		{
			maxelements.insert(*minelements.begin());
			minelements.erase(minelements.begin());
		}
	}
public:
	Dynamedian() {}
	Dynamedian(std::vector<T>& initial_elements);
	Dynamedian(VectoCirc<T>& initial_elements);
	Dynamedian& operator=(Dynamedian&&) noexcept;
	void earse(const T num);
	void insert(const T num);
	void clear();
	T median();
	~Dynamedian() {}
};

template<typename T>
class Dyvariance
{
private:
	T squaresum = 0, sum = 0;
	size_t _size = 0;
public:
	Dyvariance() {}
	Dyvariance(std::vector<T> initial);
	Dyvariance(VectoCirc<T> initial);
	Dyvariance& operator=(Dyvariance&&) noexcept;
	void earse(const T num);
	void insert(const T num);
	T variance();
	~Dyvariance() {}
};

template<typename T>
class Dynastein
{
private:
	struct cmp
	{
		bool operator()(const T& a, const T& b) const
		{
			return fabs(a) < fabs(b);
		}
	};
	std::multiset<T, cmp> elements;
public:
	Dynastein() {}
	Dynastein(std::vector<T>& initial_elements) :elements(initial_elements.begin(), initial_elements.end()) {}
	Dynastein(VectoCirc<T>& initial_elements) :elements(initial_elements.begin(), initial_elements.end()) {}
	Dynastein& operator=(Dynastein&&) noexcept;
	void earse(const T num);
	void insert(const T num);
	void clear();
	T thrsehold();
	~Dynastein() {}
};

template<typename T>
class BayesShrink
{
private:
	std::vector<Dyvariance<T>> parthr;
	Dynamedian<T> sigma;
	T sigma2()
	{
		return pow(sigma.median() / Q4, 2);
	}
public:
	BayesShrink() {}
	BayesShrink(DWT<T, BayesShrink<T>>& dwt);
	BayesShrink& operator=(BayesShrink&&) noexcept;
	void insert(size_t j, T num);
	void earse(size_t j, T num);
	T thrshold(const size_t j);
};

template<typename T>
class VisuShrink
{
private:
	std::vector<Dynamedian<T>> parthr;
	T sigma(size_t j)
	{
		return parthr[j].median() / Q4;
	}
	size_t N;
public:
	VisuShrink() {}
	VisuShrink(DWT<T, VisuShrink<T>>& dwt);
	VisuShrink& operator=(VisuShrink&&) noexcept;
	void insert(const size_t j, T num);
	void earse(const size_t j, T num);
	T thrshold(const size_t j);
};

template<typename T>
class Rigrsure
{
private:
	std::vector<Dynastein<T>> parthr;
public:
	Rigrsure() {}
	Rigrsure(DWT<T, Rigrsure<T>>& dwt);
	Rigrsure& operator=(Rigrsure&&) noexcept;
	void insert(const size_t j, T num);
	void earse(const size_t j, T num);
	T thrshold(const size_t j);
};