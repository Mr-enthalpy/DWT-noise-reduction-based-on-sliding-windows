#include <limits>
#include <iostream>

#include "DWT.h"

template<typename T>
Dynamedian<T>::Dynamedian(std::vector<T>& initial_elements)
{
	for (auto& i : initial_elements)
		insert(i);
}

template<typename T>
Dynamedian<T>::Dynamedian(VectoCirc<T>& initial_elements)
{
	for (size_t i = 0; i < initial_elements.size(); ++i)
		insert(initial_elements[i]);
}

template<typename T>
Dynamedian<T>& Dynamedian<T>::operator=(Dynamedian&& other) noexcept
{
	maxelements = std::move(other.maxelements);
	minelements = std::move(other.minelements);
	other.maxelements.clear();
	other.minelements.clear();
	return *this;
}

template<typename T>
void Dynamedian<T>::earse(const T num)
{
	if (!maxelements.empty() && abs(num) >= abs(*maxelements.begin()))
	{
		auto it = maxelements.find(num);
		if (it != maxelements.end())
			maxelements.erase(it);
		else throw std::exception("Element not found");
	}
	else if (!minelements.empty() && abs(num) <= abs(*minelements.begin()))
	{
		auto it = minelements.find(num);
		if (it != minelements.end())
			minelements.erase(it);
		else throw std::exception("Element not found");
	}
	else throw std::exception("Element not found");
	balance();
}

template<typename T>
void Dynamedian<T>::insert(const T num)
{
	if (maxelements.empty() || abs(num) >= abs(*maxelements.begin()))
		maxelements.insert(num);
	else minelements.insert(num);
	balance();
}

template<typename T>
void Dynamedian<T>::clear()
{
	maxelements.clear();
	minelements.clear();
}

template<typename T>
T Dynamedian<T>::median()
{
	if (maxelements.size() == minelements.size())
		return (abs(*maxelements.begin()) + abs(*minelements.begin())) / 2.0;
	return abs(*maxelements.begin());
}

template<typename T>
Dyvariance<T>::Dyvariance(std::vector<T> initial) :_size(initial.size())
{
	for (auto i : initial)
	{
		squaresum += pow(i, 2);
		sum += i;
	}
}


template<typename T>
Dyvariance<T>::Dyvariance(VectoCirc<T> initial) :_size(initial.size())
{
	for(size_t i=0;i<initial.size();++i)
	{
		squaresum += pow(initial[i], 2);
		sum += initial[i];
	}
}

template<typename T>
Dyvariance<T>& Dyvariance<T>::operator=(Dyvariance&& other) noexcept
{
	_size = other._size;
	squaresum = other.squaresum;
	sum = other.sum;
	other._size = 0;
	other.squaresum = 0;
	other.sum = 0;
	return *this;
}

template<typename T>
void Dyvariance<T>::earse(const T num)
{
	--_size;
	squaresum -= num * num;
	sum -= num;
}

template<typename T>
void Dyvariance<T>::insert(const T num)
{
	++_size;
	squaresum += num * num;
	sum += num;
}

template<typename T>
T Dyvariance<T>::variance()
{
	return squaresum / _size - pow((sum / _size), 2);
}

template<typename T>
Dynastein<T>& Dynastein<T>::operator=(Dynastein&& other) noexcept
{
	elements = std::move(other.elements);
	other.elements.clear();
	return *this;
}

template<typename T>
void Dynastein<T>::earse(const T num)
{
	elements.erase(elements.find(num));
}

template<typename T>
void Dynastein<T>::insert(const T num)
{
	elements.insert(num);
}

template<typename T>
void Dynastein<T>::clear()
{
	elements.clear();
}

template<typename T>
T Dynastein<T>::thrsehold()
{
	T sum, risk_min, Si, risk, thr = 0.0;
	size_t num = 0;
	auto it = elements.begin();
	sum = Si = pow((*it), 2);
	risk_min = 1 + (sum - 2 + (elements.size() - 2) * Si) / elements.size();
	thr = abs(*it);
	for (size_t i = 1; it != elements.end(); ++it, ++i)
	{
		Si = pow((*it), 2);
		sum += Si;
		risk = 1 + (sum - 2 * (i + 1) + (elements.size() - i - 2) * Si) / elements.size();
		if (risk_min > risk)
		{
			num = i;
			risk_min = risk;
			thr = abs(*it);
		}
	}
	return thr;
}

template<typename T>
BayesShrink<T>::BayesShrink(DWT<T, BayesShrink<T>>& dwt) :parthr(dwt.j), sigma(dwt.dwt[dwt.j - 1].detailCoeffs)
{
	for (size_t i = 0; i < parthr.size(); ++i)
		parthr[i] = Dyvariance(dwt.dwt[i].detailCoeffs);
}


template<typename T>
BayesShrink<T>& BayesShrink<T>::operator=(BayesShrink&& other) noexcept
{
	parthr = std::move(other.parthr);
	sigma = std::move(other.sigma);
	other.parthr.clear();
	other.sigma.clear();
	return *this;
}


template<typename T>
void BayesShrink<T>::insert(const size_t j, T num)
{
	parthr[j].insert(num);
	if (j == parthr.size() - 1)
		sigma.insert(num);
}

template<typename T>
void BayesShrink<T>::earse(const size_t j, T num)
{
	parthr[j].earse(num);
	if (j == parthr.size() - 1)
		sigma.earse(num);
}

template<typename T>
T BayesShrink<T>::thrshold(const size_t j)
{
	double sig1 = sigma2(), sig2 = parthr[j].variance();
	if (sig1 > sig2)
		return std::numeric_limits<T>::max();
	return sqrt(2) * sig1 / sqrt(sig2 - sig1);
}

template<typename T>
VisuShrink<T>::VisuShrink(DWT<T, VisuShrink<T>>& dwtemp) :parthr(dwtemp.j), N(dwtemp.dwt[0].detailCoeffs.size())
{
	for (size_t i = 0; i < parthr.size(); ++i)
		parthr[i] = Dynamedian(dwtemp.dwt[i].detailCoeffs);
}

template<typename T>
VisuShrink<T>& VisuShrink<T>::operator=(VisuShrink&& other) noexcept
{
	parthr = std::move(other.parthr);
	N = other.N;
	other.parthr.clear();
	return *this;
}

template<typename T>
void VisuShrink<T>::insert(const size_t j, T num)
{
	parthr[j].insert(num);
}

template<typename T>
void VisuShrink<T>::earse(const size_t j, T num)
{
	parthr[j].earse(num);
}

template<typename T>
T VisuShrink<T>::thrshold(const size_t j)
{
	double see = sqrt(2 * log(N)) / log(j + 2);
	return sigma(j) * sqrt(2 * log(N)) / log(j + 2);
}

template<typename T>
Rigrsure<T>::Rigrsure(DWT<T, Rigrsure<T>>& dwt) :parthr(dwt.j)
{
	for (size_t i = 0; i < parthr.size(); ++i)
		parthr[i] = Dynastein(dwt.dwt[i].detailCoeffs);
}

template<typename T>
Rigrsure<T>& Rigrsure<T>::operator=(Rigrsure&& other) noexcept
{
	parthr = std::move(other.parthr);
	other.parthr.clear();
	return *this;
}

template<typename T>
void Rigrsure<T>::insert(const size_t j, T num)
{
	parthr[j].insert(num);
}

template<typename T>
void Rigrsure<T>::earse(const size_t j, T num)
{
	parthr[j].earse(num);
}

template<typename T>
T Rigrsure<T>::thrshold(const size_t j)
{
	return parthr[j].thrsehold();
}