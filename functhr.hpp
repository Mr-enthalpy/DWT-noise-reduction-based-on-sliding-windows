#include <cmath>

#include "functhr.h"

template<typename T>
T softhresh(T si, T thr)
{
	if (abs(si) < thr)
		return 0.0;
	return si - sign(si) * thr;
}

template<typename T>
T threshard(T si, T thr)
{
	if (abs(si) < thr)
		return 0.0;
	return si;
}

template<typename T>
T Semisofthr(T si, T thr)
{
	if (abs(si) < thr)
		return 0.0;
	return si - sign(si) * thr / 2;
}

template<typename T>
T comprothr(T si, T thr)
{
	if (abs(si) < thr)
		return 0.0;
	return si - sign(si) * thr / exp(thr * (abs(si) - thr));
}

template<typename T>
T Mosquare(T si, T thr)
{
	if (abs(si) < thr)
		return 0.0;
	return sign(si) * sqrt(pow(si, 2) - pow(thr, 2));
}