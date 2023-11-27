#pragma once

template<typename T>
int sign(T val)
{
	return (T(0) < val) - (val < T(0));
}

template<typename T>
T softhresh(T si, T thr);

template<typename T>
T threshard(T si, T thr);

template<typename T>
T Semisofthr(T si, T thr);

template<typename T>
T comprothr(T si, T thr);

template<typename T>
T Mosquare(T si, T thr);
