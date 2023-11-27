// SDWT.cpp: 定义应用程序的入口点。
//
#include "DWT.h"


template<typename T, typename THRMAKE>
SDecompose<T, THRMAKE>::SDecompose(Decompose<T, THRMAKE>& DWT, THRMAKE& s, size_t J, size_t reclac) : detailCoeffs(DWT.detailCoeffs), approxCoeffs(DWT.approxCoeffs), functhr(DWT.functhr), Rwave(DWT.Rwave), wavelet(DWT.wavelet)
{
	rescs = reclac;
	if (Rwave.size() % 2) ++rescs;
	thrshold = std::bind(&THRMAKE::thrshold, &s, J);
	insert = std::bind(&THRMAKE::insert, &s, J, std::placeholders::_1);
	earse = std::bind(&THRMAKE::earse, &s, J, std::placeholders::_1);
};

template<typename T, typename THRMAKE>
SDecompose<T, THRMAKE>::SDecompose(Decompose<T, THRMAKE>&& DWT, THRMAKE& s, size_t J, size_t reclac) : detailCoeffs(std::move(DWT.detailCoeffs)), approxCoeffs(std::move(DWT.approxCoeffs)), Rwave(std::move(DWT.Rwave)), functhr(DWT.functhr), wavelet(DWT.wavelet)
{
	rescs = reclac;
	if (Rwave.size() % 2) ++rescs;
	thrshold = std::bind(&THRMAKE::thrshold, &s, J);
	insert = std::bind(&THRMAKE::insert, &s, J, std::placeholders::_1);
	earse = std::bind(&THRMAKE::earse, &s, J, std::placeholders::_1);
	DWT.detailCoeffs.clear();
	DWT.approxCoeffs.clear();
	DWT.Rwave.clear();
	DWT.functhr = nullptr;
}

template<typename T, typename THRMAKE>
SDecompose<T, THRMAKE>& SDecompose<T, THRMAKE>::operator =(SDecompose&& other) noexcept
{
	using std::move;
	detailCoeffs = move(other.detailCoeffs);
	approxCoeffs = TailVectoCirc<T>(move(other.approxCoeffs));
	Rwave = TailVectoCirc<T>(move(other.Rwave));
	wavelet = other.wavelet;
	functhr = move(other.functhr);
	thrshold = move(other.thrshold);
	insert = move(other.insert);
	earse = move(other.earse);
	mrastyle = other.mrastyle;
	rescs = other.rescs;
	rers = other.rers;
	other.functhr = nullptr;
	other.thrshold = nullptr;
	other.earse = nullptr;
	other.insert = nullptr;
	other.mrastyle = other.rescs = other.rers = 0;
	other.detailCoeffs.clear();
	other.approxCoeffs.clear();
	other.Rwave.clear();
	return *this;
}

template<typename T, typename THRMAKE>
TailVectoCirc<T>& SDecompose<T, THRMAKE>::mra()
{
	return approxCoeffs;
}

template<typename T, typename THRMAKE>
void SDecompose<T, THRMAKE>::SweepConvoSub(TailVectoCirc<T>& WAVE)
{
	if (!WAVE.size()) return;
	bool recflag = true;
	size_t idx, sfront = (WAVE.frontmove() + 1 - mrastyle) / 2;
	if (WAVE.frontmove() % 2)
		mrastyle = 1 - mrastyle;
	for (size_t i = rescs; i < detailCoeffs.size(); ++i)
	{
		earse(detailCoeffs[i]);
		approxCoeffs[i] = detailCoeffs[i] = 0;
	}
	for (size_t i = 0; i < sfront; ++i)
		earse(detailCoeffs.pop_front());
	approxCoeffs.shrink_front(sfront);
	size_t srescs = rescs - sfront;
	Rwave.shrink_front(WAVE.frontmove());
	Rwave.extend_back(WAVE.endmove());
	auto sweep = [&](auto& Coeffs)->void
		{
			if ((WAVE.size() - WAVE.endmove() + wavelet.size() - 1) % 2)
				Coeffs.extend_back((WAVE.endmove() + mrastyle) / 2);
			else Coeffs.extend_back((WAVE.endmove() + 1 - mrastyle) / 2);
		};
	sweep(detailCoeffs);
	sweep(approxCoeffs);
	WAVE.Strelief();
	detailCoeffs.Strelief();
	for (long long i = srescs; i < detailCoeffs.size(); ++i)
	{
		for (long long j = 0; j < wavelet.size(); ++j)
		{
			idx = wave_edge(WAVE.size(), 2 * i - j + mrastyle);
			detailCoeffs[i] += WAVE[idx] * wavelet.HI_D[j];
			approxCoeffs[i] += WAVE[idx] * wavelet.LO_D[j];
		}
		insert(detailCoeffs[i]);
	}
	rers = 2 * srescs + mrastyle;
}

template<typename T, typename THRMAKE>
TailVectoCirc<T>& SDecompose<T, THRMAKE>::RSDWT()
{
	return RSDWT(approxCoeffs);
}

template<typename T, typename THRMAKE>
TailVectoCirc<T>& SDecompose<T, THRMAKE>::RSDWT(TailVectoCirc<T>& approxCoeffs)
{
	size_t idx = 0;
	T thr = thrshold();
	std::fill(Rwave.begin() + rers, Rwave.end(), 0);
	for (size_t i = rers; i < Rwave.size(); ++i)
	{
		for (size_t j = 0; j < wavelet.size(); ++j)
		{
			idx = i - j + wavelet.size() - 1;
			if (idx % 2 == mrastyle && idx / 2 < detailCoeffs.size())
				Rwave[i] += functhr(detailCoeffs[idx / 2], thr) * wavelet.HI_R[j] + approxCoeffs[idx / 2] * wavelet.LO_R[j];
		}
	}
	Rwave.Strelief();
	return Rwave;
}

template<typename T, typename THRMAKE>
TailVectoCirc<T>& SDecompose<T, THRMAKE>::csignalband()
{
	return Rwave;
}

template<typename T, typename THRMAKE, size_t stweep>
SDWT<T, THRMAKE, stweep>::SDWT(SDWT& other) noexcept
{
	j = other.j;
	sdwt = other.sdwt;
	functhr = other.functhr;
	wavelet = other.wavelet;
	S = other.S;
}

template<typename T, typename THRMAKE, size_t stweep>
SDWT<T, THRMAKE, stweep>::SDWT(SDWT&& other) noexcept
{
	using std::move;
	j = other.j;
	sdwt = move(other.sdwt);
	functhr = other.functhr;
	wavelet = other.wavelet;
	S = other.S;
	other.wave.clear();
	other.sdwt.clear();
	other.functhr = nullptr;
	other.S = nullptr;
}

template<typename T, typename THRMAKE, size_t stweep>
SDWT<T, THRMAKE, stweep>& SDWT<T, THRMAKE, stweep>::operator =(SDWT& other) noexcept
{
	j = other.j;
	sdwt = other.sdwt;
	functhr = other.functhr;
	wavelet = other.wavelet;
	S = other.S;
	return *this;
}

template<typename T, typename THRMAKE, size_t stweep>
SDWT<T, THRMAKE, stweep>& SDWT<T, THRMAKE, stweep>::operator =(SDWT&& other) noexcept
{
	using std::move;
	j = other.j;
	sdwt = move(other.sdwt);
	S = move(other.S);
	functhr = other.functhr;
	wavelet = other.wavelet;
	other.wave.clear();
	other.sdwt.clear();
	other.functhr = nullptr;
	return *this;
}

template<typename T, typename THRMAKE, size_t stweep>
SDWT<T, THRMAKE, stweep>::SDWT(std::vector<T>& WAVE, size_t J, Wavelet wlet, std::function<T(T, T)> func) : j(J), sdwt(J), functhr(func), wavelet(wlet)
{
	if (stweep > WAVE.size() / 3)
		throw std::invalid_argument("stweep is too large");
	DWT<T, THRMAKE> dwtemp(WAVE, j, wavelet, functhr);
	std::vector<size_t> max_tail_clac(j);
	max_tail_clac[0] = (stweep + 1) / 2;
	for (size_t i = 1; i < j; ++i)
		max_tail_clac[i] = (max_tail_clac[i - 1] + 1) / 2;
	S = std::move(dwtemp.S);
	size_t reclac = (WAVE.size() + 1) / 2;
	for (size_t i = 0; i < j; ++i)
	{
		sdwt[i] = SDecompose<T, THRMAKE>(std::move(dwtemp.dwt[i]), S, i, reclac);
		sdwt[i].approxCoeffs.virtual_transform(reclac - max_tail_clac[i]);
		if (i < j - 1) reclac = (reclac + 1) / 2;
	}
	max_tail_clac[j - 1] = 2 * (reclac - max_tail_clac[j - 1]) - 2;
	for (size_t i = j - 1; i; --i)
	{
		sdwt[i].Rwave.virtual_transform(sdwt[i].Rwave.size() - max_tail_clac[i]);
		max_tail_clac[i - 1] = 2 * max_tail_clac[i] - wavelet.size();
	}
}

template<typename T, typename THRMAKE, size_t stweep>
int SDWT<T, THRMAKE, stweep>::level() const
{
	return j;
}

template<typename T, typename THRMAKE, size_t stweep>
VectoCirc<T>& SDWT<T, THRMAKE, stweep>::csignal()
{
	return sdwt[0].csignalband();
}

template<typename T, typename THRMAKE, size_t stweep>
void SDWT<T, THRMAKE, stweep>::slide(TailVectoCirc<T>& newave)
{
	if (newave.endmove() > stweep)
		throw std::invalid_argument("newave is too large");
	sdwt[0].SweepConvoSub(newave);
	for (size_t i = 1; i < j; ++i)
		sdwt[i].SweepConvoSub(sdwt[i - 1].mra());
	sdwt[j - 1].RSDWT();
	for (size_t i = j - 2; i; --i)
	{
		sdwt[i].rers = 2 * sdwt[i + 1].rers - wavelet.size();
		sdwt[i].RSDWT(sdwt[i + 1].csignalband());
	}
	sdwt[0].rers = 2 * sdwt[1].rers - wavelet.size();
	sdwt[0].RSDWT(sdwt[1].csignalband());
}

template<typename T, typename THRMAKE, size_t stweep>
void SDWT<T, THRMAKE, stweep>::print(std::string prefilename)
{
	std::string filename;
	for (int i = 0; i < j; i++)
	{
		std::cout << "test SDWT thr:" << ' ' << S.thrshold(i) << '\n';
		filename = prefilename + "d" + std::to_string(i) + ".txt";
		std::ofstream outfile(filename);
		for (size_t j = 0; j < sdwt[i].detailCoeffs.size(); ++j)
			outfile << sdwt[i].detailCoeffs[j] << '\n';
		filename = prefilename + "a" + std::to_string(i) + ".txt";
		outfile.close();
		outfile.open(filename);
		for (size_t j = 0; j < sdwt[i].approxCoeffs.size(); ++j)
			outfile << sdwt[i].approxCoeffs[j] << '\n';
		outfile.close();
		filename = prefilename + "r" + std::to_string(i) + ".txt";
		outfile.open(filename);
		for (size_t j = 0; j < sdwt[i].Rwave.size(); ++j)
			outfile << sdwt[i].Rwave[j] << '\n';
		outfile.close();
	}
}
/*
class Sdenoeval//降噪效果评价
{
private:
	vector<double> csignal, signal, dsignal;//csignal为干净信号，signal为含噪信号,dsignal为降噪信号
	int J, N;//分解层数，信号总长度
	double C;//阈值调参
	SDWT dwt;
public:
	Sdenoeval(vector<double> Csignal, vector<double> Signal, int n, int j, double c) :N(n), J(j), C(c), dwt(signal, J, C)
	{
		csignal = Csignal;
		signal = Signal;
		dsignal = dwt.csignal();
	};
	Sdenoeval()
	{
		csignal = vector<double>{ 0 };
		signal = vector<double>{ 0 };
		dwt = SDWT();
		N = J = 0;
		C = 1.0;
		dsignal = vector<double>{ 0 };
	}
	double SNR()//信噪比
	{
		double sum = 0, sum1 = 0;
		for (int i = 0; i < csignal.size(); ++i)
		{
			sum += csignal[i] * csignal[i];
			sum1 += (csignal[i] - dsignal[i]) * (csignal[i] - dsignal[i]);
		}
		return 10 * log10(sum / sum1);
	}
	double MSE()//均方误差
	{
		double sum = 0;
		for (int i = 0; i < csignal.size(); ++i)
			sum += (csignal[i] - dsignal[i]) * (csignal[i] - dsignal[i]);
		return sum / csignal.size();
	}
	double SD()//信号失真度
	{
		double sum = 0;
		for (int i = 0; i < csignal.size(); ++i)
			sum += (csignal[i] - dsignal[i]) * (csignal[i] - dsignal[i]) / csignal[i] * csignal[i];
		return sum;
	}
	double Timcomplex() const
	{
		return (2 - pow(0.5, J)) * N;
	}
};*/

