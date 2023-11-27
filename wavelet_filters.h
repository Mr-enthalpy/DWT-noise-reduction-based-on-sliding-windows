#pragma once
#include <span>

enum class waveletype
{
	db1, db2, db3, db4, db5
	, db6, db7, db8, db9, db10
	, db11, db12, db13, db14, db15
	, db16, db17, db18, db19, db20
	, db21, db22, db23, db24, db25
	, db26, db27, db28, db29, db30
	, db31, db32, db33, db34, db35
	, db36, db37, db38, db39, db40
	, db41, db42, db43, db44, db45
	, coif1, coif2, coif3, coif4, coif5
	, fk4, fk6, fk8, fk10, fk12
	, fk14, fk16, fk18, fk20, fk22
	, sym2, sym3, sym4, sym5
	, sym6, sym7, sym8, sym9, sym10
	, sym11, sym12, sym13, sym14, sym15
	, sym16, sym17, sym18, sym19, sym20
	, sym21, sym22, sym23, sym24, sym25
	, sym26, sym27, sym28, sym29, sym30
	, sym31, sym32, sym33, sym34, sym35
	, sym36, sym37, sym38, sym39, sym40
	, sym41, sym42, sym43, sym44, sym45
	, dmey
	, bior1_1, bior1_3, bior1_5
	, bior2_2, bior2_4, bior2_6, bior2_8
	, bior3_1, bior3_3, bior3_5, bior3_7, bior3_9
	, bior4_4, bior5_5, bior6_8
	, rbio1_1, rbio1_3, rbio1_5
	, rbio2_2, rbio2_4, rbio2_6, rbio2_8
	, rbio3_1, rbio3_3, rbio3_5, rbio3_7
	, rbio3_9, rbio4_4, rbio5_5, rbio6_8
};



struct Wavelet
{
	std::span<const double> LO_D;
	std::span<const double> HI_D;
	std::span<const double> LO_R;
	std::span<const double> HI_R;
	Wavelet() = default;
	Wavelet(std::span<const double> lo_d, std::span<const double> hi_d, std::span<const double> lo_r, std::span<const double> hi_r) :LO_D(lo_d), HI_D(hi_d), LO_R(lo_r), HI_R(hi_r) {}
	constexpr size_t size() const { return LO_D.size(); }
};