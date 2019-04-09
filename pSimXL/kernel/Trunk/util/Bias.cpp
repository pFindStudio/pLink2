#include "../include/sdk.h"
#include "../include/bias.h"

using namespace std;

namespace sdk
{
Bias::Bias() :
		m_lfPeptideTol(10), m_lfFragmentTol(20), haveNextNextGaussian(false), nextNextGaussian(0.0)
{
	//srand(time(NULL));
}

Bias::~Bias()
{

}
//标准正太分布N(0,1)：https://www.cnblogs.com/tsingke/p/6194737.html
double Bias:: nextGaussian()
{
	if (haveNextNextGaussian) {
		haveNextNextGaussian = false;
		return nextNextGaussian;
	} else {
		double v1, v2, s;
		do {
			double u1 = (double)rand() / RAND_MAX;
			double u2 = (double)rand() / RAND_MAX;
			v1 = 2 * u1 - 1;
			v2 = 2 * u2 - 1;
			s = v1 * v1 + v2 * v2;
		} while (s >= 1 || s == 0);
		double multiplier = sqrt(-2 * log(s)/s);
		nextNextGaussian = v2 * multiplier;
		haveNextNextGaussian = true;
		return v1 * multiplier;
	}
}

double Bias:: nextGaussian(double a, double b)
{
	if (haveNextNextGaussian) {
		haveNextNextGaussian = false;
		return nextNextGaussian;
	} else {
		double v1, v2, s;
		do {
			double u1 = (double)rand() / RAND_MAX;
			double u2 = (double)rand() / RAND_MAX;
			v1 = 2 * u1 - 1;
			v2 = 2 * u2 - 1;
			s = v1 * v1 + v2 * v2;
		} while (s >= 1 || s == 0);
		double multiplier = sqrt(-2 * log(s)/s);
		nextNextGaussian = v2 * multiplier;
		nextNextGaussian = nextNextGaussian * b + a;
		haveNextNextGaussian = true;
		return v1 * multiplier * b + a; //产生均值为a,方差为b的随机信号
	}
}

double Bias::genFragmentBias()
{
	double lfBias = 3 * nextGaussian();
	return lfBias * 0.000001;
	//return m_lfFragmentTol * nextGaussian();

}

double Bias::genPeptideBias()
{
	//double lfBias = m_lfPeptideTol * nextGaussian();
	double lfBias = 3 * nextGaussian(); //乘以3，变成N(0,3^2)
	if(lfBias > 9.9)
		lfBias = 9.9;
	return lfBias * 0.000001;
}

double Bias::genIntensityBias()
{
	double lfBias = nextGaussian();
	return lfBias;
}

double Bias::generateNoise()
{
	double lfBias = nextGaussian();
	return lfBias;
}

}
