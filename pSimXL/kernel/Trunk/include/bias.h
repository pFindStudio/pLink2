
#ifndef BIAS_H_
#define BIAS_H_
#include "util.h"
#include "sdk.h"
namespace sdk
{
class Bias {
private:
	double m_lfPeptideTol;
	double m_lfFragmentTol;
	bool haveNextNextGaussian;
	double nextNextGaussian;
public:
	Bias();
	virtual ~Bias();
	double  nextGaussian();
	double  nextGaussian(double a, double b);
	virtual double genPeptideBias();
	virtual double genFragmentBias();
	virtual double genIntensityBias();
	virtual double generateNoise();
};

} /* namespace sdk */

#endif /* BIAS_H_ */
