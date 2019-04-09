#include "../include/scorers.h"

using namespace std;

namespace sdk
{

BMIonSimpleGenerator::BMIonSimpleGenerator(const SearchParameter *pParameter) :
		RefinedScorer(pParameter), m_pBias(0)
{
	m_vBasicAA.push_back('K');
	m_vBasicAA.push_back('R');
	m_vBasicAA.push_back('H');
	srand(time(NULL));
	if(m_pBias == NULL) {
		m_pBias = new Bias;
	}
	genMatchRatio();
	genIntensityMean();
}

const string BMIonSimpleGenerator::strTitlePrefix("Simulation");
const string BMIonSimpleGenerator::strTitleSuffix("dta");
size_t BMIonSimpleGenerator::m_tScanNo = 1; // scan从1开始编号

BMIonSimpleGenerator::~BMIonSimpleGenerator()
{
	if(m_pBias)
		delete m_pBias;
	m_pBias = NULL;
}

size_t BMIonSimpleGenerator::generateCharge()
{
	double arrPro[] = {0, 0, 0, 1257, 801, 165, 66};
	vector<double> vPro(arrPro, arrPro + sizeof(arrPro)/sizeof(double) );
	for(size_t i = 1; i < vPro.size(); i++)
	{
		vPro[i] += vPro[i-1];
	}
	for(size_t i = 1; i < vPro.size(); i++)
	{
		vPro[i] /= vPro[vPro.size()-1];
	}
	double lfRandom = (double)rand() / RAND_MAX;
	for(size_t i = 1; i < vPro.size(); i++)
	{
		if(lfRandom < vPro[i])
		{
			return i;
		}
	}
	return vPro.size();
}

void BMIonSimpleGenerator::insertRatio(map<pair<char, int>, double> &mpStatistic, char cType, int nCharge, double lfRatio) {
	mpStatistic.insert(pair<pair<char, int>, double>(pair<char, int>(cType, nCharge), lfRatio));
}

void BMIonSimpleGenerator::genMatchRatio()
{
	insertRatio(m_mpIonRatio, 'b', 1, 0.5);
	insertRatio(m_mpIonRatio, 'b', 2, 0.5);
	insertRatio(m_mpIonRatio, 'y', 1, 0.8);
	insertRatio(m_mpIonRatio, 'y', 2, 0.5);

}

void BMIonSimpleGenerator::genIntensityMean()
{
	insertRatio(m_mpIonIntensity, 'b', 1, 0.11);
	insertRatio(m_mpIonIntensity, 'b', 2, 0.07);
	insertRatio(m_mpIonIntensity, 'y', 1, 0.21);
	insertRatio(m_mpIonIntensity, 'y', 2, 0.12);
}

double BMIonSimpleGenerator::getMatchRatio(IonType stIonType)
{
	double lfMatchRatio = 0.0;
	if(m_mpIonRatio.count(pair<char, int>(stIonType.cType, stIonType.nCharge)) > 0)
		lfMatchRatio = m_mpIonRatio.at(pair<char, int>(stIonType.cType, stIonType.nCharge));
	return lfMatchRatio;
}

double BMIonSimpleGenerator::getIntensityMean(IonType stIonType)
{
	double lfRatio = 0.0;
	if(m_mpIonIntensity.count(pair<char, int>(stIonType.cType, stIonType.nCharge)) > 0)
		lfRatio = m_mpIonIntensity.at(pair<char, int>(stIonType.cType, stIonType.nCharge));
	return lfRatio;
}

void BMIonSimpleGenerator::mergePeaks()
{
	sort(m_vPeaks.begin(), m_vPeaks.end(), Spectrum::mzLesser);
	for(size_t i = 0; i < m_vPeaks.size(); i++) {
		for(size_t j = i+1; j < m_vPeaks.size(); j++) {
			if(m_vPeaks[j].m_lfMz - m_vPeaks[i].m_lfMz >= 0.00001) {
				m_pSpec->m_vPeaks.push_back(m_vPeaks[j-1]);
				break;
			}
			else {
				m_vPeaks[j].m_lfIntensity += m_vPeaks[i].m_lfIntensity;
				i = j;
			}
		}
	}
	if(0 != m_vPeaks.size()) {
		m_pSpec->m_vPeaks.push_back(m_vPeaks[m_vPeaks.size()-1]);
	}
}

double BMIonSimpleGenerator::score()
{
	m_pSpec->m_tCharge = generateCharge();
	if(!m_bComputeMz) {
		computeMz();
		m_bComputeMz = true;
	}
	m_pSpec->m_tScanNo = m_tScanNo++;
	m_pSpec->m_lfMZ = (m_pPeptide->m_lfCalcMH + (m_pSpec->m_tCharge-1) * PROTON_MASS) / m_pSpec->m_tCharge;
	double lfError = 0.0;
	lfError = m_pSpec->m_lfMZ * m_pBias->genPeptideBias();
	m_pSpec->m_lfMZ += lfError;
	ostringstream oss;
	oss<<strTitlePrefix<<"."
		<<m_pSpec->m_tScanNo<<"."
		<<m_pSpec->m_tScanNo<<"."
		<<m_pSpec->m_tCharge<<"."
		<<strTitleSuffix;
	m_pSpec->m_strTitle = oss.str();
	match();
	vector<IonPeak>().swap(m_pSpec->m_vPeaks);
	mergePeaks();

	return 0.0;
}

void BMIonSimpleGenerator::match()
{
	vector<IonPeak>().swap(m_vPeaks);
	InstrumentDict *pDict = InstrumentDict::getInstance();
	const Instrument &stInst = pDict->getInstrument(m_pParameter->m_strRefinedInstrument);
	const vector<IonType> &vIonTypes = stInst.vIonTypes;
	std::vector<std::vector<bool> > bContinuity;
	bContinuity.resize(vIonTypes.size(), vector<bool>(MAX_PEPTIDE_LEN * 2 + 2, false));
	int AlphaPepLen = m_pPeptide->m_stAlphaPep.m_stPep.m_strSq.length();
	for(size_t tPepPos = 0; tPepPos < m_vIonMz.size(); ++tPepPos) {
		if((m_pSpec->m_tCharge == 1 && (size_t)vIonTypes[m_vIonMz[tPepPos].yIonTypeOrder].nCharge > m_pSpec->m_tCharge)
				|| (m_pSpec->m_tCharge > 1 && (size_t)vIonTypes[m_vIonMz[tPepPos].yIonTypeOrder].nCharge >= m_pSpec->m_tCharge))
			continue;

		int nOrder = m_vIonMz[tPepPos].yIonTypeOrder;

		int nRandom = rand() % 100;
		if(nRandom > getMatchRatio(vIonTypes[nOrder]) * 100) continue; // 以不同概率缺峰

		IonPeak stIonPeak;
		stIonPeak.m_lfIntensity = INTENSITY_MULTIPLIER * getIntensityMean(vIonTypes[nOrder]);
		stIonPeak.m_lfMz = m_vIonMz[tPepPos].lfMz;
		stIonPeak.m_lfMz += stIonPeak.m_lfMz * m_pBias->genFragmentBias();
		m_vPeaks.push_back(stIonPeak);

		double lfMonoNeuMass = (m_vIonMz[tPepPos].lfMz - PROTON_MASS) * vIonTypes[nOrder].nCharge;
		double lfMonoIntensity = stIonPeak.m_lfIntensity;

		//生成第一同位素峰
		stIonPeak.m_lfMz = m_vIonMz[tPepPos].lfMz + 1.0 * PROTON_MASS/vIonTypes[nOrder].nCharge;
		stIonPeak.m_lfMz += stIonPeak.m_lfMz * m_pBias->genFragmentBias();
		stIonPeak.m_lfIntensity = INTENSITY_MULTIPLIER * getIntensityMean(vIonTypes[nOrder]);
		stIonPeak.m_lfIntensity = lfMonoIntensity * lfMonoNeuMass * 5.43 * 1e-4; // Table S-1 of Park et al., AC 2008
		m_vPeaks.push_back(stIonPeak);

		double lfI1Intensity = stIonPeak.m_lfIntensity;

		//生成第二同位素峰
		stIonPeak.m_lfMz = m_vIonMz[tPepPos].lfMz + 2.0 * PROTON_MASS/vIonTypes[nOrder].nCharge;
		stIonPeak.m_lfMz += stIonPeak.m_lfMz * m_pBias->genFragmentBias();
		stIonPeak.m_lfIntensity = INTENSITY_MULTIPLIER * getIntensityMean(vIonTypes[nOrder]);
		stIonPeak.m_lfIntensity = lfI1Intensity * (lfMonoNeuMass * 2.71 * 1e-4 + 8.17 * 1e-2);
		m_vPeaks.push_back(stIonPeak);

		if(rand()%5 == 0){ // 20%概率生成噪声峰
			stIonPeak.m_lfMz = m_vIonMz[tPepPos].lfMz;
			stIonPeak.m_lfMz += 6.5 * m_pBias->generateNoise(); // (0, 6.5^2) ppm 3倍方差为20ppm
			stIonPeak.m_lfIntensity = INTENSITY_MULTIPLIER * 0.01;
			m_vPeaks.push_back(stIonPeak);
		}
	}
}

} /* namespace sdk */



