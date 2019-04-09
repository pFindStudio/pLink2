#include "../include/scorers.h"

using namespace std;

namespace sdk
{

BMIonPerfectGenerator::BMIonPerfectGenerator(const SearchParameter *pParameter) :
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
	genExistFlag();
}

const string BMIonPerfectGenerator::strTitlePrefix("Simulation");
const string BMIonPerfectGenerator::strTitleSuffix("dta");
//size_t BMIonGenerator::m_tScanNo = 1000;
size_t BMIonPerfectGenerator::m_tScanNo = 1; // scan从1开始编号

BMIonPerfectGenerator::~BMIonPerfectGenerator()
{
	if(m_pBias)
		delete m_pBias;
	m_pBias = NULL;
}

size_t BMIonPerfectGenerator::generateCharge()
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

void BMIonPerfectGenerator::insertRatio(map<pair<char, int>, double> &mpStatistic, char cType, int nCharge, double lfRatio) {
	mpStatistic.insert(pair<pair<char, int>, double>(pair<char, int>(cType, nCharge), lfRatio));
}

void BMIonPerfectGenerator::insertRatio(map<pair<char, int>, bool> &mpStatistic, char cType, int nCharge, double lfRatio) {
	double lfTotal = 1141.0;
	double lfRandom = (double)rand() / RAND_MAX;
	mpStatistic.insert(pair<pair<char, int>, bool>(pair<char, int>(cType, nCharge), lfRandom < lfRatio/lfTotal));
}

void BMIonPerfectGenerator::genMatchRatio()
{

	if(m_pParameter->m_vLinkers[0] == "SS") { // SS ion gain ratio均值

		Trace::getInstance()->alert("Use SS ion gain ratio.");

		insertRatio(m_mpIonRatio, 'y', 1, 0.49);
		insertRatio(m_mpIonRatio, 'b', 1, 0.3);
		insertRatio(m_mpIonRatio, 'y', 2, 0.3);
		insertRatio(m_mpIonRatio, 'a', 1, 0.19);
		insertRatio(m_mpIonRatio, 'b', 2, 0.25);
		insertRatio(m_mpIonRatio, 'p', 1, 0.61);
		insertRatio(m_mpIonRatio, 'p', 2, 0.59);
		insertRatio(m_mpIonRatio, 'p', 3, 0.55);
		insertRatio(m_mpIonRatio, 'y', 3, 0.17);
		insertRatio(m_mpIonRatio, 'a', 2, 0.17);

	} else { // BS3 ion gain ratio均值

		Trace::getInstance()->alert("Use BS3 ion gain ratio.");

		insertRatio(m_mpIonRatio, 'b', 1, 0.19);
		insertRatio(m_mpIonRatio, 'b', 2, 0.23);
		insertRatio(m_mpIonRatio, 'y', 1, 0.71);
		insertRatio(m_mpIonRatio, 'y', 2, 0.37);
		insertRatio(m_mpIonRatio, 'y', 3, 0.08);
		insertRatio(m_mpIonRatio, 'a', 1, 0.08);
		insertRatio(m_mpIonRatio, 'a', 2, 0.09);

	}

}

void BMIonPerfectGenerator::genIntensityMean()
{

	if(m_pParameter->m_vLinkers[0] == "SS") { // SS谱峰强度均值

		Trace::getInstance()->alert("Use SS average intensity.");

		insertRatio(m_mpIonIntensity, 'y', 1, 0.20);
		insertRatio(m_mpIonIntensity, 'b', 1, 0.11);
		insertRatio(m_mpIonIntensity, 'y', 2, 0.09);
		insertRatio(m_mpIonIntensity, 'a', 1, 0.20);
		insertRatio(m_mpIonIntensity, 'b', 2, 0.06);
		insertRatio(m_mpIonIntensity, 'p', 1, 0.06);
		insertRatio(m_mpIonIntensity, 'p', 2, 0.08);
		insertRatio(m_mpIonIntensity, 'p', 3, 0.07);
		insertRatio(m_mpIonIntensity, 'y', 3, 0.06);
		insertRatio(m_mpIonIntensity, 'a', 2, 0.06);

	} else { // BS3谱峰强度均值

		Trace::getInstance()->alert("Use BS3 average intensity.");

		insertRatio(m_mpIonIntensity, 'b', 1, 0.12);
		insertRatio(m_mpIonIntensity, 'b', 2, 0.08);
		insertRatio(m_mpIonIntensity, 'y', 1, 0.22);
		insertRatio(m_mpIonIntensity, 'y', 2, 0.14);
		insertRatio(m_mpIonIntensity, 'y', 3, 0.05);
		insertRatio(m_mpIonIntensity, 'a', 1, 0.18);
		insertRatio(m_mpIonIntensity, 'a', 2, 0.05);

	}

}

void BMIonPerfectGenerator::genExistFlag()
{
	insertRatio(m_mpIonExsist, 'y', 1, 1130);
	insertRatio(m_mpIonExsist, 'b', 1, 1084);
	insertRatio(m_mpIonExsist, 'y', 2, 1048);
	insertRatio(m_mpIonExsist, 'a', 1, 1057);
	insertRatio(m_mpIonExsist, 'b', 2, 969);
	insertRatio(m_mpIonExsist, 'p', 1, 512);
	insertRatio(m_mpIonExsist, 'p', 2, 450);
	insertRatio(m_mpIonExsist, 'p', 3, 84);
	insertRatio(m_mpIonExsist, 'y', 3, 406);
	insertRatio(m_mpIonExsist, 'a', 2, 910);
}

double BMIonPerfectGenerator::getMatchRatio(IonType stIonType)
{
	double lfMatchRatio = 0.1;
	if(m_mpIonRatio.count(pair<char, int>(stIonType.cType, stIonType.nCharge)) > 0)
		lfMatchRatio = m_mpIonRatio.at(pair<char, int>(stIonType.cType, stIonType.nCharge));
	return lfMatchRatio;
}

double BMIonPerfectGenerator::getIntensityMean(IonType stIonType)
{
	double lfRatio = 0.03;
	if(m_mpIonIntensity.count(pair<char, int>(stIonType.cType, stIonType.nCharge)) > 0)
		lfRatio = m_mpIonIntensity.at(pair<char, int>(stIonType.cType, stIonType.nCharge));
	return lfRatio;
}

bool BMIonPerfectGenerator::getExistFlag(IonType stIonType)
{
	bool bFlag = true;
	if(m_mpIonExsist.count(pair<char, int>(stIonType.cType, stIonType.nCharge)) > 0)
		bFlag = m_mpIonExsist.at(pair<char, int>(stIonType.cType, stIonType.nCharge));
	return bFlag;
}
//double BMIonGenerator::getBasicAANum(XLinkPeptideResult &stXLinkPeptideResult)
//{
//	string strAlpha(stXLinkPeptideResult.m_stAlphaPep.m_stPep.m_strSq);
//	string strBeta=(stXLinkPeptideResult.m_stBetaPep.m_stPep.m_strSq);
//	for(size_t i = 0; i < strAlpha.size(); i++){
//		for(size_t j = 0; j < m_vBasicAA.size(); j++){
//			//if()
//		}
//	}
//}

void BMIonPerfectGenerator::mergePeaks()
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

double BMIonPerfectGenerator::score()
{
	m_pSpec->m_tCharge = generateCharge();
	if(!m_bComputeMz) {
		computeMz();
		m_bComputeMz = true;
	}
	m_pSpec->m_tScanNo = m_tScanNo++;
	m_pSpec->m_lfMZ = (m_pPeptide->m_lfCalcMH + (m_pSpec->m_tCharge-1) * PROTON_MASS) / m_pSpec->m_tCharge;
//	double lfError = 0.0;
//	lfError = m_pSpec->m_lfMZ * m_pBias->genPeptideBias();
//	m_pSpec->m_lfMZ += lfError;
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

double BMIonPerfectGenerator::getPepTypeScore(PeptideType ePepType){
	double arrExtraScore[] = {0.02, 0.02, 0.03, 0.001};
	return arrExtraScore[ePepType];
}

void BMIonPerfectGenerator::match()
{
	vector<IonPeak>().swap(m_vPeaks);
	InstrumentDict *pDict = InstrumentDict::getInstance();
	const Instrument &stInst = pDict->getInstrument(m_pParameter->m_strRefinedInstrument);
	const vector<IonType> &vIonTypes = stInst.vIonTypes;
	std::vector<std::vector<bool> > bContinuity;
	bContinuity.resize(vIonTypes.size(), vector<bool>(MAX_PEPTIDE_LEN * 2 + 2, false));
	int AlphaPepLen = m_pPeptide->m_stAlphaPep.m_stPep.m_strSq.length();
	double lfExtraScore = getPepTypeScore(m_pPeptide->m_ePepType);
	for(size_t tPepPos = 0; tPepPos < m_vIonMz.size(); ++tPepPos) {
		if((m_pSpec->m_tCharge == 1 && (size_t)vIonTypes[m_vIonMz[tPepPos].yIonTypeOrder].nCharge > m_pSpec->m_tCharge)
				|| (m_pSpec->m_tCharge > 1 && (size_t)vIonTypes[m_vIonMz[tPepPos].yIonTypeOrder].nCharge >= m_pSpec->m_tCharge))
			continue;

		int nRandom = rand() % 100;
		if(nRandom < m_pParameter->m_nLostPeakPercentage) continue; // 以m_nLostPeakPercentage%的概率丢峰

		int nOrder = m_vIonMz[tPepPos].yIonTypeOrder;

//		byte yPos = m_vIonMz[tPepPos].yPepPosOrder1;
//		byte yPosSingle = yPos;
//		if(yPosSingle >= AlphaPepLen)
//		{
//			yPosSingle -= AlphaPepLen;
//		}
//		double lfThreshold = 0.0;
//		if(!getExistFlag(vIonTypes[nOrder]))
//			continue;
//		lfThreshold += getMatchRatio(vIonTypes[nOrder]);
//
//		if(vIonTypes[nOrder].bContinuity)
//		{
//			bool bFlag= false;
//			if(yPosSingle > 0)
//			{
//				if(bContinuity[nOrder][yPos-1])
//					bFlag = true;
//			}
//
//			if(bContinuity[nOrder][yPos+1])
//				bFlag = true;
//
//			if(bFlag)
//			{
//				lfThreshold += 0.05;
//			}
//		}
//
//		//最好增加一条碱性氨基酸容易带电荷
//		lfThreshold += lfExtraScore;
//		if(lfThreshold > 0.8)
//			lfThreshold = 0.8;
//		if(lfRandom > lfThreshold)
//			continue;
//
//		if(vIonTypes[nOrder].bContinuity)
//		{
//			bContinuity[nOrder][yPos] = true;
//		}

		IonPeak stIonPeak;
		stIonPeak.m_lfIntensity = INTENSITY_MULTIPLIER * getIntensityMean(vIonTypes[nOrder]);
//		stIonPeak.m_lfIntensity = stIonPeak.m_lfIntensity + 0.4 * stIonPeak.m_lfIntensity * m_pBias->genIntensityBias();
//		if(stIonPeak.m_lfIntensity < 1.0)
//			continue;
		stIonPeak.m_lfMz = m_vIonMz[tPepPos].lfMz;
//		stIonPeak.m_lfMz += stIonPeak.m_lfMz * m_pBias->genFragmentBias();
		m_vPeaks.push_back(stIonPeak);



		//生成同位素模式
//		size_t tIsotopicNun = rand()%4+1;
//		for(size_t i = 0; i < tIsotopicNun; i++) {
//			stIonPeak.m_lfMz += 1.0 * PROTON_MASS/vIonTypes[nOrder].nCharge;
//			stIonPeak.m_lfMz += stIonPeak.m_lfMz * m_pBias->genFragmentBias();
//			stIonPeak.m_lfIntensity = 10000.0 * getIntensityMean(vIonTypes[nOrder]);
//			stIonPeak.m_lfIntensity = (stIonPeak.m_lfIntensity + 0.4 * stIonPeak.m_lfIntensity * m_pBias->genIntensityBias())/2.0;
//			if(stIonPeak.m_lfIntensity < 1.0)
//				continue;
//			m_vPeaks.push_back(stIonPeak);
//		}
//
//		if(rand()%4 == 0){
//			stIonPeak.m_lfMz = m_vIonMz[tPepPos].lfMz;
//			stIonPeak.m_lfMz += ((double)rand() / RAND_MAX * (rand() % 20) + 1) * m_pBias->generateNoise();
//			stIonPeak.m_lfIntensity = 10000.0 * 0.03;
//			stIonPeak.m_lfIntensity = stIonPeak.m_lfIntensity + 0.1 * stIonPeak.m_lfIntensity * m_pBias->genIntensityBias();
//			if(stIonPeak.m_lfMz * (size_t)vIonTypes[m_vIonMz[tPepPos].yIonTypeOrder].nCharge < 50 ||
//					stIonPeak.m_lfIntensity < 1.0) {
//				continue;
//			}
//			else {
//				m_vPeaks.push_back(stIonPeak);
//			}
//
//		}

	}
}

} /* namespace sdk */



