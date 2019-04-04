
#include "../include/scorers.h"

using namespace std;

namespace sdk
{


RefinedScorer::RefinedScorer(const SearchParameter *pParameter) :
		m_pParameter(pParameter), m_pSpec(NULL), m_pPeptide(NULL),
		m_pCalculator(NULL), m_bComputeMz(false)
{
	m_vIonMz.reserve(10 * MAX_IONTYPE_NUM * 2 * MAX_PEPTIDE_LEN);
}

RefinedScorer::~RefinedScorer()
{
}

void RefinedScorer::setSpectrum(Spectrum *pSpec)
{
	m_pSpec = pSpec;
}

void RefinedScorer::setPeptide(XLinkPeptideResult *pPeptide)
{
	m_pPeptide = pPeptide;
	m_bComputeMz = false;
}

void RefinedScorer::computeMz()
{
	if(m_pCalculator == NULL) {
		m_pCalculator = new MzCalculator();
		m_pCalculator->init(m_pPeptide);

		m_vIonMz.clear();
		InstrumentDict *pDict = InstrumentDict::getInstance();
		const Instrument &stInst = pDict->getInstrument(m_pParameter->m_strRefinedInstrument);

		for(size_t i = 0; i < stInst.vIonTypes.size(); ++i) {
			switch(stInst.vIonTypes[i].cType) {
			case 'a': case 'b': case 'c':
				m_pCalculator->attachNTermIons(m_vIonMz, stInst.vIonTypes[i], i);
				break;
			case 'x': case 'y': case 'z':
				m_pCalculator->attachCTermIons(m_vIonMz, stInst.vIonTypes[i], i);
				break;
			default:
				break;
			}
		}

#ifdef DEBUG_SCORE
//		Test for thereotical ions
		Trace::getInstance()->info("%s", m_pPeptide->toString().c_str());
		Trace::getInstance()->info("ion number: %d", m_vIonMz.size());
		for(size_t tIdx = 0; tIdx < m_vIonMz.size(); ++tIdx) {
			Trace::getInstance()->info("%c: %s", stInst.vIonTypes[(int)m_vIonMz[tIdx].yIonTypeOrder].cType, m_vIonMz[tIdx].toString().c_str());
		}
#endif
	}

	if(m_pCalculator) {
		delete m_pCalculator;
		m_pCalculator = NULL;
	}
}

void RefinedScorer::match()
{
	InstrumentDict *pDict = InstrumentDict::getInstance();
	const Instrument &stInst = pDict->getInstrument(m_pParameter->m_strRefinedInstrument);
	const vector<IonType> &vIonTypes = stInst.vIonTypes;

	m_vMatchedInfo.clear();
	m_vMatchedInfo.resize(m_pSpec->m_vPeaks.size(), vector<int>());

	for(size_t tPepPos = 0; tPepPos < m_vIonMz.size(); ++tPepPos) {
		if((m_pSpec->m_tCharge == 1 && (size_t)vIonTypes[m_vIonMz[tPepPos].yIonTypeOrder].nCharge > m_pSpec->m_tCharge)
				|| (m_pSpec->m_tCharge > 1 && (size_t)vIonTypes[m_vIonMz[tPepPos].yIonTypeOrder].nCharge >= m_pSpec->m_tCharge))
			continue;

		int nMin, nMax;
		getMassBorder(nMin, nMax, m_vIonMz[tPepPos].nMz);

		int nTemp = nMin / MZ_MULTIPLIER;

		if(nTemp >= (int)m_pSpec->m_vHashIndex.size() || nTemp < 0)
			continue;

		size_t tSpecPos = m_pSpec->m_vHashIndex[nTemp];
		while(tSpecPos < m_pSpec->m_vPeaks.size() && m_pSpec->m_vPeaks[tSpecPos].m_nMz < nMin)
			++tSpecPos;

		while(tSpecPos < m_pSpec->m_vPeaks.size() && m_pSpec->m_vPeaks[tSpecPos].m_nMz <= nMax) {
			if(m_pSpec->m_vPeaks[tSpecPos].m_tCharge == 0 ||
					m_pSpec->m_vPeaks[tSpecPos].m_tCharge == (size_t)vIonTypes[m_vIonMz[tPepPos].yIonTypeOrder].nCharge) {
				m_vMatchedInfo[tSpecPos].push_back(tPepPos);

				// Test for matched peaks
//				Trace::getInstance()->info("matched peaks: %.4lf, %c, %s",
//						m_pSpec->m_vPeaks[tSpecPos].m_lfMz,
//						vIonTypes[(int)m_vIonMz[tPepPos].yIonTypeOrder].cType,
//						m_vIonMz[tPepPos].toString().c_str());
			}
			++tSpecPos;
		}
	}
}



void RefinedScorer::getMassBorder(int &nMin, int &nMax, int nMz)
{
	double lfTol = m_pParameter->m_lfFragTol;
	if(m_pParameter->m_eFragTolType == TT_PPM)
		lfTol = m_pParameter->m_lfFragTol * nMz * PART_PER_MILLION;

	nMin = nMz - lfTol;
	nMax = nMz + lfTol;
}

double RefinedScorer::getError(size_t tThrMz, size_t tExpMz)
{
	double lfError = (tExpMz >= tThrMz) ? (tExpMz - tThrMz) : (-1.0 * (tThrMz - tExpMz));
	return MILLION * lfError / tThrMz;
}

RefinedScorerFactory::RefinedScorerFactory(const SearchParameter *pParameter) :
		m_pParameter(pParameter)
{
}

RefinedScorer *RefinedScorerFactory::getScorer(ScoreType eScoreType) const
{
	switch(eScoreType) {
	case ST_XLINK_LEGACY:
		return new BMIonGenerator(m_pParameter);
	case ST_XLINK_PERFECT:
		return new BMIonPerfectGenerator(m_pParameter);
	case ST_XLINK_SIMPLE:
		return new BMIonSimpleGenerator(m_pParameter);
	default:
        ErrorInfo err("RefinedScorerFactory", "getScorer", "unkown RefinedScorer type.");
        throw runtime_error(err.get());
        return NULL;
	}
}

} // end of namespace




