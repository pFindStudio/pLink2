#include "../include/index.h"

using namespace std;

namespace sdk
{

/* begin: ProteinIndex */
ProteinIndex::ProteinIndex(const SearchParameter *pParameter) :
		m_pParameter(pParameter), m_tProteinNum(0), m_bAttachFinished(false)
{
	 m_vPrtnEntr.reserve(pParameter->m_tProteinsBatchSize + 5);
	 m_strSq.reserve(pParameter->m_tProteinsBatchSize * AVERAGE_SQ_LEN);
	 m_vSqMass.reserve(pParameter->m_tProteinsBatchSize * AVERAGE_SQ_LEN);
	 m_vSiteEntr.reserve(pParameter->m_tProteinsBatchSize + 5);
	 m_vSite.reserve(pParameter->m_tProteinsBatchSize * 50);
}

ProteinIndex::~ProteinIndex()
{
}

void ProteinIndex::attachProtein(const ProteinItem &stItem)
{
	if(m_bAttachFinished)
		return;

	size_t tCurStart = 0;
	if(m_vPrtnEntr.size() == 0) {
		m_vPrtnEntr.push_back(0);
		m_vSqMass.push_back(0.0);
		m_vSiteEntr.push_back(0);
	}
	tCurStart = m_vPrtnEntr.back();

	size_t tCurLen = stItem.strSq.length();
	m_vPrtnEntr.push_back(tCurStart+tCurLen); /* protein index */
	m_strSq += stItem.strSq; /* protein sequences */

	AminoAcidDict *pAADict = AminoAcidDict::getInstance();
	EnzymeDict *pEnzymeDict = EnzymeDict::getInstance();
	const Enzyme &stEnzyme = pEnzymeDict->getEnzyme(m_pParameter->m_strEnzymeName);
	for(size_t tSite = 0; tSite < tCurLen; ++tSite) {
		double lfMass = m_vSqMass.back();
		lfMass += pAADict->getMonoMass(stItem.strSq[tSite]);
		m_vSqMass.push_back(lfMass); /* cumulative mass */

		if((m_vSite.empty() || m_vSite.back() != tCurStart+tSite) &&
				stEnzyme.strNCleaveSite.find(stItem.strSq[tSite]) != string::npos &&
				(tSite == 0 || stEnzyme.strNExceptSite.find(stItem.strSq[tSite-1]) == string::npos)) {
			m_vSite.push_back(tCurStart+tSite);
		}
		if(stEnzyme.strCCleaveSite.find(stItem.strSq[tSite]) != string::npos &&
				(tSite+1 == tCurLen || stEnzyme.strCExceptSite.find(stItem.strSq[tSite+1]) == string::npos)) {
			m_vSite.push_back(tCurStart+tSite+1); /* available cleave sites */
		}
	}
	m_vSiteEntr.push_back(m_vSite.size()); /* sites index */

	++m_tProteinNum;
}

void ProteinIndex::endAttach()
{
	if(m_bAttachFinished)
		return;

	m_bAttachFinished = true;
}

void ProteinIndex::clear()
{
	m_bAttachFinished = false;
	m_tProteinNum = 0;
	m_vPrtnEntr.clear();
	m_vSqMass.clear();
	m_strSq.clear();
	m_vSiteEntr.clear();
	m_vSite.clear();
}

size_t ProteinIndex::getPeptideSequence(std::string &strSq, size_t tStart, size_t tLen)
{
	strSq = m_strSq.substr(tStart, tLen);
	return tLen;
}

size_t ProteinIndex::getPeptideSequence(char *szSequence, size_t tStart, size_t tLen)
{
	size_t t = 0;
	if(szSequence != NULL) {
		for(; t < tLen; ++t)
			*(szSequence + t) = m_strSq.at(tStart+t);
		*(szSequence +t) = '\0';
	} else {
		ErrorInfo err("ProteinIndex", "getPeptideSequence", "caught a null pointer.");
		throw runtime_error(err.get());
	}
	return t;
}
double ProteinIndex::getPeptideMass(size_t tStart, size_t tLen)
{
	return m_vSqMass[tStart+tLen] - m_vSqMass[tStart];
}

size_t ProteinIndex::getMissSiteNum(size_t tStart, size_t tLen)
{
	size_t tNum = 0;

	try {
		size_t tIdx = atProtein(tStart);
		if(tIdx == ALL_BIT_ON) {
			throw runtime_error("Fail to locate the protein.");
		}
		int nLow = m_vSiteEntr[tIdx];
		int nHigh = m_vSiteEntr[tIdx+1]-1;
		int nMid = nLow + (nHigh - nLow) / 2;

		/* binary search for start position of first cleave site in this peptide*/
		while(nLow <= nHigh) {
			nMid = nLow + (nHigh - nLow) / 2;
			if(nMid == nLow && nMid == nHigh) {
				if(tStart >= m_vSite[nMid])
					break;
				else
					return 0;
			} else if(tStart >= m_vSite[nMid] && tStart < m_vSite[nMid+1]) {
				break;
			} else if(tStart < m_vSite[nMid]) {
				nHigh = nMid - 1;
			} else {
				nLow = nMid + 1;
			}
		}

		for(int i = nMid; i <= nHigh && m_vSite[i] <= tStart + tLen; ++i) {
			++tNum;
		}
	} catch(exception &e) {
		ErrorInfo err("IonIndexer", "getMissSiteNum", "caught an exception.", e);
		throw runtime_error(err.get());
	} catch(...) {
		ErrorInfo err("IonIndexer", "getMissSiteNum", "caught an unknown exception.");
		throw runtime_error(err.get());
	}
	return tNum;
}

bool ProteinIndex::isM(size_t tProteinStart) const
{
	if(tProteinStart >= m_strSq.size())
		return false;
	return m_strSq[tProteinStart] == 'M';
}

size_t ProteinIndex::at(size_t tProteinId) const
{
	return m_vPrtnEntr[tProteinId];
}

size_t ProteinIndex::atProtein(size_t tPos) const
{
	int nLow = 0;
	int nHigh = m_tProteinNum - 1;

	while(nLow <= nHigh) {
		int nMid = nLow + (nHigh - nLow) / 2;
		if(tPos >= m_vPrtnEntr[nMid] && tPos < m_vPrtnEntr[nMid+1]) {
			return (size_t)nMid;
		} else if(tPos < m_vPrtnEntr[nMid]) {
			nHigh = nMid - 1;
		} else {
			nLow = nMid + 1;
		}
	}
	ostringstream oss;
	oss<<"index out of range "<<tPos<<".";
	ErrorInfo error("ProteinIndex", "atProtein", oss.str());
	throw runtime_error(error.get());
	return ALL_BIT_ON;
}

size_t ProteinIndex::findEndSite(size_t tSiteIdxStart, size_t tSiteIdxEnd, size_t tPos) const
{
	if(tSiteIdxStart >= tSiteIdxEnd)
		return ALL_BIT_ON;

	int nLow = 0;
	int nHigh = int(tSiteIdxEnd - tSiteIdxStart) - 1;
	while(nLow <= nHigh) {
		int nMid = nLow + (nHigh - nLow) / 2;
		if(m_vSite[nMid+tSiteIdxStart] > tPos && (nMid == 0 || m_vSite[nMid-1+tSiteIdxStart] <= tPos)) {
			return m_vSite[nMid+tSiteIdxStart];
		} else if(m_vSite[nMid+tSiteIdxStart] > tPos && m_vSite[nMid-1+tSiteIdxStart] > tPos) {
			nHigh = nMid - 1;
		} else {
			nLow = nMid + 1;
		}
	}
	return ALL_BIT_ON;
}

bool ProteinIndex::isCleaveSite(size_t tSiteIdxStart, size_t tSiteIdxEnd, size_t tPos) const
{
	if(tSiteIdxStart >= tSiteIdxEnd)
		return false;

	int nLow = 0;
	int nHigh = tSiteIdxEnd - tSiteIdxStart - 1;
	while(nLow <= nHigh) {
		int nMid = nLow + (nHigh - nLow) / 2;
		if(tPos == m_vSite[nMid+tSiteIdxStart]) {
			return true;
		} else if(tPos < m_vSite[nMid+tSiteIdxStart]) {
			nHigh = nMid - 1;
		} else {
			nLow = nMid + 1;
		}
	}
	return false;
}

/* begin: for debug */
void ProteinIndex::showIndex()
{
	if(!m_bAttachFinished)
		return;

	printf("begin to show the protein index information......\n");
	printf("This load contains "FORMAT_SIZE_T" proteins.\n", m_tProteinNum);
	for(size_t i = 0; i < m_tProteinNum; i++) {
		printf("protein "FORMAT_SIZE_T" (start, length): ("FORMAT_SIZE_T", "FORMAT_SIZE_T")\n", i,
				m_vPrtnEntr[i], m_vPrtnEntr[i+1] - m_vPrtnEntr[i]);
		printf("sequence: \n%s\n",
				m_strSq.substr(m_vPrtnEntr[i],
						m_vPrtnEntr[i+1]-m_vPrtnEntr[i]).c_str());
		printf("cumulative mass: \n");
		for(size_t j = m_vPrtnEntr[i]; j < m_vPrtnEntr[i+1]; ++j) {
			printf("%lf\t", m_vSqMass[j+1]-m_vSqMass[m_vPrtnEntr[i]]);
		}
		printf("\n");

		printf("enzyme site map:\n");
		printf("protein "FORMAT_SIZE_T" (start, length): ("FORMAT_SIZE_T", "FORMAT_SIZE_T")\n", i,
				m_vSiteEntr[i],
				m_vSiteEntr[i+1]-m_vSiteEntr[i]);
		printf("\t sites: \n");
		for(size_t tSite = m_vSiteEntr[i];
				tSite < m_vSiteEntr[i+1]; ++tSite) {
			printf("\t "FORMAT_SIZE_T"\t", m_vSite[tSite]);
		}
		printf("\n\n");
	}
	printf("end................................................\n\n");
}
/* end: for debug */
/* end: ProteinIndex */

}
