#include "../include/decorators.h"
#include "../include/index.h"

using namespace std;

namespace sdk
{

EnumPeptideIndex::EnumPeptideIndex(const SearchParameter *pParameter) :
		m_pParameter(pParameter), m_pModDecorator(NULL), m_pXLinkDecorator(NULL),
		m_tMassRange(size_t(
				(pParameter->m_lfMaxPepMass - pParameter->m_lfMinPepMass) *
				MZ_INDEX_SCALE) ), m_vCnt(m_tMassRange+1, 0)
{
	try {
		m_pModDecorator = new ModificationDecorator<EnumPeptideIndex, EnumPeptideIndexDispatch>(m_pParameter);
		m_pXLinkDecorator = new XLinkDecorator(m_pParameter);
	} catch(exception &e) {
		ErrorInfo err("EnumPeptideIndex", "EnumPeptideIndex", "caught an exception", e);
		throw runtime_error(err.get());
	} catch(...) {
		ErrorInfo err("EnumPeptideIndex", "EnumPeptideIndex", "caught an unknown exception");
		throw runtime_error(err.get());
	}
}

EnumPeptideIndex::~EnumPeptideIndex()
{
	if(m_pModDecorator) {
		delete m_pModDecorator;
		m_pModDecorator = NULL;
	}

	if(m_pXLinkDecorator) {
		delete m_pXLinkDecorator;
		m_pXLinkDecorator = NULL;
	}
}

void EnumPeptideIndex::clear()
{
	fill(m_vCnt.begin(), m_vCnt.end(), 0);
	vector<Peptide>().swap(m_vPepTable);
	m_vPepTable.clear();
}

void EnumPeptideIndex::getPeptides(std::vector<Peptide*> &vPeps)
{
	for(size_t i = 0; i < m_vPepTable.size(); i++)
	{
		vPeps.push_back(&m_vPepTable[i]);
	}
}

void EnumPeptideIndex::getXLinkPeptides(std::vector<Peptide*> &vPeps)
{
	for(size_t i = 0; i < m_vPepTable.size(); i++)
	{
		if(m_pXLinkDecorator->isXLinked(m_vPepTable[i]))
		{
			vPeps.push_back(&m_vPepTable[i]);
		}
	}
}

void EnumPeptideIndex::attachPeptides(size_t tStartId, size_t tEndId, const ProteinIndex *pProIndex)
{
	try {
		EnzymeDict *pEnzymeDict = EnzymeDict::getInstance();
		EnzymeType eType = pEnzymeDict->getEnzymeType(m_pParameter->m_strEnzymeName);
		std::vector<Peptide> vPepTable;
		for(size_t tId = tStartId; tId < tEndId; ++tId) {
			switch(eType) {
			case ET_SPECIFIC:
				attachSpecific(tId, pProIndex, vPepTable);
				break;
			default:
				throw runtime_error("unknown cleave way.\n");
			}
		}

		getUsefulPep(vPepTable);

	} catch(exception &e) {
		ErrorInfo err("EnumPeptideIndex", "attachPeptides", "caught an exception.", e);
		throw runtime_error(err.get());
	} catch(...) {
		ErrorInfo err("EnumPeptideIndex", "attachPeptides", "caught an unknown exception.");
		throw runtime_error(err.get());
	}
}

void EnumPeptideIndex::attachSpecific(size_t tProteinId, const ProteinIndex *pProIndex, std::vector<Peptide> &vPepTable)
{
	size_t tPepMinLen = m_pParameter->m_tMinPepLen;
	size_t tPepMaxLen = m_pParameter->m_tMaxPepLen;
	double lfPepMinMass = m_pParameter->m_lfMinPepMass;
	double lfPepMaxMass = m_pParameter->m_lfMaxPepMass;
	size_t tMaxMissSite = m_pParameter->m_tMaxMissSite;

	size_t tProStart = pProIndex->m_vPrtnEntr[tProteinId];
	size_t tProEnd = pProIndex->m_vPrtnEntr[tProteinId+1];
	size_t tSiteIdxStart = pProIndex->m_vSiteEntr[tProteinId];
	size_t tSiteIdxEnd = pProIndex->m_vSiteEntr[tProteinId+1];

	size_t tPepStart(tProStart);
	size_t tPepEnd(0);

	PeptideEntry stPepEntry;
	stPepEntry.tProStart = tProStart;
	stPepEntry.tProLen = tProEnd - tProStart;

	size_t tCurSiteIdx(tSiteIdxStart);
	size_t tMissSite(0);

	while(true) {
		if(tCurSiteIdx < tSiteIdxEnd) {
			tPepEnd = pProIndex->m_vSite[tCurSiteIdx];
		} else if(tCurSiteIdx == tSiteIdxEnd) {
			if(tSiteIdxStart < tSiteIdxEnd && pProIndex->m_vSite[tSiteIdxEnd-1] == tProEnd) {
				if(tMissSite <= 1) {
                    break;
                } else {
                    if(tPepStart == tProStart && pProIndex->m_strSq[tPepStart] == 'M') {
				        ++tPepStart;
				        tCurSiteIdx = tCurSiteIdx - tMissSite;
			        } else if(tCurSiteIdx - tMissSite < tSiteIdxEnd) {
				        tPepStart = pProIndex->m_vSite[tCurSiteIdx-tMissSite];
				        tCurSiteIdx = tCurSiteIdx + 1 - tMissSite;
			        } else {
				        tPepStart = tProEnd;
				        tCurSiteIdx = tCurSiteIdx + 1;
			        }
                    tMissSite = 0;
                }
            }
			tPepEnd = tProEnd;
		} else {
			break;
		}

		size_t tPepLen = tPepEnd - tPepStart;
		double lfMass = pProIndex->m_vSqMass[tPepEnd] - pProIndex->m_vSqMass[tPepStart];
		if(tPepLen >= tPepMinLen && tPepLen <= tPepMaxLen &&
				lfMass >= lfPepMinMass && lfMass <= lfPepMaxMass) {
			stPepEntry.tPepStart = tPepStart;
			stPepEntry.tPepLen = tPepLen;

			m_stCurPep.clear();
			string strSq = pProIndex->m_strSq.substr(tPepStart, tPepLen);
			m_stCurPep.setPeptideInfo(strSq, lfMass, tMissSite, getPeptideType(stPepEntry));

			m_stCurPep.m_lfTag = m_stCurPep.getGodelCode();
			vPepTable.push_back(m_stCurPep);
		}

		if(tCurSiteIdx < tSiteIdxEnd && tMissSite < tMaxMissSite) {
			++tMissSite;
			++tCurSiteIdx;
		} else {
			if(tPepStart == tProStart && pProIndex->m_strSq[tPepStart] == 'M') {
				++tPepStart;
				tCurSiteIdx = tCurSiteIdx - tMissSite;
			} else if(tCurSiteIdx - tMissSite < tSiteIdxEnd) {
				tPepStart = pProIndex->m_vSite[tCurSiteIdx-tMissSite];
				tCurSiteIdx = tCurSiteIdx + 1 - tMissSite;
			} else {
				tPepStart = tProEnd;
				tCurSiteIdx = tCurSiteIdx + 1;
			}
			tMissSite = 0;
		}
	}
}


void EnumPeptideIndex::attachPTMPeptides()
{
	if(m_stCurPep.m_lfMass >= m_pParameter->m_lfMinPepMass &&
			m_stCurPep.m_lfMass <= m_pParameter->m_lfMaxPepMass) {
		m_vPepTable.push_back(m_stCurPep);
	}
}

void EnumPeptideIndex::getUsefulPep(vector<Peptide> &vPepTable)
{
	vector<Peptide *> vTempPepTable;
	for (size_t i = 0; i < vPepTable.size(); ++i){
		vTempPepTable.push_back(&vPepTable[i]);
	}
	std::sort(vTempPepTable.begin(), vTempPepTable.end(), Peptide::massTagLesser);
	if(vTempPepTable.size() > 0)
	{
		m_stCurPep = *vTempPepTable[0];
		m_pModDecorator->addModifications(m_stCurPep, this, &EnumPeptideIndex::attachPTMPeptides);
	}
	for (size_t i = 1; i < vTempPepTable.size(); ++i){
		if(fabs(vTempPepTable[i]->m_lfMass - vTempPepTable[i-1]->m_lfMass) < PRECISION &&
				fabs(vTempPepTable[i]->m_lfTag - vTempPepTable[i-1]->m_lfTag) < PRECISION)
		{
			continue;
		}
		m_stCurPep = *vTempPepTable[i];
		m_pModDecorator->addModifications(m_stCurPep, this, &EnumPeptideIndex::attachPTMPeptides);
	}
}

unsigned char EnumPeptideIndex::getPeptideType(const PeptideEntry &stEntry)
{
	if(stEntry.tPepStart == stEntry.tProStart &&
			stEntry.tPepStart + stEntry.tPepLen == stEntry.tProStart + stEntry.tProLen)
		return 3;
	else if(stEntry.tPepStart == stEntry.tProStart) {
		return 1;
	} else if(stEntry.tPepStart + stEntry.tPepLen == stEntry.tProStart + stEntry.tProLen){
		return 2;
	} else {
		return 0;
	}
}

}

