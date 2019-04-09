
#include "../include/decorators.h"

using namespace std;

namespace sdk
{

XLinkDecorator::XLinkDecorator(const SearchParameter *pParameter) :
		m_pParameter(pParameter), m_lfLeastNegativeLinker(0.0)
{
	initLinkers();
	setOnLinkers();
}

bool XLinkDecorator::isXLinked(const Peptide &stPep)
{

    for (size_t j = 0; j < m_pParameter->m_vLinkers.size(); ++j) {
		if (m_bMonoLinked[j][PEPTIDE_N_SITE-'A']) {
			return true;
		}

		if (m_bMonoLinked[j][PEPTIDE_C_SITE-'A']) {
			return true;
		}

		if (m_bMonoLinked[j][PROTEIN_N_SITE-'A'] && stPep.isProNTerm()) {
			return true;
		}

		if (m_bMonoLinked[j][PROTEIN_C_SITE-'A'] && stPep.isProCTerm()) {
			return true;
		}

		for (size_t i = 0; i < stPep.m_strSq.length(); ++i) {
			if(m_bMonoLinked[j][stPep.m_strSq[i]-'A'] &&
					(i+1 < stPep.m_strSq.length() ||
					(i+1 == stPep.m_strSq.length() && stPep.isProCTerm()))) {
				return true;
			}
		}
	}
	
	return false;
}

bool XLinkDecorator::isXLinkSiteLegal(const Peptide &stPep, size_t tIdx, size_t tLinkerId)
{
    if(tIdx == 0) {
		if(stPep.isProNTerm() && m_bMonoLinked[tLinkerId][PROTEIN_N_SITE-'A']) {
			return true;
		} else if(m_bMonoLinked[tLinkerId][PEPTIDE_N_SITE-'A']) {
			return true;
		} else {
			return m_bMonoLinked[tLinkerId][stPep.m_strSq[tIdx]-'A'];
		}
	} else if(tIdx == stPep.m_strSq.length()-1) {
		if(stPep.isProCTerm()
				&& m_bMonoLinked[tLinkerId][PROTEIN_C_SITE-'A']) {
			return true;
		} else if(m_bMonoLinked[tLinkerId][PEPTIDE_C_SITE-'A']) {
			return true;
		} else if(stPep.isProCTerm() && m_bMonoLinked[tLinkerId][stPep.m_strSq[tIdx]-'A']) {
			return true;
		} else {
			return false; // m_bMonoLinked[tLinkerId][stPep.m_strSq[tIdx]-'A'];
		}
	} else {
		return m_bMonoLinked[tLinkerId][stPep.m_strSq[tIdx]-'A'];
	}
}

bool XLinkDecorator::isXLinkSiteLegal(const Peptide &stPep1, size_t tIdx1,
		const Peptide &stPep2, size_t tIdx2, size_t tLinkerId)
{
	size_t tSite1 = getCandXLinkSiteIndex(stPep1, tIdx1, tLinkerId);
	size_t tSite2 = getCandXLinkSiteIndex(stPep2, tIdx2, tLinkerId);
	return m_bLinked[tLinkerId][tSite1][tSite2];
}

double XLinkDecorator::getLeastNegativeLinker() const
{
	return m_lfLeastNegativeLinker;
}

double XLinkDecorator::getLinkerMass(bool bMono, size_t tLinkerId)
{
	XLinkerDict *pDict = XLinkerDict::getInstance();
	const XLinker &stLinker = pDict->getXLinker(m_pParameter->m_vLinkers[tLinkerId]);
	if(bMono) {
		return stLinker.lfMLMonoMassDiff;
	} else {
		return stLinker.lfMonoMassDiff;
	}
}


size_t XLinkDecorator::getCandXLinkSiteIndex(const Peptide &stPep, size_t tIdx, size_t tLinkerId)
{
	if(tIdx == 0) {
		if(stPep.isProNTerm() && m_bMonoLinked[tLinkerId][PROTEIN_N_SITE-'A']) {
			return PROTEIN_N_SITE-'A';
		} else if(m_bMonoLinked[tLinkerId][PEPTIDE_N_SITE-'A']) {
			return PEPTIDE_N_SITE-'A';
		} else {
			return stPep.m_strSq[tIdx]-'A';
		}
	} else if(tIdx == stPep.m_strSq.length()-1) {
		if(stPep.isProCTerm()
				&& m_bMonoLinked[tLinkerId][PROTEIN_C_SITE-'A']) {
			return PROTEIN_C_SITE-'A';
		} else if(m_bMonoLinked[tLinkerId][PEPTIDE_C_SITE-'A']) {
			return PEPTIDE_C_SITE-'A';
		} else {
			return stPep.m_strSq[tIdx]-'A';
		}
	} else {
		return stPep.m_strSq[tIdx]-'A';
	}
}

void XLinkDecorator::initLinkers()
{
	for (size_t k = 0; k < MAX_LINKER_NUM; ++k) {
		for (size_t i = 0; i < LINK_SITE_NUM; ++i) {
			m_bMonoLinked[k][i] = false;
			for (size_t j = 0; j < LINK_SITE_NUM; ++j) {
				m_bLinked[k][i][j] = false;
			}
		}
	}
	m_lfLeastNegativeLinker = 0.0;
}

void XLinkDecorator::setOnLinkers()
{
	XLinkerDict *pDict = XLinkerDict::getInstance();
	const vector<string> &vLinkers = m_pParameter->m_vLinkers;
	for (size_t k = 0; k < vLinkers.size(); ++k) {
		const XLinker &linker = pDict->getXLinker(vLinkers[k]);
		for (size_t i = 0; i < linker.strAlphaAA.length(); ++i) {
			char cSite1 = linker.strAlphaAA[i];

			if (cSite1 == PROTEIN_N_CHAR)
				cSite1 = PROTEIN_N_SITE;
			else if (cSite1 == PROTEIN_C_CHAR)
				cSite1 = PROTEIN_C_SITE;
			else if (cSite1 == PEPTIDE_N_CHAR)
				cSite1 = PEPTIDE_N_SITE;
			else if (cSite1 == PEPTIDE_C_CHAR)
				cSite1 = PEPTIDE_C_SITE;
			else if (cSite1 > 'Z' || cSite1 < 'A') {
				ErrorInfo err("IonIndexFlow", "initLinkers", "unknown linker site symbol " + cSite1);
				throw runtime_error(err.get());
			}

			m_bMonoLinked[k][cSite1 - 'A'] = true;

			for (size_t j = 0; j < linker.strBetaAA.length(); ++j) {
				char cSite2 = linker.strBetaAA[j];

				if (cSite2 == PROTEIN_N_CHAR)
					cSite2 = PROTEIN_N_SITE;
				else if (cSite2 == PROTEIN_C_CHAR)
					cSite2 = PROTEIN_N_SITE;
				else if (cSite2 == PEPTIDE_N_CHAR)
					cSite2 = PEPTIDE_N_SITE;
				else if (cSite2 == PEPTIDE_C_CHAR)
					cSite2 = PEPTIDE_C_SITE;
				else if (cSite2 > 'Z' || cSite2 < 'A') {
					ErrorInfo err("IonIndexFlow", "initLinkers", "unknown linker site symbol " + cSite1);
					throw runtime_error(err.get());
				}

				m_bMonoLinked[k][cSite2 - 'A'] = true;
				m_bLinked[k][cSite1 - 'A'][cSite2 - 'A'] = true;
				m_bLinked[k][cSite2 - 'A'][cSite1 - 'A'] = true;
			}
		}

		// initialize least negative linker
		if(linker.lfMonoMassDiff < m_lfLeastNegativeLinker) {
			m_lfLeastNegativeLinker = linker.lfMonoMassDiff;
		}
		if(linker.lfMLMonoMassDiff < m_lfLeastNegativeLinker) {
			m_lfLeastNegativeLinker = linker.lfMLMonoMassDiff;
		}
	}
}


}

