
#include "../include/index.h"

using namespace std;

namespace sdk {

IonIndexer::IonIndexer(const SearchParameter *pParameter) :
		m_pParameter(pParameter), m_pFr(NULL), m_pProIndex(NULL), m_tProteinsNum(
				0), m_tProteinIdKeeper(0) {
	try {
		if (m_pFr == NULL) {
			m_pFr = new FastaReader(m_pParameter->m_strFastaFilePath.c_str());
		}

		countProteinsNum();

		if (m_pFr != NULL) {
			m_pFr->getConnection();
			m_pFr->first();
		}

		if (m_pProIndex == NULL) {
			m_pProIndex = new ProteinIndex(m_pParameter);
		}
	} catch (exception &e) {
		ErrorInfo err("IonIndexer", "IonIndexer", "caught an exception.", e);
		throw runtime_error(err.get());
	} catch (...) {
		ErrorInfo err("IonIndexer", "IonIndexer",
				"caught an unknown exception.");
		throw runtime_error(err.get());
	}
}

IonIndexer::~IonIndexer() {
	if (m_pFr != NULL) {
		delete m_pFr;
		m_pFr = NULL;
	}

	if (m_pProIndex != NULL) {
		delete m_pProIndex;
		m_pProIndex = NULL;
	}
}

bool IonIndexer::hasNextLoad() {
	return (m_tProteinIdKeeper < m_tProteinsNum);
}

size_t IonIndexer::loadNext() {
	size_t tOneLoad = m_pParameter->m_tProteinsBatchSize;
	size_t tProteinId = m_tProteinIdKeeper;
	m_tProteinIdKeeper =
			(tProteinId + tOneLoad < m_tProteinsNum) ?
					(tProteinId + tOneLoad) : m_tProteinsNum;
	try {
		m_pProIndex->clear();
		for (size_t tId = tProteinId; tId < m_tProteinIdKeeper; ++tId) {
			m_pProIndex->attachProtein(m_pFr->currentItem());
			m_pFr->next();
		}
		m_pProIndex->endAttach();
//		m_pProIndex->showIndex();
	} catch (exception &e) {
		ErrorInfo err("IonIndexer", "loadNext", "caught an exception.", e);
		throw runtime_error(err.get());
	} catch (...) {
		ErrorInfo err("IonIndexer", "loadNext", "caught an unknown exception.");
		throw runtime_error(err.get());
	}
	return m_tProteinIdKeeper - tProteinId;
}

size_t IonIndexer::getPeptideSequence(string &strSq, size_t tStart,
		size_t tLen) {
	return m_pProIndex->getPeptideSequence(strSq, tStart, tLen);
}

size_t IonIndexer::getPeptideSequence(char *szSequence, size_t tStart,
		size_t tLen) {
	return m_pProIndex->getPeptideSequence(szSequence, tStart, tLen);
}

double IonIndexer::getPeptideMass(size_t tStart, size_t tLen) {
	return m_pProIndex->getPeptideMass(tStart, tLen);
}

unsigned char IonIndexer::getPeptideType(size_t tStart, size_t tLen) {
	unsigned char cType = 0;
	try {
		size_t tIdx = m_pProIndex->atProtein(tStart);
		if (tIdx == ALL_BIT_ON) {
			throw runtime_error("Fail to locate the protein.");
		}
		size_t tProN = m_pProIndex->at(tIdx);
		size_t tProC = m_pProIndex->at(tIdx + 1);
		if(tStart == tProN || (tStart == tProN + 1 && m_pProIndex->isM(tProN))) {
			cType = (tStart + tLen == tProC) ? 3 : 1;
		} else if (tStart + tLen == tProC) {
			cType = 2;
		} else {
			cType = 0;
		}
	} catch (exception &e) {
		ErrorInfo err("IonIndexer", "getPeptideType", "caught an exception.",
				e);
		throw runtime_error(err.get());
	} catch (...) {
		ErrorInfo err("IonIndexer", "getPeptideType",
				"caught an unknown exception.");
		throw runtime_error(err.get());
	}
	return cType;
}

size_t IonIndexer::getMissSiteNum(size_t tStart, size_t tLen) {
	size_t tNum = 0;
	try {
		tNum = m_pProIndex->getMissSiteNum(tStart, tLen);
	} catch (exception &e) {
		ErrorInfo err("IonIndexer", "getMissSiteNum", "caught an exception.",
				e);
		throw runtime_error(err.get());
	} catch (...) {
		ErrorInfo err("IonIndexer", "getMissSiteNum",
				"caught an unknown exception.");
		throw runtime_error(err.get());
	}
	return tNum;
}

size_t IonIndexer::getProteinsNum() const {
	return m_tProteinsNum;
}

/* begin: internal functions */

void IonIndexer::countProteinsNum() {
	try {
		m_pFr->getConnection();
		for (m_pFr->first(); m_pFr->hasNext(); m_pFr->next())
			++m_tProteinsNum;
		m_pFr->endConnection();
	} catch (exception &e) {
		ErrorInfo err("IonIndexer", "_countProteinsNum", "caught an exception.",
				e);
		throw runtime_error(err.get());
	} catch (...) {
		ErrorInfo err("IonIndexer", "_countProteinsNum",
				"caught an unknown exception.");
		throw runtime_error(err.get());
	}
}


EnumPeptideIndexer::EnumPeptideIndexer(const SearchParameter *pParameter) :
		IonIndexer(pParameter), m_pPepIndex(NULL) {
	try {
		if (m_pPepIndex == NULL) {
			m_pPepIndex = new EnumPeptideIndex(m_pParameter);
		}
	} catch (exception &e) {
		ErrorInfo err("EnumPeptideIndexer", "EnumPeptideIndexer",
				"caught an exception.", e);
		throw runtime_error(err.get());
	} catch (...) {
		ErrorInfo err("EnumPeptideIndexer", "EnumPeptideIndexer",
				"caught an unknown exception.");
		throw runtime_error(err.get());
	}
}

EnumPeptideIndexer::~EnumPeptideIndexer() {
	if (m_pPepIndex != NULL) {
		delete m_pPepIndex;
		m_pPepIndex = NULL;
	}
}

size_t EnumPeptideIndexer::loadNext() {
	size_t tLoadedNum(0);
	try {
		tLoadedNum = IonIndexer::loadNext();
		m_pPepIndex->clear();
		m_pPepIndex->attachPeptides(0, tLoadedNum, m_pProIndex);
		//m_pPepIndex->showPeptideTable(m_pProIndex);
	} catch (exception &e) {
		ErrorInfo err("EnumPeptideIndexer", "loadNext", "caught an exception.", e);
		throw runtime_error(err.get());
	} catch (...) {
		ErrorInfo err("EnumPeptideIndexer", "loadNext",
				"caught an unknown exception.");
		throw runtime_error(err.get());
	}
	return tLoadedNum;
}

void EnumPeptideIndexer::getPeptides(std::vector<Peptide*> &vPeps){
	try {
		m_pPepIndex->getPeptides(vPeps);
	} catch (exception &e) {
		ErrorInfo err("EnumPeptideIndexer", "getPeptides",
				"caught an exception.", e);
		throw runtime_error(err.get());
	} catch (...) {
		ErrorInfo err("EnumPeptideIndexer", "getPeptides",
				"caught an unknown exception.");
		throw runtime_error(err.get());
	}
}

void EnumPeptideIndexer::getXLinkPeptides(std::vector<Peptide*> &vPeps){
	try {
		m_pPepIndex->getXLinkPeptides(vPeps);
	} catch (exception &e) {
		ErrorInfo err("EnumPeptideIndexer", "getPeptides",
				"caught an exception.", e);
		throw runtime_error(err.get());
	} catch (...) {
		ErrorInfo err("EnumPeptideIndexer", "getPeptides",
				"caught an unknown exception.");
		throw runtime_error(err.get());
	}
}

/* end: internal function */

}
