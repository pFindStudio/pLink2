#include "../include/sdk.h"
#include "../include/index.h"

using namespace std;

namespace sdk
{

map<string, AccessionIndexer *> AccessionIndexer::m_mpACIndexers;

AccessionIndexer *AccessionIndexer::getAccessionIndex(const std::string &strFastaFilePath)
{
	if(m_mpACIndexers.find(strFastaFilePath) == m_mpACIndexers.end()) {
		AccessionIndexer *pACIndexer = new AccessionIndexer(strFastaFilePath);
		m_mpACIndexers.insert(pair<string, AccessionIndexer *>(strFastaFilePath, pACIndexer));
	}
	return m_mpACIndexers[strFastaFilePath];
}

void AccessionIndexer::destroyAccessionIndex()
{
	map<string, AccessionIndexer *>::iterator it = m_mpACIndexers.begin();
	while(it != m_mpACIndexers.end()) {
		if(it->second != NULL) {
			delete it->second;
			it->second = NULL;
		}
		++it;
	}
	m_mpACIndexers.clear();
}

AccessionIndexer::AccessionIndexer(const string &strFastaFilePath) :
		m_tProteinNum(0)
{
	createAccessionIndex(strFastaFilePath);
}

void AccessionIndexer::createAccessionIndex(const string &strFastaFilePath)
{
	FastaReader *pFr = new FastaReader(strFastaFilePath.c_str());
	pFr->getConnection();

	size_t tCurStart(0);
	for(pFr->first(); pFr->hasNext(); pFr->next()) {
		ProteinItem stItem = pFr->currentItem();
		m_strAcs += stItem.strAc;
		m_vPrtnEntr.push_back(tCurStart);
		tCurStart += stItem.strAc.length();
		++m_tProteinNum;
	}
	m_vPrtnEntr.push_back(tCurStart);

	pFr->endConnection();
	delete pFr;
}

string AccessionIndexer::getACByID(size_t tProteinId)
{
	if(tProteinId >= m_tProteinNum)
		return "";

	return m_strAcs.substr(m_vPrtnEntr[tProteinId],
			m_vPrtnEntr[tProteinId+1]-m_vPrtnEntr[tProteinId]);
}

size_t AccessionIndexer::getProteinNum() const
{
	return m_tProteinNum;
}

}
