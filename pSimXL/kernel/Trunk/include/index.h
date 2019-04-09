
#ifndef IONINDEXER_H_
#define IONINDEXER_H_
#include "util.h"
#include "sdk.h"

namespace sdk
{

/* begin: structures, to save space, attention: not to contain methods */

struct PeptideEntry {
	size_t tProStart;
	size_t tProLen;
	size_t tPepStart;
	size_t tPepLen;
};

struct ProteinItem {
	std::string strAc;	/* Accession number*/
	std::string strDe;	/* Description */
	std::string strSq;	/* Protein sequence */
};

/* end: structures */

/* begin: classes */

class FastaReader {
	std::string m_strFastaFilePath;
	std::ifstream m_fin;
	ProteinItem m_current;
	ProteinItem m_next;

public:
	FastaReader();
	FastaReader(const char *szFastaFilePath);
	~FastaReader();
	void setFastaFilePath(const char *szFastaFilePath);
	void getConnection() throw(std::runtime_error);
	void endConnection();
	void first();
	void next();
	bool hasNext();
	ProteinItem currentItem() const;

	static void testCase();
};



class ProteinIndex {
	const SearchParameter *m_pParameter;
	size_t m_tProteinNum;
	std::vector<size_t> m_vPrtnEntr; /* protein index */
	std::string m_strSq; /* protein sequences */
	std::vector<double> m_vSqMass; /* cumulative mass */
	std::vector<size_t> m_vSiteEntr; /* sites index */
	std::vector<size_t> m_vSite; /* available cleave sites */
	bool m_bAttachFinished;

public:
	friend class EnumPeptideIndex;
	friend class PeptideIndex;
	friend class FragmentIndex;
	ProteinIndex(const SearchParameter *pParameter);
	~ProteinIndex();

	void attachProtein(const ProteinItem &stItem);
	void endAttach();
	void clear();
	bool isAttachFinished() const;
	bool isM(size_t tProteinStart) const;
	size_t getPeptideSequence(std::string &strSq, size_t tStart, size_t tLen);
	size_t getPeptideSequence(char *szSequence, size_t tStart, size_t tLen);
	double getPeptideMass(size_t tStart, size_t tLen);
	size_t getMissSiteNum(size_t tStart, size_t tLen);
	size_t at(size_t tProteinId) const;
	size_t atProtein(size_t ) const;
	size_t findEndSite(size_t tSiteIdxStart, size_t tSiteIdxEnd, size_t tPos) const;
	bool isCleaveSite(size_t tSiteIdxStart, size_t tSiteIdxEnd, size_t tPos) const;
	void showIndex();
};

class EnumPeptideIndex {
	const SearchParameter *m_pParameter;

	ModificationDecorator<EnumPeptideIndex, EnumPeptideIndexDispatch> *m_pModDecorator;
	XLinkDecorator *m_pXLinkDecorator;
	size_t m_tMassRange;
	std::vector<size_t> m_vCnt;
	std::vector<Peptide> m_vPepTable;

	Peptide m_stCurPep;

public:
	EnumPeptideIndex(const SearchParameter *pParameter);
	~EnumPeptideIndex();
	void attachPeptides(size_t tStartId, size_t tEndId,
			const ProteinIndex *pProIndex);
	void clear();
	void getPeptides(std::vector<Peptide*> &vPeps);
	void getXLinkPeptides(std::vector<Peptide*> &vPeps);

protected:
	void attachSpecific(size_t tProteinId, const ProteinIndex *pProIndex, std::vector<Peptide> &vPepTable);
	void attachSemiSpecific(size_t tProteinId, const ProteinIndex *pProIndex, std::vector<Peptide> &vPepTable);
	void attachNonSpecific(size_t tProteinId, const ProteinIndex *pProIndex, std::vector<Peptide> &vPepTable);
	void attachPTMPeptides();
	void getUsefulPep(std::vector<Peptide> &vPepTable);
	unsigned char getPeptideType(const PeptideEntry &stEntry);
};


/* This is the whole work-flow of ion indexer; it is just a proxy */
class IonIndexer {
protected:
	const SearchParameter *m_pParameter;
	FastaReader *m_pFr;
	ProteinIndex *m_pProIndex;
	size_t m_tProteinsNum;
	size_t m_tProteinIdKeeper;

public:
	IonIndexer(const SearchParameter *pParameter);
	virtual ~IonIndexer();

	virtual bool hasNextLoad();
	virtual size_t loadNext();
	virtual size_t getProteinsNum() const;

	virtual size_t getPeptideSequence(std::string &strSq, size_t tStart, size_t tLen);
	virtual size_t getPeptideSequence(char *szSequence, size_t tStart, size_t tLen);
	virtual double getPeptideMass(size_t tStart, size_t tLen);
	virtual unsigned char getPeptideType(size_t tStart, size_t tLen);
	virtual size_t getMissSiteNum(size_t tStart, size_t tLen);

protected:
	virtual void countProteinsNum();
};



class EnumPeptideIndexer : public IonIndexer {
	EnumPeptideIndex *m_pPepIndex;

public:
	EnumPeptideIndexer(const SearchParameter *pParameter);
	virtual ~EnumPeptideIndexer();

	virtual size_t loadNext();
	void getPeptides(std::vector<Peptide*> &vPeps);
	void getXLinkPeptides(std::vector<Peptide*> &vPeps);
};

/* end: classes */

}
#endif /* IONINDEXER_H_ */
