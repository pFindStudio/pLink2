
#ifndef SDK_H_
#define SDK_H_
#include "util.h"

namespace sdk
{

/* begin: classes */

struct IonPeak {
	double m_lfMz;
	double m_lfIntensity;
	int m_nMz;
	size_t m_tCharge;

	IonPeak(double lfMz = 0.0, double lfIntensity = 0.0, size_t tCharge = 0);
	IonPeak(const IonPeak &stPeak);
	IonPeak &operator=(const IonPeak &stPeak);
};

struct QueryPeak {
	size_t m_tLowIdx;
	size_t m_tHighIdx;

	QueryPeak(size_t tLowIdx = 0, size_t tHighIdx = 0);
};

struct Spectrum {
	std::string m_strTitle;
	double m_lfMH;
	double m_lfMZ;
	double m_lfTotalInt;
	double m_lfSqrtMaxInt;
	size_t m_tCharge;
    size_t m_tScanNo;
	std::vector<IonPeak> m_vPeaks;
	std::vector<int> m_vHashIndex;

	static bool intesityGreater(const IonPeak &stPeak1, const IonPeak &stPeak2);
	static bool intesityLesser(const IonPeak &stPeak1, const IonPeak &stPeak2);
	static bool mzGreater(const IonPeak &stPeak1, const IonPeak &stPeak2);
	static bool mzLesser(const IonPeak &stPeak1, const IonPeak &stPeak2);

	Spectrum();

	friend
	std::ostream &operator<<(std::ostream &stOut, const Spectrum &stSpec);
	friend
	std::fstream &operator<<(std::fstream &stOut, const Spectrum &stSpec);
	friend
	std::istream &operator>>(std::istream &stIn, Spectrum &stSpec);
	friend
	std::fstream &operator>>(std::fstream &stIn, Spectrum &stSpec);
	void copyBasicItems(const Spectrum &stSpec);
	void createHashIndex();
};

struct QuerySpectrum {
	std::vector<QueryPeak> m_vQueryBPeaks;
	std::vector<QueryPeak> m_vQueryYPeaks;

	void swapOut();
};

// This structure can be used to keep modification and linker information
struct ModificationEntry {
	int m_nSite; // position
	int m_nId; // which modification or which linker

	ModificationEntry();
	ModificationEntry(int nSite, int nId);

	friend
	bool operator==(const ModificationEntry &stEntry1, const ModificationEntry &stEntry2)
	{
		return stEntry1.m_nSite == stEntry2.m_nSite && stEntry1.m_nId == stEntry2.m_nId;
	}
	friend
	bool operator!=(const ModificationEntry &stEntry1, const ModificationEntry &stEntry2)
	{
		return !(stEntry1.m_nSite == stEntry2.m_nSite && stEntry1.m_nId == stEntry2.m_nId);
	}
	friend
	std::ostream &operator<<(std::ostream &stOut, const ModificationEntry &stEntry);
	friend
	std::fstream &operator<<(std::fstream &stOut, const ModificationEntry &stEntry);
	friend
	std::istream &operator>>(std::istream &stIn, ModificationEntry &stEntry);
	friend
	std::fstream &operator>>(std::fstream &stIn, ModificationEntry &stEntry);

	void clear();
	std::string toDebugString();
};

struct PeptideItem {
	size_t m_tCount;
	size_t m_tPepId;

	PeptideItem();
	PeptideItem(size_t, size_t);

	friend
	bool operator<(const PeptideItem &item1, const PeptideItem &item2)
	{
		return item1.m_tCount < item2.m_tCount;
	}
	friend
	bool operator>(const PeptideItem &stItem1, const PeptideItem &stItem2)
	{
		return stItem1.m_tCount > stItem2.m_tCount;
	}
};

struct Peptide {
	std::string m_strSq;
	double m_lfMass; // residue mass
	std::vector<ModificationEntry> m_vMods;
	unsigned char m_ucMiss;
	unsigned char m_ucEnd;
	double m_lfTag; //GodelCode

	static bool massGreater(const Peptide *pPep1, const Peptide *pPep2);
	static bool massLesser(const Peptide *pPep1, const Peptide *pPep2);

	static bool massTagGreater(const Peptide *pPep1, const Peptide *pPep2);
	static bool massTagLesser(const Peptide *pPep1, const Peptide *pPep2);

	Peptide();

	friend bool operator==(const Peptide &stPep1, const Peptide &stPep2);
	friend bool operator!=(const Peptide &stPep1, const Peptide &stPep2);
	friend std::ostream &operator<<(std::ostream &stOut, const Peptide &stPeptide);
	friend std::fstream &operator<<(std::fstream &stOut, const Peptide &stPeptide);
	friend std::istream &operator>>(std::istream &stIn, Peptide &stPeptide);
	friend std::fstream &operator>>(std::fstream &stIn, Peptide &stPeptide);

	bool isProNTerm() const;
	bool isProCTerm() const;
	void setPeptideInfo(const char *szSq, double lfMass, unsigned char ucMiss, unsigned char ucEnd);
	void setPeptideInfo(const std::string &strSq, double lfMass, unsigned char ucMiss, unsigned char ucEnd);
	void clear();
	std::string toString();
	std::string toDebugString();
	double getGodelCode() const;
	double calculatePeptideMass();
};

struct ProteinInformation {
	std::vector<size_t> m_tProteinID;
	std::vector<size_t> m_tProteinSite;

	ProteinInformation() {}

	friend std::ostream &operator<<(std::ostream &stOut, const ProteinInformation &stInfo);
	friend std::fstream &operator<<(std::fstream &stOut, const ProteinInformation &stInfo);
	friend std::istream &operator>>(std::istream &stIn, ProteinInformation &stInfo);
	friend std::fstream &operator>>(std::fstream &stIn, ProteinInformation &stInfo);

	std::string toString();
};

struct PeptideResult {
	Peptide m_stPep; // peptide sequence information
	size_t m_tWinId; // peptide matched window's id
	double m_lfOpenMass; // extra open mass of whole peptide
	double m_lfScore; // pre-score value
	ModificationEntry m_stLinkMod; // crosslinking site information
	ProteinInformation *m_pProteinInfo;

	static bool scoreGreater(const PeptideResult &stResult1, const PeptideResult &stResult2);
	static bool scoreLesser(const PeptideResult &stResult1, const PeptideResult &stResult2);

	PeptideResult();
	PeptideResult(const PeptideResult &stResult);
	~PeptideResult();
	PeptideResult &operator=(const PeptideResult &stResult);

	friend bool operator==(const PeptideResult &stResult1, const PeptideResult &stResult2);
	friend bool operator!=(const PeptideResult &stResult1, const PeptideResult &stResult2);
	friend std::ostream &operator<<(std::ostream &stOut, const PeptideResult &stResult);
	friend std::fstream &operator<<(std::fstream &stOut, const PeptideResult &stResult);
	friend std::istream &operator>>(std::istream &stIn, PeptideResult &stResult);
	friend std::fstream &operator>>(std::fstream &stIn, PeptideResult &stResult);

	bool isTarget(const char *szDecoySign) const;
	void clear();
	std::string toString();
	std::string toDebugString();
};

struct RatioFeature {
	int m_nLength;
	int m_nTagLen;
	double m_lfMatchedPeakRatio;
	double m_lfMatchedIonRatio;

	RatioFeature() :
		m_nLength(0), m_nTagLen(0), m_lfMatchedPeakRatio(0.0),
		m_lfMatchedIonRatio(0.0)
	{
	}

	friend std::ostream &operator<<(std::ostream &stOut, const RatioFeature &stInfo);
	friend std::fstream &operator<<(std::fstream &stOut, const RatioFeature &stInfo);
	friend std::istream &operator>>(std::istream &stIn, RatioFeature &stInfo);
	friend std::fstream &operator>>(std::fstream &stIn, RatioFeature &stInfo);

    void clear();
};

struct XLinkFeature {
	RatioFeature m_stAlpha;
	RatioFeature m_stBeta;

	friend std::ostream &operator<<(std::ostream &stOut, const XLinkFeature &stInfo);
	friend std::fstream &operator<<(std::fstream &stOut, const XLinkFeature &stInfo);
	friend std::istream &operator>>(std::istream &stIn, XLinkFeature &stInfo);
	friend std::fstream &operator>>(std::fstream &stIn, XLinkFeature &stInfo);

    void clear();
};

struct FeatureInformation {
	RatioFeature m_stRatio;
	double m_lfErrorMean;
	double m_lfErrorDeviation;
	XLinkFeature m_stXLink;

	FeatureInformation();

	friend std::ostream &operator<<(std::ostream &stOut, const FeatureInformation &stInfo);
	friend std::fstream &operator<<(std::fstream &stOut, const FeatureInformation &stInfo);
	friend std::istream &operator>>(std::istream &stIn, FeatureInformation &stInfo);
	friend std::fstream &operator>>(std::fstream &stIn, FeatureInformation &stInfo);

	std::string toString();
	std::string toDebugString();
    void clear();
};

struct ReportInformation {
	double m_lfMatchedPercentage;
	double m_lfTotalIntensity;
	size_t m_tMatchedNumber;
	size_t m_tTotalNumber;
	std::vector<bool> m_bNContinuity;
	std::vector<bool> m_bCContinuity;
	std::vector<int> m_nNCoexistence;
	std::vector<int> m_nCCoexistence;
	std::vector<int> m_vIonNumber;
	std::vector<int> m_vIonMatchedNumber;

	ReportInformation();

	friend std::ostream &operator<<(std::ostream &stOut, const ReportInformation &stInfo);
	friend std::fstream &operator<<(std::fstream &stOut, const ReportInformation &stInfo);
	friend std::istream &operator>>(std::istream &stIn, ReportInformation &stInfo);
	friend std::fstream &operator>>(std::fstream &stIn, ReportInformation &stInfo);

	std::string toString();
	std::string toDebugString();
};

struct XLinkPeptideItem {
	std::string m_strSq;
	size_t m_tSpecIdx;
	size_t m_tRankIdx;
	bool m_bAlpha;

	XLinkPeptideItem(std::string strSq = "", size_t tSpecIdx = 0,
			size_t tRankIdx = 0, bool bAlpha = true);
};

struct XLinkPeptideResult {
	PeptideResult m_stAlphaPep;
	PeptideResult m_stBetaPep;
	PeptideType m_ePepType;
	ProteinType m_eProType;
	TDType m_eTDType;
	double m_lfScore;
	double m_lfEvalue;
	double m_lfQvalue;
	double m_lfCalcMH;
	FeatureInformation m_stFeatureInfo;

	XLinkPeptideResult();
	~XLinkPeptideResult();

	friend std::ostream &operator<<(std::ostream &stOut, const XLinkPeptideResult &stResult);
	friend std::fstream &operator<<(std::fstream &stOut, const XLinkPeptideResult &stResult);
	friend std::istream &operator>>(std::istream &stIn, XLinkPeptideResult &stResult);
	friend std::fstream &operator>>(std::fstream &stIn, XLinkPeptideResult &stResult);
	friend bool operator==(const XLinkPeptideResult &stResult1, const XLinkPeptideResult &stResult2);
	friend bool operator!=(const XLinkPeptideResult &stResult1, const XLinkPeptideResult &stResult2);
	friend bool operator<(const XLinkPeptideResult &stResult1, const XLinkPeptideResult &stResult2);
	friend bool operator>(const XLinkPeptideResult &stResult1, const XLinkPeptideResult &stResult2);

	void setProteinTDType(const char *szDecoySign);
	int unifiedPosition(int nIdx, bool bAlpha);
	void getModifications(std::vector<std::string> &vMods);
	void getModifications(std::string &strMods);
	void getModifications(std::set<size_t> &vModTypes);
	void setPeptideSq(std::string &strSq);
	void setPeptideLinkedId(std::string &strLinker);
	void setPeptideType(std::string &strType);
	void setModifications(std::string &strMods);
	void clear();
	std::string toString();
	std::string toDebugString();
	double getGodelCode() const;
	double calculateXLinkPeptideMass();

private:
	ProteinType getProteinType(const char *szDecoySign) const;
	int getModificationId(std::string &strMod);
};

struct OpenMatchResult {
	std::vector<PeptideItem> m_vPepItemsQ;
	std::vector<PeptideResult> m_vPepsQ;

	OpenMatchResult(size_t tSize);
	~OpenMatchResult();

	void buildMinHeap(size_t tCount, const PeptideResult &stResult);
	void heapify(size_t tCount, const PeptideResult &stResult);
	void heapify(int nIdx);
	size_t getSize() const
	{
		return m_vPepItemsQ.size();
	}

	size_t topCount() const
	{
		return m_vPepItemsQ[0].m_tCount;
	}
	void removeRedundant();
	void swapOut();
};



class AccessionIndexer {
	static std::map<std::string, AccessionIndexer *> m_mpACIndexers;

	size_t m_tProteinNum;
	std::vector<size_t> m_vPrtnEntr; /* protein index */
	std::string m_strAcs; /* protein AC */

public:
	static AccessionIndexer *getAccessionIndex(const std::string &strFastaFilePath);
	static void destroyAccessionIndex();

	// get a protein AC by IDs
	std::string getACByID(size_t tProteinId);
	size_t getProteinNum() const;

private:
	// need m_strFastaFilePath
	AccessionIndexer(const std::string &strFastaFilePath);
	void createAccessionIndex(const std::string &strFastaFilePath);
};

/* struct GroupCounter
 *   count psm for each group
 */
struct GroupCounter {
	size_t m_vPepCounts[PEPTIDE_TYPE_NUM];
	size_t m_vProCounts[PROTEIN_TYPE_NUM];
	size_t m_vChargeCounts[CHARGE_RANGE];
	size_t m_vErrorCounts[PRECURSOR_ERROR_INTERVAL_NUM];
	size_t m_vModCounts[MAX_MODIFY_NUM];
	double m_lfSystemaicError;

	GroupCounter();
	void clear();
};

/* struct XLinkPSMFeature
 *   features extracted from PSMs
 */
struct XLinkPSMFeature {
	double m_lfScore; // SVM score
	double m_lfEvalue; // 1 / (1+exp(SVM Score))
	std::vector<double> m_vFeatures; // features for SVM

	XLinkPSMFeature();

	friend bool operator<(const XLinkPSMFeature &stItem1, const XLinkPSMFeature &stItem2);
	friend bool operator<=(const XLinkPSMFeature &stItem1, const XLinkPSMFeature &stItem2);
	friend bool operator>(const XLinkPSMFeature &stItem1, const XLinkPSMFeature &stItem2);
	friend bool operator>=(const XLinkPSMFeature &stItem1, const XLinkPSMFeature &stItem2);
};

/* struct PSMSortItem
 *   top-1 PSM of each spectrum
 */
struct PSMSortItem {
	double m_lfScore; // SVM score
	double m_lfQvalue; // PSM q-value
	size_t m_tIndex; // PSM identity

	PSMSortItem(double lfScore = 0.0, size_t tIndex = 0) :
		m_lfScore(lfScore), m_lfQvalue(0.0), m_tIndex(tIndex)
	{
	}

	friend bool operator<(const PSMSortItem &stItem1, const PSMSortItem &stItem2);
	friend bool operator<=(const PSMSortItem &stItem1, const PSMSortItem &stItem2);
	friend bool operator>(const PSMSortItem &stItem1, const PSMSortItem &stItem2);
	friend bool operator>=(const PSMSortItem &stItem1, const PSMSortItem &stItem2);
};

struct ReportSummary {
	size_t m_tSpectraTotal; // all spectra
	size_t m_vSpectraNum[PEPTIDE_TYPE_NUM]; // all spectra have identification results, classified by peptide type
	size_t m_tSpecTotal; // filtered by FDR, positive
	size_t m_vSpecNum[PEPTIDE_TYPE_NUM]; // spectra filtered by FDR, positive, classified by peptide type
	size_t m_vPepNum[PEPTIDE_TYPE_NUM]; // peptides filtered by FDR, positive, classified by peptide type
	size_t m_vSitesNum[PEPTIDE_TYPE_NUM]; // sites filtered by FDR, classified by peptide type
	size_t m_vGroupNum[PEPTIDE_TYPE_NUM];

	ReportSummary();
};

template<typename Invoke, typename Func>
class ModificationDecorator;

class XLinkDecorator;

/* end: classes */

} /* end of namespace sdk */

#endif /* SDK_H_ */
