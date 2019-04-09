#ifndef SCORES_H_
#define SCORES_H_
#include "util.h"
#include "sdk.h"
#include "bias.h"
namespace sdk
{

struct MzTriple {
	int nMz;
	double lfMz;
	byte yIonTypeOrder;
	byte yPepPosOrder1;
	byte yPepPosOrder2;
	byte yNTerm1;
	byte yNTerm2;
	byte yAddToTag;
	byte yContainLinker;

	MzTriple();

	friend bool operator<(const MzTriple &a, const MzTriple &b)
	{
		return a.lfMz < b.lfMz;
	}

	friend bool operator<=(const MzTriple &a, const MzTriple &b)
	{
		return a.lfMz <= b.lfMz;
	}

	friend bool operator>(const MzTriple &a, const MzTriple &b)
	{
		return a.lfMz > b.lfMz;
	}

	friend bool operator>=(const MzTriple &a, const MzTriple &b)
	{
		return a.lfMz >= b.lfMz;
	}

	friend std::ostream &operator<<(std::ostream &stOut, const MzTriple &stTriple);

	void clear();
	std::string toString();
};

struct PeptideAAMass {
	std::vector<double> m_nPepAAMass;
	std::vector<double> m_bPepAAMass;
	std::vector<double> m_yPepAAMass;
	double m_lfPepMass;
	int m_nPepLen;

	PeptideAAMass();
	~PeptideAAMass();

	void resize(size_t tSize);
	void clear();
};

class MzCalculator {
	double *m_pAAMass;
	SearchParameter *m_pParameter;
	std::vector<Modification> m_vMods;
	Peptide *m_pPeptide;
	PeptideResult *m_pPeptideResult;
	XLinkPeptideResult *m_pXLinkPeptideResult;
	PeptideAAMass *m_pAlphaPep;
	PeptideAAMass *m_pBetaPep;

public:
	MzCalculator();
	virtual ~MzCalculator();

	virtual void init(Peptide *pPeptide);
	virtual void init(PeptideResult *pPeptide);
	virtual void init(XLinkPeptideResult *pPeptide);

	virtual void attachNTermIons(std::vector<MzTriple> &vIonMz, const IonType &stType, byte yOrder);
	virtual void attachCTermIons(std::vector<MzTriple> &vIonMz, const IonType &stType, byte yOrder);

protected:
	virtual void close();
};

class MzCalculatorWithBias : public MzCalculator {
	double *m_pAAMass;
	SearchParameter *m_pParameter;
	std::vector<Modification> m_vMods;
	Peptide *m_pPeptide;
	PeptideResult *m_pPeptideResult;
	XLinkPeptideResult *m_pXLinkPeptideResult;
	PeptideAAMass *m_pAlphaPep;
	PeptideAAMass *m_pBetaPep;

public:
	MzCalculatorWithBias();
	virtual ~MzCalculatorWithBias();
	virtual void attachNTermIons(std::vector<MzTriple> &vIonMz, const IonType &stType, byte yOrder);
	virtual void attachCTermIons(std::vector<MzTriple> &vIonMz, const IonType &stType, byte yOrder);
	virtual void attachInternalIons(std::vector<MzTriple> &vIonMz, const IonType &stType, byte yOrder);
	virtual void attachPrecursorIons(std::vector<MzTriple> &vIonMz, const IonType &stType, byte yOrder);
	virtual void attachKLIons(std::vector<MzTriple> &vIonMz, const IonType &stType, byte yOrder);
	virtual void attachYPIons(std::vector<MzTriple> &vIonMz, const IonType &stType, byte yOrder);
	virtual void attachBPIons(std::vector<MzTriple> &vIonMz, const IonType &stType, byte yOrder);

protected:
	virtual void close();
};


/* class RefinedScorer
 *   this class implements theoretical peaks generation and
 *   peptide pair spectrum matching.
 *   m_vIonMz: stores generated theoretical peaks
 *   m_vMatchedInfo: stores matched information,
 *     for each spectrum peak, the structure keeps those matched theoretical peaks' index in m_vIonMz
 */
class RefinedScorer {
protected:
	const SearchParameter *m_pParameter;
	Spectrum *m_pSpec;
	XLinkPeptideResult *m_pPeptide;
	MzCalculator *m_pCalculator;
	std::vector<MzTriple> m_vIonMz;
	std::vector<std::vector<int> > m_vMatchedInfo;
	bool m_bComputeMz;

public:
	RefinedScorer(const SearchParameter *pParameter);
	virtual ~RefinedScorer();

	virtual void setSpectrum(Spectrum *pSpec);
	virtual void setPeptide(XLinkPeptideResult *pPeptide);
	virtual double score() = 0;

protected:
	virtual void computeMz();
	virtual void match();
	virtual void getMassBorder(int &nMin, int &nMax, int nMz);
	virtual double getError(size_t tThrMz, size_t tExpMz);
};

class RefinedScorerFactory {
	const SearchParameter *m_pParameter;

public:
	RefinedScorerFactory(const SearchParameter *pParameter);
	virtual ~RefinedScorerFactory() {}

	RefinedScorer *getScorer(ScoreType eScoreType) const;

};


class BMIonGenerator : public RefinedScorer {
private:
	const static std::string strTitlePrefix;
	const static std::string strTitleSuffix;
    static size_t m_tScanNo;
    Bias *m_pBias;
    std::vector<char> m_vBasicAA;
    std::vector<IonPeak> m_vPeaks;
    std::map<std::pair<char, int>, bool> m_mpIonExsist;
    std::map<std::pair<char, int>, double> m_mpIonRatio;
    std::map<std::pair<char, int>, double> m_mpIonIntensity;
public:
	BMIonGenerator(const SearchParameter *pParameter);
	virtual ~BMIonGenerator();
	virtual double score();
protected:
	virtual void match();
	virtual size_t generateCharge();
	virtual double getPepTypeScore(PeptideType ePepType);
	void insertRatio(std::map<std::pair<char, int>, double> &mpStatistic, char cSymbol, int nCharge, double lfRatio);
	void insertRatio(std::map<std::pair<char, int>, bool> &mpStatistic, char cSymbol, int nCharge, double lfRatio);
	void genMatchRatio();
	void genIntensityMean();
	void genExistFlag();
	double getMatchRatio(IonType stIonType);
	double getIntensityMean(IonType stIonType);
	bool getExistFlag(IonType stIonType);
//	virtual double getBasicAANum(XLinkPeptideResult &stXLinkPeptideResult);
	virtual void mergePeaks();
};

class BMIonPerfectGenerator : public RefinedScorer {
private:
	const static std::string strTitlePrefix;
	const static std::string strTitleSuffix;
    static size_t m_tScanNo;
    Bias *m_pBias;
    std::vector<char> m_vBasicAA;
    std::vector<IonPeak> m_vPeaks;
    std::map<std::pair<char, int>, bool> m_mpIonExsist;
    std::map<std::pair<char, int>, double> m_mpIonRatio;
    std::map<std::pair<char, int>, double> m_mpIonIntensity;
public:
    BMIonPerfectGenerator(const SearchParameter *pParameter);
	virtual ~BMIonPerfectGenerator();
	virtual double score();
protected:
	virtual void match();
	virtual size_t generateCharge();
	virtual double getPepTypeScore(PeptideType ePepType);
	void insertRatio(std::map<std::pair<char, int>, double> &mpStatistic, char cSymbol, int nCharge, double lfRatio);
	void insertRatio(std::map<std::pair<char, int>, bool> &mpStatistic, char cSymbol, int nCharge, double lfRatio);
	void genMatchRatio();
	void genIntensityMean();
	void genExistFlag();
	double getMatchRatio(IonType stIonType);
	double getIntensityMean(IonType stIonType);
	bool getExistFlag(IonType stIonType);
//	virtual double getBasicAANum(XLinkPeptideResult &stXLinkPeptideResult);
	virtual void mergePeaks();
};

class BMIonSimpleGenerator : public RefinedScorer {
private:
	const static std::string strTitlePrefix;
	const static std::string strTitleSuffix;
    static size_t m_tScanNo;
    Bias *m_pBias;
    std::vector<char> m_vBasicAA;
    std::vector<IonPeak> m_vPeaks;
//    std::map<std::pair<char, int>, bool> m_mpIonExsist;
    std::map<std::pair<char, int>, double> m_mpIonRatio;
    std::map<std::pair<char, int>, double> m_mpIonIntensity;
public:
    BMIonSimpleGenerator(const SearchParameter *pParameter);
	virtual ~BMIonSimpleGenerator();
	virtual double score();
protected:
	virtual void match();
	virtual size_t generateCharge();
	void insertRatio(std::map<std::pair<char, int>, double> &mpStatistic, char cSymbol, int nCharge, double lfRatio);
	void genMatchRatio();
	void genIntensityMean();
	double getMatchRatio(IonType stIonType);
	double getIntensityMean(IonType stIonType);
	virtual void mergePeaks();
};

}

#endif /* SCORES_H_ */
