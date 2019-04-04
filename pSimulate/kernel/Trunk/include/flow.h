#ifndef FLOW_H_
#define FLOW_H_
#include "util.h"
#include "bias.h"
#include "sdk.h"
#include "decorators.h"
#include "io.h"
#include "index.h"
#include "scorers.h"

namespace sdk
{

#define EXPERIMENT

class Flow {
protected:
	const SearchParameter *m_pParameter;
	std::vector<Spectrum> m_vSpectra;

public:
	Flow();
	virtual ~Flow();

	virtual void init() = 0;
	virtual void run() = 0;
	virtual void close() = 0;

protected:
    virtual bool isModLinkSiteConflict(const Peptide &stPep, size_t tIdx);
    virtual void getPrecursorBorder(PrecursorWindow &stWindow, size_t tSpecIndex);

};

class FlowFactory {
    const SearchParameter *m_pParameter;

public:
    FlowFactory();
    virtual ~FlowFactory();

    Flow *getFlow(FlowType eType) const;
};

class BenchMarkFlow : public Flow {
	IonIndexState *m_pIonState;
	Trace *m_pTrace;
	IonIndexer *m_pIndexer;
	RefinedScorer *m_pRefinedScorer;
	XLinkDecorator *m_pXLinkDecorator;
	std::vector<PrecursorWindow> m_vSpecInfo;
	XLinkPeptideResult m_stCurXLinkPep;
	WindowsHash m_stWndHash;
	MS2Output *m_pMS2Output;
	SimplePSMIO *m_pSimplePSMIO;
	Bias *m_pBias;
	std::vector<size_t> m_vSpectraNum;
	std::vector<size_t> m_vCurSpecNo;
	size_t m_tRound;


public:
	BenchMarkFlow();
	virtual ~BenchMarkFlow();

	virtual void init();
	virtual void run();
	virtual void close();

private:
	void gotoRefineStage(size_t tSpecId);
	void initMatchResults();
	void searchPeptide();
	void searchCommonPeptide(std::vector<Peptide*> &vPeps);
	void searchMonoPeptide(std::vector<Peptide*> &vPeps);
	void searchLoopPeptide(std::vector<Peptide*> &vPeps);
	void searchXLinkPeptide(std::vector<Peptide*> &vPeps);

	void startWrite(const std::string &strTmpFilePath);
	void writeSpectrum(size_t tSpecId);
	void writePSM(XLinkPeptideResult &stXLinkPep, size_t tSpecId);
	size_t generateBigRand();
};

class BenchMarkPerfectFlow : public Flow {
	IonIndexState *m_pIonState;
	Trace *m_pTrace;
	IonIndexer *m_pIndexer;
	RefinedScorer *m_pRefinedScorer;
	XLinkDecorator *m_pXLinkDecorator;
	std::vector<PrecursorWindow> m_vSpecInfo;
	XLinkPeptideResult m_stCurXLinkPep;
	WindowsHash m_stWndHash;
	MS2Output *m_pMS2Output;
	SimplePSMIO *m_pSimplePSMIO;
	Bias *m_pBias;
	std::vector<size_t> m_vSpectraNum;
	std::vector<size_t> m_vCurSpecNo;
	size_t m_tRound;

public:
	BenchMarkPerfectFlow();
	virtual ~BenchMarkPerfectFlow();

	virtual void init();
	virtual void run();
	virtual void close();

private:
	void gotoRefineStage(size_t tSpecId);
	void initMatchResults();
	void searchPeptide();
	void searchCommonPeptide(std::vector<Peptide*> &vPeps);
	void searchMonoPeptide(std::vector<Peptide*> &vPeps);
	void searchLoopPeptide(std::vector<Peptide*> &vPeps);
	void searchXLinkPeptide(std::vector<Peptide*> &vPeps);

	void startWrite(const std::string &strTmpFilePath);
	void writeSpectrum(size_t tSpecId);
	void writePSM(XLinkPeptideResult &stXLinkPep, size_t tSpecId);
	size_t generateBigRand();
};

class BenchMarkSimpleFlow : public Flow {
	IonIndexState *m_pIonState;
	Trace *m_pTrace;
	IonIndexer *m_pIndexer;
	RefinedScorer *m_pRefinedScorer;
	XLinkDecorator *m_pXLinkDecorator;
	std::vector<PrecursorWindow> m_vSpecInfo;
	XLinkPeptideResult m_stCurXLinkPep;
	WindowsHash m_stWndHash;
	MS2Output *m_pMS2Output;
	SimplePSMIO *m_pSimplePSMIO;
	Bias *m_pBias;
	std::vector<size_t> m_vSpectraNum;
	std::vector<size_t> m_vCurSpecNo;
	size_t m_tRound;


public:
	BenchMarkSimpleFlow();
	virtual ~BenchMarkSimpleFlow();

	virtual void init();
	virtual void run();
	virtual void close();

private:
	void gotoRefineStage(size_t tSpecId);
	void initMatchResults();
	void searchPeptide();
	void searchCommonPeptide(std::vector<Peptide*> &vPeps);
	void searchMonoPeptide(std::vector<Peptide*> &vPeps);
	void searchLoopPeptide(std::vector<Peptide*> &vPeps);
	void searchXLinkPeptide(std::vector<Peptide*> &vPeps);

	void startWrite(const std::string &strTmpFilePath);
	void writeSpectrum(size_t tSpecId);
	void writePSM(XLinkPeptideResult &stXLinkPep, size_t tSpecId);
	size_t generateBigRand();
};

} /* namespace sdk */

#endif /* FLOW_H_ */
