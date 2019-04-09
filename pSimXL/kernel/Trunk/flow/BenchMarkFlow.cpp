#include "../include/flow.h"

using namespace std;

namespace sdk
{

BenchMarkFlow::BenchMarkFlow() :
		Flow(),   m_pIonState(NULL),
		m_pTrace(NULL), m_pIndexer(NULL),
		m_pRefinedScorer(NULL), m_pXLinkDecorator(NULL), m_pMS2Output(0), m_pSimplePSMIO(0),
		m_pBias(0), m_vSpectraNum(m_pParameter->m_tSpectraNo), m_tRound(3)
{
	for(size_t i = 0; i < 4; i++)
	{
		m_vCurSpecNo.push_back(0);
	}
}

BenchMarkFlow::~BenchMarkFlow()
{
	close();
}

void BenchMarkFlow::init()
{
	try {

		if(m_pIonState == NULL) {
			m_pIonState = new IonIndexState();
		}

		if(m_pTrace == NULL) {
			m_pTrace = Trace::getInstance();
		}


		if(m_pRefinedScorer == NULL) {

			RefinedScorerFactory stRefinedScorerFactory(m_pParameter);
			m_pRefinedScorer = stRefinedScorerFactory.getScorer(m_pParameter->m_eRefinedScoreType);

			//m_pRefinedScorer = new BMIonGenerator(m_pParameter);
		}

		if(m_pXLinkDecorator == NULL) {
			m_pXLinkDecorator = new XLinkDecorator(m_pParameter);
		}

		if(m_pMS2Output == NULL) {
			m_pMS2Output = new MGFOutput;
		}

		if(m_pSimplePSMIO == NULL) {
			m_pSimplePSMIO = new SimplePSMIO;
		}

		if(m_pBias == NULL) {
			m_pBias = new Bias;
		}

	} catch(exception &e) {
		ErrorInfo err("BenchMarkFlow", "init", "initialize BenchMark flow failed!", e);
		throw runtime_error(err.get());
	} catch(...) {
		ErrorInfo err("BenchMarkFlow", "init", "caught an unknown exception.");
		throw runtime_error(err.get());
	}
}

void BenchMarkFlow::searchCommonPeptide(vector<Peptide*> &vPeps)
{
	PeptideResult stOne;
	m_stCurXLinkPep.clear();
	m_stCurXLinkPep.m_ePepType = PT_COMMON;
	if(m_vCurSpecNo[m_stCurXLinkPep.m_ePepType] >= m_vSpectraNum[m_stCurXLinkPep.m_ePepType])
	{
		return;
	}
	stOne.m_lfOpenMass = 0.0;
	stOne.m_stLinkMod.m_nSite = -1;
	stOne.m_stLinkMod.m_nId = -1;
	for(size_t idx = 0; idx < vPeps.size(); ++idx) {
		size_t i = generateBigRand() % vPeps.size();
		stOne.m_stPep = *vPeps[i];
		double lfMass = vPeps[i]->m_lfMass;
		m_stCurXLinkPep.m_lfCalcMH = lfMass + PROTON_MASS + H2O_MONO_MASS ;
		m_stCurXLinkPep.m_stAlphaPep = stOne;
		// refinement score
		gotoRefineStage(m_vSpecInfo[0].m_tIndex);
		if(m_vCurSpecNo[m_stCurXLinkPep.m_ePepType] >= m_vSpectraNum[m_stCurXLinkPep.m_ePepType])
		{
			return;
		}
	}
}

void BenchMarkFlow::searchMonoPeptide(vector<Peptide*> &vPeps)
{
	PeptideResult stOne;
	m_stCurXLinkPep.clear();
	m_stCurXLinkPep.m_ePepType = PT_MONO;
	if(m_vCurSpecNo[m_stCurXLinkPep.m_ePepType] >= m_vSpectraNum[m_stCurXLinkPep.m_ePepType])
	{
		return;
	}
	for(size_t tLinkerId = 0; tLinkerId < m_pParameter->m_vLinkers.size(); ++tLinkerId){
		double lfMonoLinkerMass = m_pXLinkDecorator->getLinkerMass(true, tLinkerId);
		stOne.m_stLinkMod.m_nId = tLinkerId;
		stOne.m_lfOpenMass = lfMonoLinkerMass;
		for(size_t idx = 0; idx < vPeps.size(); ++idx) {
			size_t i = generateBigRand() % vPeps.size();
			stOne.m_stPep = *vPeps[i];
			double lfMass = vPeps[i]->m_lfMass + lfMonoLinkerMass;

			m_stCurXLinkPep.m_lfCalcMH = lfMass + PROTON_MASS + H2O_MONO_MASS ;
			for(size_t tSite1 = 0; tSite1 < vPeps[i]->m_strSq.length(); ++tSite1) {
				if(m_pXLinkDecorator->isXLinkSiteLegal(*vPeps[i], tSite1, tLinkerId) &&
						!isModLinkSiteConflict(*vPeps[i], tSite1)){
					stOne.m_stLinkMod.m_nSite = tSite1;
					m_stCurXLinkPep.m_stAlphaPep = stOne;
					// refinement score
					gotoRefineStage(m_vSpecInfo[0].m_tIndex);
					if(m_vCurSpecNo[m_stCurXLinkPep.m_ePepType] >= m_vSpectraNum[m_stCurXLinkPep.m_ePepType])
					{
						return;
					}
				}
			}
		}
	}
}

void BenchMarkFlow::searchLoopPeptide(vector<Peptide*> &vPeps)
{
	PeptideResult stOne;
	PeptideResult stTheOther;
	m_stCurXLinkPep.clear();
	m_stCurXLinkPep.m_ePepType = PT_LOOP;
	if(m_vCurSpecNo[m_stCurXLinkPep.m_ePepType] >= m_vSpectraNum[m_stCurXLinkPep.m_ePepType])
	{
		return;
	}
	vector<int> vXLinkSite;
	vXLinkSite.reserve(m_pParameter->m_tMaxPepLen + 1);
	for(size_t tLinkerId = 0; tLinkerId < m_pParameter->m_vLinkers.size(); ++tLinkerId){
		double lfInterMass = m_pXLinkDecorator->getLinkerMass(false,tLinkerId);
		stOne.m_stLinkMod.m_nId = tLinkerId;
		stOne.m_lfOpenMass = lfInterMass;
		for(size_t idx = 0; idx < vPeps.size(); ++idx) {
			size_t i = generateBigRand() % vPeps.size();
			stOne.m_stPep = *vPeps[i];
			double lfMass = vPeps[i]->m_lfMass + lfInterMass;
			m_stCurXLinkPep.m_lfCalcMH = lfMass + PROTON_MASS + H2O_MONO_MASS ;
			vXLinkSite.clear();
			for(size_t tSite = 0; tSite < vPeps[i]->m_strSq.length(); ++tSite){
				if(m_pXLinkDecorator->isXLinkSiteLegal(*vPeps[i], tSite, tLinkerId) &&
						!isModLinkSiteConflict(*vPeps[i], tSite)){
					vXLinkSite.push_back(tSite);
				}
			}
			for(size_t tSiteIdx1 = 0; tSiteIdx1 < vXLinkSite.size(); ++tSiteIdx1) {
				size_t tSite1 = vXLinkSite[tSiteIdx1];
				for(size_t tSiteIdx2 = tSiteIdx1 + 1; tSiteIdx2 < vXLinkSite.size(); ++tSiteIdx2) {
					size_t tSite2 = vXLinkSite[tSiteIdx2];
					if(m_pXLinkDecorator->isXLinkSiteLegal(*vPeps[i], tSite1, *vPeps[i], tSite2, tLinkerId)){
						stOne.m_stLinkMod.m_nSite = tSite1;
						m_stCurXLinkPep.m_stAlphaPep = stOne;
						m_stCurXLinkPep.m_stBetaPep.m_stLinkMod.m_nId = stOne.m_stLinkMod.m_nId;
						m_stCurXLinkPep.m_stBetaPep.m_stLinkMod.m_nSite = tSite2;

						// refinement score
						gotoRefineStage(m_vSpecInfo[0].m_tIndex);
						if(m_vCurSpecNo[m_stCurXLinkPep.m_ePepType] >= m_vSpectraNum[m_stCurXLinkPep.m_ePepType])
						{
							return;
						}
					}
				}
			}
		}
	}
}

void BenchMarkFlow::searchXLinkPeptide(vector<Peptide*> &vPeps)
{
	PeptideResult stOne;
	PeptideResult stTheOther;
	m_stCurXLinkPep.clear();
	m_stCurXLinkPep.m_ePepType = PT_XLINK;
	if(m_vCurSpecNo[m_stCurXLinkPep.m_ePepType] >= m_vSpectraNum[m_stCurXLinkPep.m_ePepType])
	{
		return;
	}
	for(size_t tLinkerId = 0; tLinkerId < m_pParameter->m_vLinkers.size(); ++tLinkerId){
		double lfInterMass = m_pXLinkDecorator->getLinkerMass(false,tLinkerId);
		stOne.m_stLinkMod.m_nId = tLinkerId;
		stTheOther.m_stLinkMod.m_nId = tLinkerId;
		for(size_t idx1 = 0; idx1 < vPeps.size(); ++idx1) {
			size_t i = generateBigRand() % vPeps.size();
			stOne.m_stPep = *vPeps[i];
			for(size_t idx2 = idx1; idx2 < vPeps.size(); ++idx2) {
				if(generateBigRand() % 30 == 0)
					break;
				size_t j = generateBigRand() % vPeps.size();
				double lfPairMass = vPeps[i]->m_lfMass + vPeps[j]->m_lfMass + H2O_MONO_MASS + lfInterMass;
				stTheOther.m_stPep = *vPeps[j];
				stTheOther.m_lfOpenMass = stOne.m_stPep.m_lfMass + H2O_MONO_MASS + lfInterMass;
				stOne.m_lfOpenMass = stTheOther.m_stPep.m_lfMass + H2O_MONO_MASS + lfInterMass;
				m_stCurXLinkPep.m_lfCalcMH = lfPairMass + PROTON_MASS + H2O_MONO_MASS;
				for(size_t tSite1 = 0; tSite1 < vPeps[i]->m_strSq.length(); ++tSite1) {
					if(m_pXLinkDecorator->isXLinkSiteLegal(*vPeps[i], tSite1, tLinkerId) &&
							!isModLinkSiteConflict(*vPeps[i], tSite1)){
						stOne.m_stLinkMod.m_nSite = tSite1;
						for(size_t tSite2 = 0; tSite2 < vPeps[j]->m_strSq.length(); ++tSite2){
							if(m_pXLinkDecorator->isXLinkSiteLegal(*vPeps[j], tSite2,tLinkerId ) &&
								!isModLinkSiteConflict(*vPeps[j], tSite2) &&
								m_pXLinkDecorator->isXLinkSiteLegal(*vPeps[i], tSite1, *vPeps[j], tSite2, tLinkerId)) {
								stTheOther.m_stLinkMod.m_nSite = tSite2;
								if(vPeps[i]->m_lfMass >= vPeps[j]->m_lfMass) {
									m_stCurXLinkPep.m_stAlphaPep = stOne;
									m_stCurXLinkPep.m_stBetaPep = stTheOther;
								} else {
									m_stCurXLinkPep.m_stAlphaPep = stTheOther;
									m_stCurXLinkPep.m_stBetaPep = stOne;
								}
								// refinement score
								gotoRefineStage(m_vSpecInfo[0].m_tIndex);
								if(m_vCurSpecNo[m_stCurXLinkPep.m_ePepType] >= m_vSpectraNum[m_stCurXLinkPep.m_ePepType])
								{
									return;
								}
							}
						}
					}
				}
			}
		}
	}
}

void BenchMarkFlow::searchPeptide()
{
	try {
		/* search peptide using EnumPeptide indexer */
		// initialize ion-indexer
		if(m_pIndexer == NULL) {
			m_pIndexer = new EnumPeptideIndexer(m_pParameter);
		}

		m_pIonState->setClock();
		m_pIonState->setTotalProteins(m_pIndexer->getProteinsNum());
		while(m_pIndexer->hasNextLoad()) {
			m_pIonState->setCurrentProteins(m_pIndexer->loadNext());
			m_pTrace->info(m_pIonState->getProgress());
			std::vector<Peptide *> vPeps;
			size_t tMaxRound = m_tRound * 500;
			dynamic_cast<EnumPeptideIndexer *>(m_pIndexer)->getXLinkPeptides(vPeps);
			if(vPeps.size() != 0) {
				for(size_t i = 0; i < tMaxRound; i++) {
					if(m_vCurSpecNo[PT_XLINK] < m_vSpectraNum[PT_XLINK])
					{
						searchXLinkPeptide(vPeps);
					}

					if(m_vCurSpecNo[PT_MONO] < m_vSpectraNum[PT_MONO])
					{
						searchMonoPeptide(vPeps);
					}

					if(m_vCurSpecNo[PT_LOOP] < m_vSpectraNum[PT_LOOP])
					{
						searchLoopPeptide(vPeps);
					}
				}
			}
			vPeps.clear();
			dynamic_cast<EnumPeptideIndexer *>(m_pIndexer)->getPeptides(vPeps);
			if(vPeps.size() != 0) {
				for(size_t i = 0; i < tMaxRound; i++) {
					if(m_vCurSpecNo[PT_COMMON] < m_vSpectraNum[PT_COMMON])
					{
						searchCommonPeptide(vPeps);
					}

				}
			}
		}

		vector<PrecursorWindow>().swap(m_vSpecInfo);

		if(m_pIndexer != NULL) {
			delete m_pIndexer;
			m_pIndexer = NULL;
		}

	} catch(exception &e) {
		ErrorInfo err("EnumFlow", "searchPeptide", "search peptide failed!", e);
		throw runtime_error(err.get());
	} catch(...) {
		ErrorInfo err("EnumFlow", "searchPeptide", "caught an unknown exception");
		throw runtime_error(err.get());
	}
}

void BenchMarkFlow::run()
{
	try {

		srand(time(NULL)); // 放在这里有效，放在init无效

		// initialize match results
		initMatchResults();

		ostringstream oss;
		oss<<m_pParameter->m_strOutputDirPath<<SLASH<<"results"<<SLASH;

		string strTmpFilePath = oss.str();

		startWrite(strTmpFilePath);

		// search the peptide
		m_pTrace->info("Generating spectra...");
		searchPeptide();

		m_pMS2Output->endWrite();
		m_pSimplePSMIO->endWrite();

		m_pTrace->info("pSimulate completed.");

	} catch(exception &e) {
		ErrorInfo err("BenchMarkFlow", "run", "run enum flow failed!", e);
		throw runtime_error(err.get());
	} catch(...) {
		ErrorInfo err("BenchMarkFlow", "run", "caught an unknown exception.");
		throw runtime_error(err.get());
	}
}

void BenchMarkFlow::close()
{


	if(m_pIonState) {
		delete m_pIonState;
		m_pIonState = NULL;
	}

	if(m_pIndexer) {
		delete m_pIndexer;
		m_pIndexer = NULL;
	}

	if(m_pRefinedScorer) {
		delete m_pRefinedScorer;
		m_pRefinedScorer = NULL;
	}

	if(m_pXLinkDecorator) {
		delete m_pXLinkDecorator;
		m_pXLinkDecorator = NULL;
	}

	if(m_pMS2Output) {
		delete m_pMS2Output;
	}

	if(m_pSimplePSMIO) {
		delete m_pSimplePSMIO;
	}

	if(m_pBias){
		delete m_pBias;
	}

}

void BenchMarkFlow::initMatchResults()
{
	try {
		m_vSpecInfo.clear();

		m_vSpecInfo.reserve(m_vSpectra.size());

		for(size_t i = 0;i < m_vSpectra.size(); ++i) {

			PrecursorWindow window;
			window.m_tIndex = i;
			getPrecursorBorder(window, i);
			m_vSpecInfo.push_back(window);


//			// compute each precursor window
//			vector<PrecursorWindow> vWindow(m_pParameter->m_vPeptideTols.size(), PrecursorWindow());
//			for(size_t j = 0; j < vWindow.size(); ++j) {
//				vWindow[j].m_tIndex = i;
//				getPrecursorBorder(vWindow[j], j, i);
//				m_vSpecInfo.push_back(vWindow[j]);
//			}

		}

		sort(m_vSpecInfo.begin(), m_vSpecInfo.end());
		m_stWndHash.construct(m_vSpecInfo);
	} catch(exception &e) {
		ErrorInfo err("BenchMarkFlow", "initMatchResults", "initialize match results failed!", e);
		throw runtime_error(err.get());
	} catch(...) {
		ErrorInfo err("BenchMarkFlow", "initMatchResults", "caught an unknown exception");
		throw runtime_error(err.get());
	}
}



void BenchMarkFlow::startWrite(const std::string &strTmpFilePath)
{
	try {
		string strFilePath = getFilePath(strTmpFilePath);
		string strFileName1("Simulation.mgf");
		int nTotal = 0;
		m_pMS2Output->startWrite(strFilePath + strFileName1, nTotal);
		string strFileName2("Simulation.csv");
		m_pSimplePSMIO->startWrite(strFilePath + strFileName2);


	} catch(exception &e) {
		ErrorInfo err("BenchMarkFlow", "startWrite", "startWrite failed!", e);
		throw runtime_error(err.get());
	} catch(...) {
		ErrorInfo err("BenchMarkFlow", "startWrite", "caught an unknown exception.");
		throw runtime_error(err.get());
	}
}

void BenchMarkFlow::writeSpectrum(size_t tSpecId)
{
	int idx = 0;
	m_pMS2Output->writeNext(m_vSpectra[tSpecId], idx);
}

void BenchMarkFlow::writePSM(XLinkPeptideResult &stXLinkPep, size_t tSpecId)
{
	m_pSimplePSMIO->writeNext(m_vSpectra[tSpecId], stXLinkPep);
}

void BenchMarkFlow::gotoRefineStage(size_t tSpecId)
{
	try {
		if(m_vCurSpecNo[m_stCurXLinkPep.m_ePepType] >= m_vSpectraNum[m_stCurXLinkPep.m_ePepType])
		{
			return;
		}

		size_t tModNum1 = m_stCurXLinkPep.m_stAlphaPep.m_stPep.m_vMods.size();
		size_t tModNum2 = m_stCurXLinkPep.m_stBetaPep.m_stPep.m_vMods.size();
		while(tModNum1--) {
			if(rand() % 2 == 0)
				return;
		}
		while(tModNum2--) {
			if(rand() % 2 == 0)
				return;
		}
		if(rand() % m_tRound != 0) {
			return;
		}
		size_t tLen = m_stCurXLinkPep.m_stAlphaPep.m_stPep.m_strSq.size();
		if(tLen < 7 && rand()%3 != 0) {
			return;
		}
		else if(tLen > 20 && rand()%(tLen*tLen) >= 40) {
			return;
		}

		if(m_stCurXLinkPep.m_ePepType == PT_XLINK) {
			int nTempRand = rand() % 4;
			if(1 != nTempRand) {
				return;
			}
		}
		m_pRefinedScorer->setPeptide(&m_stCurXLinkPep);
		m_pRefinedScorer->setSpectrum(&m_vSpectra[tSpecId]);

//		size_t tSamePepNum = rand()%3 + 1;
//		for(size_t i = 0; i < tSamePepNum; i++)
//		{
			m_stCurXLinkPep.m_lfScore = m_pRefinedScorer->score();
			if(m_vSpectra[tSpecId].m_vPeaks.size() > 3) {
				writeSpectrum(tSpecId);
				writePSM(m_stCurXLinkPep, tSpecId);
				m_vCurSpecNo[m_stCurXLinkPep.m_ePepType] ++;
				if(m_vCurSpecNo[m_stCurXLinkPep.m_ePepType] >= m_vSpectraNum[m_stCurXLinkPep.m_ePepType])
				{
					return;
				}
			}
//		}

	} catch(exception &e) {
		ErrorInfo err("BenchMarkFlow", "gotoRefineStage", "refinement scoring failed!", e);
		throw runtime_error(err.get());
	} catch(...) {
		ErrorInfo err("BenchMarkFlow", "gotoRefineStage", "caught an unknown exception.");
		throw runtime_error(err.get());
	}
}

size_t BenchMarkFlow::generateBigRand()
{
	 return (
	     (((size_t)rand()<<24)&0xFF000000u)
	    |(((size_t)rand()<<12)&0x00FFF000u)
	    |(((size_t)rand()    )&0x00000FFFu));
}

}



