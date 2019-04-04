#include "../include/util.h"

using namespace std;

namespace sdk
{

PrecursorWindow::PrecursorWindow() :
		m_tIndex(0), m_lfMinMass(0.0), m_lfMaxMass(0.0)
{
}

PrecursorWindow::PrecursorWindow(size_t tIndex, double lfMinMass, double lfMaxMass) :
		m_tIndex(tIndex), m_lfMinMass(lfMinMass), m_lfMaxMass(lfMaxMass)
{
}

PeptideTol::PeptideTol() :
		m_lfPeptideTol(20.0), m_lfPeptideTolBase(0.0),
		m_ePeptideTolType(TT_PPM), m_ePeptideTolBaseType(TT_DA)
{
}

DoubleItem::DoubleItem() :
		m_tOrder(0), m_lfNumber(0.0)
{
}

DoubleItem::DoubleItem(size_t tOrder, double lfNumber) :
		m_tOrder(tOrder), m_lfNumber(lfNumber)
{
}

bool operator<(const DoubleItem &stItem1, const DoubleItem &stItem2)
{
	return stItem1.m_lfNumber < stItem2.m_lfNumber;
}

bool operator<=(const DoubleItem &stItem1, const DoubleItem &stItem2)
{
	return stItem1.m_lfNumber <= stItem2.m_lfNumber;
}

bool operator>(const DoubleItem &stItem1, const DoubleItem &stItem2)
{
	return stItem1.m_lfNumber > stItem2.m_lfNumber;
}

bool operator>=(const DoubleItem &stItem1, const DoubleItem &stItem2)
{
	return stItem1.m_lfNumber >= stItem2.m_lfNumber;
}

bool operator==(const DoubleItem &stItem1, const DoubleItem &stItem2)
{
	return stItem1.m_lfNumber == stItem2.m_lfNumber;
}

bool operator!=(const DoubleItem &stItem1, const DoubleItem &stItem2)
{
	return stItem1.m_lfNumber != stItem2.m_lfNumber;
}

SearchParameter::SearchParameter() :
		m_strElementFilePath("element.ini"),
		m_strAAListFilePath("aa.ini"),
		m_strEnzymeListFilePath("enzyme.ini"),
		m_strModificationFilePath("modification.ini"), m_strLinkerFilePath("xlink.ini"),
		m_strInstrumentFilePath("instrument.ini"), m_strDBName("database"),
		m_strFastaFilePath(""), m_strEnzymeName("Trypsin"),
		m_tMaxMissSite(2), m_tMinPepLen(4), m_tMaxPepLen(60), m_tMinFragLen(2),
		m_tMaxFragLen(60), m_lfMinPepMass(400), m_lfMaxPepMass(6000),
		m_lfMinFragMass(2*57.021464), m_lfMaxFragMass(6000), m_bAutoReverse(0),
		m_tMaxMemSize(500*1024*1024),
		m_tMaxModsNo(MAX_MODIFY_NUM),
		m_strRefinedInstrument("HCD"), m_lfFragTol(20.0), m_eFragTolType(TT_PPM), m_lfPepTol(20.0), m_ePepTolType(TT_PPM),
		m_eRefinedScoreType(ST_XLINK_SIMPLE),
		m_eLogRank(LRT_INFO), m_eFlowType(FT_SIMULATION),
		m_tProcessorNo(1),
		m_tProteinsBatchSize(1500),
		m_strDecoySign("REVERSE_"),
		m_strSpecTitle("example"), m_eSpecType(MFT_MGF), m_nLostPeakPercentage(0)
{

	char szBuffer[STR_BUF_SIZE];
	getcwd(szBuffer, STR_BUF_SIZE);
	m_strWorkDir = szBuffer;
	if(m_strWorkDir[m_strWorkDir.length()-1] != SLASH) {
		m_strWorkDir += SLASH;
	}

	m_tSpectraNo.resize(4, 0);

}

SearchParameter::~SearchParameter()
{
}

ParameterReader *ParameterReader::m_pInstance = NULL;

ParameterReader *ParameterReader::getInstance(const string &strFilePath)
{
	if(m_pInstance == NULL)
		m_pInstance = new ParameterReader(strFilePath);
	return m_pInstance;
}

ParameterReader *ParameterReader::getInstance()
{
	if(!m_pInstance) {
		ErrorInfo err("ParameterReader", "getInstance",
				"no instance been created, try to invoke getInstance(const string&)");
		throw runtime_error(err.get());
	}
	return m_pInstance;
}

void ParameterReader::destroyInstance()
{
	if(m_pInstance) {
		delete m_pInstance;
		m_pInstance = NULL;
	}
}

ParameterReader::~ParameterReader()
{
	if(m_pParameter) {
		delete m_pParameter;
		m_pParameter = NULL;
	}
}

void ParameterReader::loadParameter()
{
	try {
		ConfigParser conf(m_strFilePath.c_str());

		// config
		if(conf.hasSection("config") && conf.hasOption("config", "element_list"))
			m_pParameter->m_strElementFilePath = conf.get("config", "element_list");

		if(conf.hasSection("config") && conf.hasOption("config", "aa_list"))
			m_pParameter->m_strAAListFilePath = conf.get("config", "aa_list");

		if(conf.hasSection("config") && conf.hasOption("config", "enzyme_list"))
			m_pParameter->m_strEnzymeListFilePath = conf.get("config", "enzyme_list");

		if(conf.hasSection("config") && conf.hasOption("config", "modification_list"))
			m_pParameter->m_strModificationFilePath = conf.get("config", "modification_list");

		if(conf.hasSection("config") && conf.hasOption("config", "xlink_list"))
			m_pParameter->m_strLinkerFilePath = conf.get("config", "xlink_list");

		if(conf.hasSection("config") && conf.hasOption("config", "instrument_list"))
			m_pParameter->m_strInstrumentFilePath = conf.get("config", "instrument_list");

		// database
		if(conf.hasOption("database", "db_name"))
			m_pParameter->m_strDBName = conf.get("database", "db_name");

		if(conf.hasOption("database", "db_path"))
			m_pParameter->m_strFastaFilePath = conf.get("database", "db_path");

		if(conf.hasOption("database", "enzyme_name"))
			m_pParameter->m_strEnzymeName = conf.get("database", "enzyme_name");

		if(conf.hasOption("database", "max_miss_site"))
			m_pParameter->m_tMaxMissSite = conf.getInt("database", "max_miss_site");

		if(conf.hasOption("database", "min_pep_len"))
			m_pParameter->m_tMinPepLen = conf.getInt("database", "min_pep_len");

		if(conf.hasOption("database", "max_pep_len"))
			m_pParameter->m_tMaxPepLen = conf.getInt("database", "max_pep_len");


		if(conf.hasOption("database", "min_pep_mass"))
			m_pParameter->m_lfMinPepMass = conf.getDouble("database", "min_pep_mass");

		if(conf.hasOption("database", "max_pep_mass"))
			m_pParameter->m_lfMaxPepMass = conf.getDouble("database", "max_pep_mass");



		if(conf.hasOption("database", "auto_reverse"))
			m_pParameter->m_bAutoReverse = conf.getBool("database", "auto_reverse");


		// modification
		if(conf.hasOption("modification", "fix_total")) {
			size_t tSize = conf.getInt("modification", "fix_total");
			for(size_t i = 0; i < tSize; ++i) {
				ostringstream oss;
				oss<<"fix_mod"<<(i+1);
				string strMod = conf.get("modification", oss.str());
				m_pParameter->m_vFixedMods.push_back(strMod);
			}
		}

		if(conf.hasOption("modification", "var_total")) {
			size_t tSize = conf.getInt("modification", "var_total");
			for(size_t i = 0; i < tSize; ++i) {
				ostringstream oss;
				oss<<"var_mod"<<(i+1);
				string strMod = conf.get("modification", oss.str());
				m_pParameter->m_vVariableMods.push_back(strMod);
			}
		}

		if(conf.hasOption("modification", "max_number")) {
			m_pParameter->m_tMaxModsNo = conf.getInt("modification", "max_number");
		}

		// linker
		if(conf.hasOption("linker", "linker_total")) {
			size_t tSize = conf.getInt("linker", "linker_total");
			for(size_t i = 0; i < tSize; ++i) {
				ostringstream oss;
				oss<<"linker"<<(i+1);
				string strLinker = conf.get("linker", oss.str());
				m_pParameter->m_vLinkers.push_back(strLinker);
			}
		}

		// ions
		if(conf.hasOption("ions", "refined_instrument")) {
			m_pParameter->m_strRefinedInstrument = conf.get("ions", "refined_instrument");
		}

		if(conf.hasOption("ions", "fragment_tol")) {
			m_pParameter->m_lfFragTol = conf.getDouble("ions", "fragment_tol");
		}

		if(conf.hasOption("ions", "fragment_tol_type")) {
			string strType = conf.get("ions", "fragment_tol_type");
			toupper(strType);
			if(!strType.compare("PPM")) {
				m_pParameter->m_eFragTolType = TT_PPM;
			} else if(!strType.compare("DA")) {
				m_pParameter->m_eFragTolType = TT_DA;
			} else {
				ErrorInfo err("ParameterReader", "loadParameter", "unkown fragment tolerance type!");
				throw runtime_error(err.get());
			}
		}

		if(conf.hasOption("ions", "peptide_tol")) {
			m_pParameter->m_lfPepTol = conf.getDouble("ions", "peptide_tol");
		}

		if(conf.hasOption("ions", "peptide_tol_type")) {
			string strType = conf.get("ions", "peptide_tol_type");
			toupper(strType);
			if(!strType.compare("PPM")) {
				m_pParameter->m_ePepTolType = TT_PPM;
			} else if(!strType.compare("DA")) {
				m_pParameter->m_ePepTolType = TT_DA;
			} else {
				ErrorInfo err("ParameterReader", "loadParameter", "unkown peptide tolerance type!");
				throw runtime_error(err.get());
			}
		}


		// score
		if(conf.hasOption("score", "refined_score_type")) {
			string strType = conf.get("score", "refined_score_type");
			if(!strType.compare("ST_XLINK_OLD")) {
				m_pParameter->m_eRefinedScoreType = ST_XLINK_LEGACY;
			} else if(!strType.compare("ST_XLINK_PERFECT")) {
				m_pParameter->m_eRefinedScoreType = ST_XLINK_PERFECT;
			} else if(!strType.compare("ST_XLINK_SIMPLE")) {
				m_pParameter->m_eRefinedScoreType = ST_XLINK_SIMPLE;
			} else {
				m_pParameter->m_eRefinedScoreType = ST_XLINK_SIMPLE;
			}
		}


		// flow
		if(conf.hasOption("flow", "flow_type")) {
			string strFlowType = conf.get("flow", "flow_type");
			if(!strFlowType.compare("FT_SIMULATION")) {
				m_pParameter->m_eFlowType = FT_SIMULATION;
			}
			else {
				ErrorInfo err("ParameterReader", "loadParameter", "unknown flow type!");
				throw runtime_error(err.get());
			}
		}



		if(conf.hasOption("flow", "spectra_regular_num")) {
			m_pParameter->m_tSpectraNo[PT_COMMON] = conf.getInt("flow", "spectra_regular_num");
		}

		if(conf.hasOption("flow", "spectra_mono_num")) {
			m_pParameter->m_tSpectraNo[PT_MONO] = conf.getInt("flow", "spectra_mono_num");
		}

		if(conf.hasOption("flow", "spectra_loop_num")) {
			m_pParameter->m_tSpectraNo[PT_LOOP] = conf.getInt("flow", "spectra_loop_num");
		}

		if(conf.hasOption("flow", "spectra_xlink_num")) {
			m_pParameter->m_tSpectraNo[PT_XLINK] = conf.getInt("flow", "spectra_xlink_num");
		}

		if(conf.hasOption("flow", "processor_num")) {
			m_pParameter->m_tProcessorNo = conf.getInt("flow", "processor_num");
			if(m_pParameter->m_tProcessorNo > getCPUCoreNo()) {
				m_pParameter->m_tProcessorNo = getCPUCoreNo();
			}
		}



		if(conf.hasOption("flow", "proteins_batch_size")) {
			m_pParameter->m_tProteinsBatchSize = conf.getInt("flow", "proteins_batch_size");
		}



		if(conf.hasOption("flow", "decoy_sign")) {
			m_pParameter->m_strDecoySign = conf.get("flow", "decoy_sign");
		}

		if(conf.hasOption("flow", "result_output_path")) {
			m_pParameter->m_strOutputDirPath = conf.get("flow", "result_output_path");
			makeDir(m_pParameter->m_strOutputDirPath);
		}

		if(conf.hasOption("flow", "lost_peak_percentage")) {
			m_pParameter->m_nLostPeakPercentage = conf.getInt("flow", "lost_peak_percentage");
		}

		// spectrum
		if(conf.hasOption("spectrum", "spec_title")) {
			m_pParameter->m_strSpecTitle = conf.get("spectrum", "spec_title");
		}

		if(conf.hasOption("spectrum", "spec_type")) {
			string strType = conf.get("spectrum", "spec_type");
			toupper(strType);
			if(!strType.compare("MGF")) {
				m_pParameter->m_eSpecType = MFT_MGF;
			} else if(!strType.compare("MS2")) {
				m_pParameter->m_eSpecType = MFT_MS2;
			} else if(!strType.compare("PF")) {
				m_pParameter->m_eSpecType = MFT_PF;
			} else {
				throw runtime_error("load spectrum type failed!");
			}
		}

		if(conf.hasOption("spectrum", "spec_num")) {
			size_t tFileNum = conf.getInt("spectrum", "spec_num");
			for(size_t tFileId = 0; tFileId < tFileNum; ++tFileId) {
				ostringstream oss;
				oss<<"spec_path"<<(tFileId+1);
				string strFilepath = conf.get("spectrum", oss.str());
				m_pParameter->m_vSpecFilePath.push_back(strFilepath);
			}
		}
	} catch(exception &e) {
		ErrorInfo err("ParameterReader", "loadParameter", "load parameter failed!", e);
		throw runtime_error(err.get());
	} catch(...) {
		ErrorInfo err("ParameterReader", "loadParameter", "caught an unknown exception.");
		throw runtime_error(err.get());
	}
}

SearchParameter *ParameterReader::getParameter()
{
	if(!m_pInstance) {
		ErrorInfo err("ParameterReader", "getParameter",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}

	if(!m_pParameter) {
		ErrorInfo err("ParameterReader", "getParameter", "get parameter failed!");
		throw runtime_error(err.get());
	}

	return m_pParameter;
}

ParameterReader::ParameterReader(const string &strPath) :
		m_strFilePath(strPath), m_pParameter(NULL)
{
	m_pParameter = new SearchParameter;
}


IonIndexState::IonIndexState() :
		m_tStartClock(clock()), m_tStartTime(0),
		m_tTotalProteins(0), m_tCurrentProteins(0)
{
	time(&m_tStartTime);
}

IonIndexState::~IonIndexState()
{

}

string IonIndexState::getProgress()
{
	ostringstream oss;
	string strBlanks("       ");

	oss<<"Start Time: "<<ctime(&m_tStartTime)
	   <<strBlanks<<"Time cost: "<<(clock() - m_tStartClock) * 1.0 / CLOCKS_PER_SEC<<"s\n"
	   <<strBlanks<<"Proteins: "<<m_tCurrentProteins<<" / "<<m_tTotalProteins<<"\n";

	return oss.str();
}

void IonIndexState::setClock()
{
	m_tStartClock = clock();
	time(&m_tStartTime);
	m_tTotalProteins = m_tCurrentProteins = 0;
}

void IonIndexState::setTotalProteins(size_t tTotal)
{
	m_tTotalProteins = tTotal;
	m_tCurrentProteins = 0;
}

void IonIndexState::setCurrentProteins(size_t tCurrent)
{
	m_tCurrentProteins += tCurrent;
}

bool IonIndexState::isFinished() const
{
	return m_tCurrentProteins >= m_tTotalProteins;
}

clock_t IonIndexState::getStartClock() const
{
	return m_tStartClock;
}

time_t IonIndexState::getStartTime() const
{
	return m_tStartTime;
}


}



