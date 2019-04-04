#include "../include/util.h"

using namespace std;

namespace sdk
{


EnzymeDict *EnzymeDict::m_pInstance = NULL;

EnzymeDict *EnzymeDict::getInstance(const char *szConfigFilePath)
{
	if(m_pInstance == NULL)
		m_pInstance = new EnzymeDict(szConfigFilePath);
	return m_pInstance;
}

EnzymeDict *EnzymeDict::getInstance() throw(runtime_error)
{
	if(m_pInstance == NULL) {
		ErrorInfo err("EnzymeDict", "getInstance",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_pInstance;
}

void EnzymeDict::destroyInstance()
{
	if(m_pInstance) {
		delete m_pInstance;
		m_pInstance = NULL;
	}
}

EnzymeDict::~EnzymeDict()
{
}

const Enzyme &EnzymeDict::getEnzyme(const std::string &strEnzymeName) const
{
	if(!m_pInstance) {
		ErrorInfo err("EnzymeDict", "getEnzyme",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpEnzymeMap.at(strEnzymeName);
}

std::string EnzymeDict::getNCleaveSite(const string &strEnzymeName) const
{
	if(!m_pInstance) {
		ErrorInfo err("EnzymeDict", "getNCleaveSite",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpEnzymeMap.at(strEnzymeName).strNCleaveSite;
}

std::string EnzymeDict::getNExceptSite(const string &strEnzymeName) const
{
	if(!m_pInstance) {
		ErrorInfo err("EnzymeDict", "getNExceptSite",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpEnzymeMap.at(strEnzymeName).strNExceptSite;
}

std::string EnzymeDict::getCCleaveSite(const string &strEnzymeName) const
{
	if(!m_pInstance) {
		ErrorInfo err("EnzymeDict", "getCCleaveSite",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpEnzymeMap.at(strEnzymeName).strCCleaveSite;
}

std::string EnzymeDict::getCExceptSite(const string &strEnzymeName) const
{
	if(!m_pInstance) {
		ErrorInfo err("EnzymeDict", "getCExceptSite",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpEnzymeMap.at(strEnzymeName).strCExceptSite;
}

EnzymeType EnzymeDict::getEnzymeType(const std::string &strEnzymeName) const
{
	if(!m_pInstance) {
		ErrorInfo err("EnzymeDict", "getEnzymeType",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpEnzymeMap.at(strEnzymeName).eType;
}

EnzymeDict::EnzymeDict(const char *szConfigFilePath)
{
	m_strConfigFilePath = szConfigFilePath;
	parse();
}

void EnzymeDict::parse()
{
	ConfigParser conf(m_strConfigFilePath.c_str());
	vector<string> vOptions = conf.options("enzyme");
	ostringstream oss;
	istringstream iss;

	for(size_t i = 0; i < vOptions.size(); ++i) {
		Enzyme stEnzyme;
		int nType;

		stEnzyme.strName = vOptions[i];
		iss.str(conf.get("enzyme", vOptions[i]));
		iss>>stEnzyme.strNCleaveSite>>stEnzyme.strNExceptSite
			>>stEnzyme.strCCleaveSite>>stEnzyme.strCExceptSite
			>>nType;

		if(!stEnzyme.strNCleaveSite.compare("_")) {
			stEnzyme.strNCleaveSite.clear();
		}

		if(!stEnzyme.strNExceptSite.compare("_")) {
			stEnzyme.strNExceptSite.clear();
		}

		if(!stEnzyme.strCCleaveSite.compare("_")) {
			stEnzyme.strCCleaveSite.clear();
		}

		if(!stEnzyme.strCExceptSite.compare("_")) {
			stEnzyme.strCExceptSite.clear();
		}

		if(nType == 0)
			stEnzyme.eType = ET_SPECIFIC;
		else if(nType == 1)
			stEnzyme.eType = ET_SEMI;
		else if(nType == 2)
			stEnzyme.eType = ET_NONE;
		else {
			ErrorInfo err("EnzymeDict", "parse", "unknown enzyme type, please check enzyme.ini.");
			throw runtime_error(err.get());
		}

		m_mpEnzymeMap.insert(pair<string, Enzyme>(vOptions[i], stEnzyme));
		iss.clear();
	}
}

XLinkerDict *XLinkerDict::m_pInstance = NULL;

XLinkerDict *XLinkerDict::getInstance(const char *szConfigFilePath)
{
	if(m_pInstance == NULL)
		m_pInstance = new XLinkerDict(szConfigFilePath);
	return m_pInstance;
}

XLinkerDict *XLinkerDict::getInstance() throw(runtime_error)
{
	if(!m_pInstance) {
		ErrorInfo err("XLinkerDict", "getInstance",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_pInstance;
}

XLinkerDict::~XLinkerDict()
{
}

void XLinkerDict::destroyInstance()
{
	if(m_pInstance)
		delete m_pInstance;
	m_pInstance = NULL;
}

const XLinker &XLinkerDict::getXLinker(const std::string &strLinkerName) const
{
	if(!m_pInstance) {
		ErrorInfo err("XLinkerDict", "getXLinker",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpLinkerMap.at(strLinkerName);
}

const std::string &XLinkerDict::getAlphaAA(const std::string &strLinkerName) const
{
	if(!m_pInstance) {
		ErrorInfo err("XLinkerDict", "getAlphaAA",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpLinkerMap.at(strLinkerName).strAlphaAA;
}

const std::string &XLinkerDict::getBetaAA(const std::string &strLinkerName) const
{
	if(!m_pInstance) {
		ErrorInfo err("XLinkerDict", "getBetaAA",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpLinkerMap.at(strLinkerName).strBetaAA;
}

double XLinkerDict::getMonoDiff(const std::string &strLinkerName) const
{
	if(!m_pInstance) {
		ErrorInfo err("XLinkerDict", "getMonoDiff",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpLinkerMap.at(strLinkerName).lfMonoMassDiff;
}

double XLinkerDict::getAvrgDiff(const std::string &strLinkerName) const
{
	if(!m_pInstance) {
		ErrorInfo err("XLinkerDict", "getAvrgDiff",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpLinkerMap.at(strLinkerName).lfAvrgMassDiff;
}

double XLinkerDict::getMLMonoDiff(const std::string &strLinkerName) const
{
	if(!m_pInstance) {
		ErrorInfo err("XLinkerDict", "getMLMonoDiff",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpLinkerMap.at(strLinkerName).lfMLMonoMassDiff;
}

double XLinkerDict::getMLAvrgDiff(const std::string &strLinkerName) const
{
	if(!m_pInstance) {
		ErrorInfo err("XLinkerDict", "getMLAvrgDiff",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpLinkerMap.at(strLinkerName).lfMLAvrgMassDiff;
}

XLinkerDict::XLinkerDict(const char *szConfigFilePath)
{
	m_strConfigFilePath = szConfigFilePath;
	parse();
}

void XLinkerDict::parse()
{
	ConfigParser conf(m_strConfigFilePath.c_str());
	size_t tNum = conf.getInt("xlink", "total");
	ostringstream oss;
	istringstream iss;

	for(unsigned int i = 0; i < tNum; i++) {
		string strOpt, strVl;
		XLinker stLinker;

		oss.str("");
		oss.clear();
		oss<<"name"<<(i+1);
		strOpt = conf.get("xlink", oss.str());
		strVl = conf.get("xlink", strOpt);
		iss.str(strVl);
		iss>>stLinker.strAlphaAA>>stLinker.strBetaAA
			>>stLinker.lfMonoMassDiff>>stLinker.lfAvrgMassDiff
			>>stLinker.lfMLMonoMassDiff>>stLinker.lfMLAvrgMassDiff;
		m_mpLinkerMap.insert(pair<string, XLinker>(strOpt, stLinker));
		iss.clear();
	}
}

InstrumentDict *InstrumentDict::m_pInstance = NULL;

InstrumentDict *InstrumentDict::getInstance(const char *szConfigFilePath)
{
	if(m_pInstance == NULL)
		m_pInstance = new InstrumentDict(szConfigFilePath);
	return m_pInstance;
}

InstrumentDict *InstrumentDict::getInstance() throw(runtime_error)
{
	if(!m_pInstance) {
		ErrorInfo err("InstrumentDict", "getInstance",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_pInstance;
}

InstrumentDict::~InstrumentDict()
{
}

void InstrumentDict::destroyInstance()
{
	if(m_pInstance) {
		delete m_pInstance;
		m_pInstance = NULL;
	}
}

const Instrument &InstrumentDict::getInstrument(const std::string &strInstrumentName) const
{
	if(!m_pInstance) {
		ErrorInfo err("InstrumentDict", "getInstrument",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_pInstance->m_mpInstrumentMap.at(strInstrumentName);
}

const std::vector<IonType> &InstrumentDict::getIonTypes(const std::string &strInstrumentName) const
{
	if(!m_pInstance) {
		ErrorInfo err("InstrumentDict", "getFragTolType",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_pInstance->m_mpInstrumentMap.at(strInstrumentName).vIonTypes;
}

InstrumentDict::InstrumentDict(const char *szConfigFilePath)
{
	m_strConfigFilePath = szConfigFilePath;
	parse();
}

void InstrumentDict::parse()
{
	ConfigParser conf(m_strConfigFilePath.c_str());
	size_t tNum = conf.getInt("instrument", "total");
	ostringstream oss;
	istringstream iss;

	for(size_t i = 0; i < tNum; i++) {
		string strSctn, strVl;
		Instrument stInst;

		oss.str("");
		oss.clear();
		oss<<"name"<<(i+1);
		strSctn = conf.get("instrument", oss.str());
		size_t tTotal = conf.getInt(strSctn, "ion_type_total");
		for(size_t j = 0; j < tTotal; ++j) {
			IonType stIonType;

			oss.str("");
			oss.clear();
			oss<<"ion_type"<<(j+1);
			strVl = conf.get(strSctn, oss.str());
			iss.str(strVl);
			stIonType.lfLoss = 0.0;

			iss>>stIonType.cType;
			switch(stIonType.cType) {
			case 'a':
				stIonType.lfLoss += CO_MONO_MASS;
				break;
			case 'b':
				break;
			case 'c':
				stIonType.lfLoss -= NH3_MONO_MASS;
				break;
			case 'x':
				stIonType.lfLoss -= CO_MONO_MASS - 2 * H_MONO_MASS;
				break;
			case 'y':
				break;
			case 'z':
				stIonType.lfLoss += NH2_MONO_MASS;
				break;
			case 'm': case 'p': case 'n': case 'u':
				// m -- internal ions
				// p -- precurosr ions (one peptide)
				// n -- KL: ya-NH3 at the linker site
				// u -- yP
				break;
			case 'v':
				// v -- bP
				break;
			default:
				ErrorInfo err("InstrumentDict", "parse",
						"no such iontype, please check instrument.ini.");
				throw runtime_error(err.get());
			}

			iss>>stIonType.nCharge;

			int nValue;
			iss>>nValue;
			stIonType.lfLoss += NH3_MONO_MASS * nValue;

			iss>>nValue;
			stIonType.lfLoss += H2O_MONO_MASS * nValue;

			double lfValue;
			iss>>lfValue;
			stIonType.lfLoss += lfValue;

			iss>>nValue;
			stIonType.bContinuity = bool(nValue);

			iss>>stIonType.nNCoexistence>>stIonType.nCCoexistence>>stIonType.nHasLinker;
			stInst.vIonTypes.push_back(stIonType);
			iss.clear();
		}
		m_mpInstrumentMap.insert(pair<string, Instrument>(strSctn, stInst));
		iss.clear();
	}
}

}
