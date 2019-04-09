#include "../include/util.h"

using namespace std;

namespace sdk {

ModificationDict *ModificationDict::m_pInstance = NULL;

ModificationDict *ModificationDict::getInstance(const char *szConfigFilePath)
{
	if(m_pInstance == NULL)
		m_pInstance = new ModificationDict(szConfigFilePath);
	return m_pInstance;
}

ModificationDict *ModificationDict::getInstance() throw(runtime_error)
{
	if(!m_pInstance) {
		ErrorInfo err("ModificationDict", "getInstance",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_pInstance;
}

void ModificationDict::destroyInstance()
{
	if(m_pInstance)
		delete m_pInstance;
	m_pInstance = NULL;
}

ModificationDict::~ModificationDict()
{
}

const Modification &ModificationDict::getModification(const std::string &strModName) const
{
	if(!m_pInstance) {
		ErrorInfo err("ModificationDict", "getModification",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpModMap.at(strModName);
}

const std::string &ModificationDict::getAminoAcid(const std::string &strModName) const
{
	if(!m_pInstance) {
		ErrorInfo err("ModificationDict", "getAminoAcid",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpModMap.at(strModName).strAA;
}

ModificationType ModificationDict::getModificationType(const std::string &strModName) const
{
	if(!m_pInstance) {
		ErrorInfo err("ModificationDict", "getModificationType",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpModMap.at(strModName).eModType;
}

double ModificationDict::getMonoDiff(const std::string &strModName) const
{
	if(!m_pInstance) {
		ErrorInfo err("ModificationDict", "getMonoDiff",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpModMap.at(strModName).lfMonoDiff;
}

double ModificationDict::getAvrgDiff(const std::string &strModName) const
{
	if(!m_pInstance) {
		ErrorInfo err("ModificationDict", "getAvrgDiff",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpModMap.at(strModName).lfAvrgDiff;
}

std::vector<double> ModificationDict::getMonoNeutralLossDiff(const std::string &strModName) const
{
	if(!m_pInstance) {
		ErrorInfo err("ModificationDict", "getMonoNeutralLossDiff",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpModMap.at(strModName).vMonoNeutralLossDiff;
}

std::vector<double> ModificationDict::getAvrgNeutralLossDiff(const std::string &strModName) const
{
	if(!m_pInstance) {
		ErrorInfo err("ModificationDict", "getAvrgNeutralLossDiff",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpModMap.at(strModName).vAvrgNeutralLossDiff;
}

ModificationDict::ModificationDict(const char *szConfigFilePath)
{
	m_strConfigFilePath = szConfigFilePath;
	parse();
}

void ModificationDict::parse()
{
	NoSectionConfigParser stConf;
	stConf.read(m_strConfigFilePath.c_str());
	size_t tNum = stConf.getInt(NoSectionConfigParser::NO_SECTION, "@NUMBER_MODIFICATION");
	ostringstream oss;
	istringstream iss;

	for(unsigned int i = 0; i < tNum; i++) {
		string strOpt, strVl, strType, strNLNum;
		Modification stMod;

		// parse: name#=modification_name regular_or_not
		oss.str("");
		oss.clear();
		oss<<"name"<<(i+1);
		strOpt = stConf.get(NoSectionConfigParser::NO_SECTION, oss.str());
		iss.str(strOpt);
		iss>>stMod.strName;

		strVl = stConf.get(NoSectionConfigParser::NO_SECTION, stMod.strName);
		iss.str(strVl);
		iss>>stMod.strAA>>strType>>stMod.lfMonoDiff>>stMod.lfAvrgDiff>>strNLNum;
		if(!stMod.strAA.compare("[ALL]"))
			stMod.strAA = ALL_AA;

		// parse: modification_name=site type mono_mass average_mass
		// neutral_loss_number neutral_loss_mono_mass neutral_loss_average_mass formula
		if(!strType.compare("NORMAL")) {
			stMod.eModType = MT_PEP_NORMAL;
		} else if(!strType.compare("PEP_N")) {
			stMod.eModType = MT_PEP_N;
		} else if(!strType.compare("PEP_C")) {
			stMod.eModType = MT_PEP_C;
		} else if(!strType.compare("PRO_N")) {
			stMod.eModType = MT_PRO_N;
		} else if(!strType.compare("PRO_C")) {
			stMod.eModType = MT_PRO_C;
		} else {
			ErrorInfo err("ModificationDict", "parse",
					"error modification type");
			throw runtime_error(err.get());
		}

		int tNLNum = atoi(strNLNum.c_str());
		if(tNLNum > 0) {
			stMod.vMonoNeutralLossDiff.resize(tNLNum);
			stMod.vAvrgNeutralLossDiff.resize(tNLNum);
		}
		for(int j = 0; j < tNLNum; j++)
			iss>>stMod.vMonoNeutralLossDiff[j]>>stMod.vAvrgNeutralLossDiff[j];
		m_mpModMap.insert(pair<string, Modification>(stMod.strName, stMod));
		iss.clear();
	}
}

}
