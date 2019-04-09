#include "../include/util.h"

using namespace std;

namespace sdk
{
ElementDict *ElementDict::m_pInstance = NULL;

ElementDict *ElementDict::getInstance(const char *szConfigFilePath)
{
	if(m_pInstance == NULL)
		m_pInstance = new ElementDict(szConfigFilePath);
	return m_pInstance;
}

ElementDict *ElementDict::getInstance() throw(runtime_error)
{
	if(!m_pInstance) {
		ErrorInfo err("ElementDict", "getInstance",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_pInstance;
}

void ElementDict::destroyInstance()
{
	if(m_pInstance) {
		delete m_pInstance;
		m_pInstance = NULL;
	}
}

ElementDict::~ElementDict()
{
}

const Element &ElementDict::getElement(const std::string &strElementName) const
{
	if(!m_pInstance) {
		ErrorInfo err("ElementDict", "getElement",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpElementMap.at(strElementName);
}

double ElementDict::getMonoMass(const std::string &strElementName) const
{
	if(!m_pInstance) {
		ErrorInfo err("ElementDict", "getMonoMass",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpElementMap.at(strElementName).vIsotopes[0].lfMass;
}

double ElementDict::getMonoAbundance(const std::string &strElementName) const
{
	if(!m_pInstance) {
		ErrorInfo err("ElementDict", "getMonoAbundance",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpElementMap.at(strElementName).vIsotopes[0].lfAbundance;
}

Isotope ElementDict::getMonoIsotope(const std::string &strElementName) const
{
	if(!m_pInstance) {
		ErrorInfo err("ElementDict", "getMonoIsotope",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpElementMap.at(strElementName).vIsotopes[0];
}

double ElementDict::getMass(const std::string &strElementName, size_t tIdx) const
{
	if(!m_pInstance) {
		ErrorInfo err("ElementDict", "getMass",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpElementMap.at(strElementName).vIsotopes[tIdx].lfMass;
}

double ElementDict::getAbundance(const std::string &strElementName, size_t tIdx) const
{
	if(!m_pInstance) {
		ErrorInfo err("ElementDict", "getAbundance",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpElementMap.at(strElementName).vIsotopes[tIdx].lfAbundance;
}

Isotope ElementDict::getIsotope(const std::string &strElementName, size_t tIdx) const
{
	if(!m_pInstance) {
		ErrorInfo err("ElementDict", "getIsotope",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpElementMap.at(strElementName).vIsotopes[tIdx];
}

double ElementDict::calculateMolecularMass(const std::string &strFormula)
{
	double lfResult = 0.0;

	if(m_pInstance) {
		string strCopy(strFormula);
		vector<string> vParts;
		split(strCopy, ')', vParts);
		for(size_t n = 0; n < vParts.size(); ++n) {
			vector<string> vElement;
			split(vParts[n], '(', vElement);
			lfResult += getMonoMass(vElement[0]) * atoi(vElement[1].c_str());
		}
	}
	return lfResult;
}

ElementDict::ElementDict(const char *szConfigFilePath)
{
	m_strConfigFilePath = szConfigFilePath;
	parse();
}

void ElementDict::parse()
{
	try {
		NoSectionConfigParser stConf;
		stConf.read(m_strConfigFilePath.c_str());
		int nNum = stConf.getInt(NoSectionConfigParser::NO_SECTION, "@NUMBER_ELEMENT");
		for(int n = 1; n <= nNum; ++n) {
			ostringstream oss;
			oss<<"E"<<n;
			string strVl = stConf.get(NoSectionConfigParser::NO_SECTION, oss.str());
			vector<string> vParts;
			split(strVl, '|', vParts);
			if(vParts.size() < 3) {
				throw strVl;
			}

			Element stElement;
			stElement.strName = vParts[0];
			vector<string> vMasses;
			split(vParts[1], ',', vMasses);
			vector<string> vAbundances;
			split(vParts[2], ',', vAbundances);
			if(vMasses.size() != vAbundances.size()) {
				throw strVl;
			}

			for(size_t tIdx = 0; tIdx < vMasses.size(); ++tIdx) {
				Isotope stIsotope;
				stIsotope.lfMass = atof(vMasses[tIdx].c_str());
				stIsotope.lfAbundance = atof(vAbundances[tIdx].c_str());
				stElement.vIsotopes.push_back(stIsotope);
			}
			m_mpElementMap.insert(pair<string, Element>(stElement.strName, stElement));
		}
	} catch(exception &e) {
		ErrorInfo err("ElementDict", "parse", "format error.", e);
		throw runtime_error(err.get());
	} catch(...) {
		ErrorInfo err("ElementDict", "parse", "caught an unknown exception.");
		throw runtime_error(err.get());
	}
}

}
