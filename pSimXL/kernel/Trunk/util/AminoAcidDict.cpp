#include "../include/util.h"

using namespace std;

namespace sdk {

AminoAcidDict *AminoAcidDict::m_pInstance = NULL;

AminoAcidDict *AminoAcidDict::getInstance(const char *szConfigFilePath) {
	if (m_pInstance == NULL)
		m_pInstance = new AminoAcidDict(szConfigFilePath);
	return m_pInstance;
}

AminoAcidDict *AminoAcidDict::getInstance() throw (runtime_error) {
	if (!m_pInstance) {
		ErrorInfo err("AminoAcidDict", "getInstance",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_pInstance;
}

void AminoAcidDict::destroyInstance() {
	if (m_pInstance) {
		delete m_pInstance;
		m_pInstance = NULL;
	}
}

AminoAcidDict::~AminoAcidDict() {
}

const AminoAcid &AminoAcidDict::getAminoAcid(char cAA) const {
	if (!m_pInstance) {
		ErrorInfo err("AminoAcidDict", "getAminoAcid",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpAAMap.at(cAA);
}

double AminoAcidDict::getMonoMass(char cAA) const {
	if (!m_pInstance) {
		ErrorInfo err("AminoAcidDict", "getMonoMass",
				"no instance been created, try to invoke getInstance(const char *)");
		throw runtime_error(err.get());
	}
	return m_mpAAMap.at(cAA).lfMono;
}

double AminoAcidDict::calculatePeptideMass(const std::string &strSq)
{
	double lfResult = 0.0;

	if(m_pInstance) {
		for(size_t n = 0; n < strSq.size(); ++n) {
			lfResult += getMonoMass(strSq[n]);
		}

		lfResult += 2 * ElementDict::getInstance()->getMonoMass("H");
		lfResult += ElementDict::getInstance()->getMonoMass("O");
	}
	return lfResult;
}

AminoAcidDict::AminoAcidDict(const char *szConfigFilePath) {
	m_strConfigFilePath = szConfigFilePath;
	parse();
}

void AminoAcidDict::parse() {
	try {
		NoSectionConfigParser stConf;
		stConf.read(m_strConfigFilePath.c_str());
		int nNum = stConf.getInt(NoSectionConfigParser::NO_SECTION, "@NUMBER_RESIDUE");
		for(int n = 1; n <= nNum; ++n) {
			ostringstream oss;
			oss<<"R"<<n;
			string strVl = stConf.get(NoSectionConfigParser::NO_SECTION, oss.str());
			vector<string> vParts;
			split(strVl, '|', vParts);
			if(vParts.size() < 2) {
				throw strVl;
			}

			AminoAcid stAA;
			stAA.cAa = vParts[0][0];
			stAA.strFormula = vParts[1];
			stAA.lfMono = ElementDict::getInstance()->calculateMolecularMass(stAA.strFormula);
			m_mpAAMap.insert(pair<char, AminoAcid>(stAA.cAa, stAA));
		}
	} catch(exception &e) {
		ErrorInfo err("AminoAcidDict", "parse", "format error.", e);
		throw runtime_error(err.get());
	} catch(...) {
		ErrorInfo err("AminoAcidDict", "parse", "caught an unknown exception.");
		throw runtime_error(err.get());
	}
}
}
