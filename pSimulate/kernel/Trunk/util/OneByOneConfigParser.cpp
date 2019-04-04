
#include "../include/util.h"

using namespace std;

namespace sdk {

OneByOneConfigParser::OneByOneConfigParser() {

}

OneByOneConfigParser::OneByOneConfigParser(const char *szFilePath) :
		ConfigParser(szFilePath) {
}

OneByOneConfigParser::~OneByOneConfigParser() {
}

bool OneByOneConfigParser::hasOption(const std::string &strOpt) const {
	return ConfigParser::hasOption(section(), strOpt);
}

std::vector<std::string> OneByOneConfigParser::options() const {
	return ConfigParser::options(section());
}

const std::string &OneByOneConfigParser::section() const {
	if (m_vSctns.size() > 0)
		return m_vSctns[0];
	else
		return strEmpty;
}

std::string OneByOneConfigParser::get(const std::string &strOpt) const {
	return ConfigParser::get(section(), strOpt);
}

bool OneByOneConfigParser::getBool(const std::string &strOpt) const {
	return ConfigParser::getBool(section(), strOpt);
}

int OneByOneConfigParser::getInt(const std::string &strOpt) const {
	return ConfigParser::getInt(section(), strOpt);
}

double OneByOneConfigParser::getDouble(const std::string &strOpt) const {
	return ConfigParser::getDouble(section(), strOpt);
}

void OneByOneConfigParser::first() {
	ConfigParser::first();
	if (ConfigParser::hasNext()) {
		m_vSctns.push_back(m_stCurrent.strSctn);
		m_mpSctns.insert(
				pair<string, SectionItem>(m_stCurrent.strSctn, m_stCurrent));
	}
}

void OneByOneConfigParser::next() {
	ConfigParser::next();
	if (ConfigParser::hasNext()) {
		m_vSctns.clear();
		m_mpSctns.clear();
		m_vSctns.push_back(m_stCurrent.strSctn);
		m_mpSctns.insert(
				pair<string, SectionItem>(m_stCurrent.strSctn, m_stCurrent));
	}
}

bool OneByOneConfigParser::hasNext() {
	return ConfigParser::hasNext();
}

void OneByOneConfigParser::parse() {
	// override, this method will be invoke by read() method of class ConfigParser
}

}

