#include "../include/util.h"

using namespace std;

namespace sdk {

const string NoSectionConfigParser::NO_SECTION = "nosection";

NoSectionConfigParser::NoSectionConfigParser() :
		ConfigParser()
{
}

NoSectionConfigParser::NoSectionConfigParser(const char *szFilePath) :
		ConfigParser(szFilePath)
{
}

NoSectionConfigParser::~NoSectionConfigParser()
{
}

void NoSectionConfigParser::parse()
{
	m_vSctns.push_back(NO_SECTION);

	SectionItem stItem;
	stItem.strSctn = NO_SECTION;
	string strLine;
	while(getline(m_fin, strLine)) {
		strLine = trim(strLine);
		if(strLine.empty() || strLine[0] == ';' || strLine[0] == '#') // ; and # deemed as comments
			continue;
		size_t tPos = 0;
		if((tPos = strLine.find_first_of('=', 0)) == string::npos) {
			ErrorInfo err("NoSectionConfigParser", "parse", strLine);
			throw runtime_error(err.get());
		} else {
			string strKey = strLine.substr(0, tPos);
			strKey = trim(strKey);
			string strValue = strLine.substr(tPos+1, strLine.length()-tPos-1);
			strValue = trim(strValue);
			stItem.mpOptVl.insert(pair<string,string>(strKey, strValue));
			stItem.vOpts.push_back(strKey);
		}
	}
	m_mpSctns.insert(pair<string, SectionItem>(NO_SECTION, stItem));
}

}
