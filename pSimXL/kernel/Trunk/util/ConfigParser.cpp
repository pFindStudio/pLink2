#include "../include/util.h"

using namespace std;

namespace sdk
{

ConfigParser::ConfigParser(){}

ConfigParser::ConfigParser(const char *szFilePath)
{
	read(szFilePath);
}

ConfigParser::~ConfigParser()
{
	close();
}

void ConfigParser::read(const char *szFilePath) throw(runtime_error)
{
	m_strFilePath = szFilePath;
	m_fin.open(m_strFilePath.c_str());
	if(!m_fin) {
		ErrorInfo err("ConfigParser", "read",
				"Failed to open configuration file: " + m_strFilePath);
		throw runtime_error(err.get());
	}
	parse();
}

bool ConfigParser::hasSection(const std::string &strSctn) const
{
	return (m_mpSctns.find(strSctn) != m_mpSctns.end());
}

bool ConfigParser::hasOption(const std::string &strSctn, const std::string &strOpt) const
{
	if(m_mpSctns.find(strSctn) != m_mpSctns.end()) {
		const map<string, string> &mpOpts = m_mpSctns.at(strSctn).mpOptVl;
		return (mpOpts.find(strOpt) != mpOpts.end());
	} else {
		return false;
	}
}

vector<string> ConfigParser::sections() const
{
	vector<string> vSctns;
	try {
		for(unsigned int i = 0; i < m_vSctns.size(); i++)
			vSctns.push_back(m_mpSctns.at(m_vSctns[i]).strSctn);
	} catch(exception &e) {
		ErrorInfo err("ConfigParser", "sections", "error section", e);
		throw runtime_error(err.get());
	} catch(...) {
		ErrorInfo err("ConfigParser", "sections", "unknown exception");
		throw runtime_error(err.get());
	}
	return vSctns;
}

vector<string> ConfigParser::options(const std::string &strSctn) const
{
	return m_mpSctns.at(strSctn).vOpts;
}

string ConfigParser::get(const string &strSctn, const string &strOpt) const
{
	return m_mpSctns.at(strSctn).mpOptVl.at(strOpt);
}

bool ConfigParser::getBool(const string &strSctn, const string &strOpt) const
{
	return getInt(strSctn, strOpt) ? true : false;
}

int ConfigParser::getInt(const string &strSctn, const string &strOpt) const
{
	return atoi(get(strSctn, strOpt).c_str());
}

double ConfigParser::getDouble(const string &strSctn, const string &strOpt) const
{
	return atof(get(strSctn, strOpt).c_str());
}

void ConfigParser::parse()
{
	vector<string>::iterator itr;

	for(first(); hasNext(); next()) {
		m_vSctns.push_back(m_stCurrent.strSctn);
		m_mpSctns.insert(pair<string, SectionItem>(m_stCurrent.strSctn, m_stCurrent));
	}
}

void ConfigParser::close()
{
	if(m_fin.is_open())
		m_fin.close();
}

void ConfigParser::first()
{
	string strLine;
	bool bFirst;
	size_t tPos;

	bFirst = true;
	while(getline(m_fin, strLine)) {
		strLine = trim(strLine);
		if(strLine.empty() || strLine[0] == ';' || strLine[0] == '#')
			continue;
		else if(strLine[0] == '[' && strLine[strLine.length()-1] == ']') {
			if(bFirst) {
				m_stCurrent.strSctn = strLine.substr(1, strLine.length()-2);
				bFirst = false;
			} else {
				m_stNext.strSctn = strLine.substr(1, strLine.length()-2);
				break;
			}
		} else {
			if((tPos = strLine.find_first_of('=', 0)) == string::npos) {
				ErrorInfo err("ConfigParser", "first", strLine);
				throw runtime_error(err.get());
			}
			else {
				string strOpt, strVl;
				strOpt = strLine.substr(0, tPos);
				strOpt = trim(strOpt);
				strVl = strLine.substr(tPos+1, strLine.length()-tPos-1);
				strVl = trim(strVl);

				m_stCurrent.vOpts.push_back(strOpt);
				m_stCurrent.mpOptVl.insert(pair<string, string>(strOpt, strVl));
			}
		}
	}
}

void ConfigParser::next()
{
	string strLine;
	size_t tPos;

	if(hasNext()) {
		m_stCurrent.strSctn = m_stNext.strSctn;
		m_stCurrent.vOpts.clear();
		m_stCurrent.mpOptVl.clear();
		m_stNext.strSctn = "";

		while(getline(m_fin, strLine)) {
			strLine = trim(strLine);
			if(strLine.empty() || strLine[0] == ';' || strLine[0] == '#')
				continue;
			else if(strLine[0] == '[' && strLine[strLine.length()-1] == ']') {
				m_stNext.strSctn = strLine.substr(1, strLine.length()-2);
				break;
			} else {
				if((tPos = strLine.find_first_of('=', 0)) == string::npos) {
					ErrorInfo err("ConfigParser", "next", strLine);
					throw runtime_error(err.get());
				} else {
					string strOpt, strVl;
					strOpt = strLine.substr(0, tPos);
					strOpt = trim(strOpt);
					strVl = strLine.substr(tPos+1, strLine.length()-tPos-1);
					strVl = trim(strVl);

					m_stCurrent.vOpts.push_back(strOpt);
					m_stCurrent.mpOptVl.insert(pair<string, string>(strOpt, strVl));
				}
			}
		}
	}
}

bool ConfigParser::hasNext()
{
	return !m_stCurrent.strSctn.empty();
}

void ConfigParser::testCase()
{
	try {
		// ConfigParser conf("E:\\Projects\\PeaksMatcher\\build\\GST_sequence_small.pindex");
		ConfigParser conf("F:\\Data\\UTP-B\\result\\1.sample\\search\\UTP_B_final1.pfind");
		vector<string> vSctns = conf.sections();
		for(vector<string>::iterator it = vSctns.begin(); it != vSctns.end(); it++) {
			cout<<(string)(*it).c_str()<<endl;
			vector<string> vOpts = conf.options((*it));
			for(unsigned int i = 0; i < vOpts.size(); i++)
				cout<<vOpts[i]<<"="<<conf.get((*it), vOpts[i])<<endl;
			cout<<endl;
		}
	} catch(exception &e) {
		ErrorInfo err("ConfigParser", "testCase", "", e);
		cerr<<err.get()<<endl;
	} catch(...) {
		ErrorInfo err("ConfigParser", "testCase", "caught an unknown exception");
		cerr<<err.get()<<endl;
	}
}

}
