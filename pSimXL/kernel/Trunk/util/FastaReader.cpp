#include "../include/index.h"

using namespace std;

namespace sdk
{

FastaReader::FastaReader(){}

FastaReader::FastaReader(const char *szFastaFilePath)
{
	m_strFastaFilePath = szFastaFilePath;
}

FastaReader::~FastaReader(){
	if(m_fin.is_open())
		endConnection();
}

void FastaReader::setFastaFilePath(const char *szFastaFilePath)
{
	m_strFastaFilePath = szFastaFilePath;
}

void FastaReader::getConnection() throw(runtime_error)
{
	m_fin.open(m_strFastaFilePath.c_str());
	if(!m_fin) {
		ErrorInfo err("FastaReader", "getConnection",
				"Failed to open file " + m_strFastaFilePath);
		throw runtime_error(err.get());
	}
}

void FastaReader::endConnection()
{
	if(m_fin.is_open())
		m_fin.close();
	m_fin.clear();
}

void FastaReader::first()
{
	string strLine;
	bool bFirst;

	bFirst = true;
	while(getline(m_fin, strLine)) {
		strLine = trim(strLine);
		if(strLine.length() > 0 && strLine[0] == '>') {
			size_t tPos = strLine.find_first_of(' ', 0);
			if(tPos == string::npos && strLine.length() > 0)
				tPos = strLine.length() - 1;
			if(bFirst) {
				m_current.strAc = strLine.substr(1, tPos);
				m_current.strDe = strLine.substr(tPos+1, strLine.length()-tPos-1);
				bFirst = false;
			} else {
				m_next.strAc = strLine.substr(1, tPos);
				m_next.strDe = strLine.substr(tPos+1, strLine.length()-tPos-1);
				break;
			}
		} else {
			m_current.strSq += strLine;
		}
	}
}

void FastaReader::next()
{
	string strLine;

	if(hasNext()) {
		m_current.strAc = m_next.strAc;
		m_current.strDe = m_next.strDe;
		m_current.strSq = "";
		m_next.strAc = "";
		m_next.strDe = "";

		while(getline(m_fin, strLine)) {
			strLine = trim(strLine);
			if(strLine.length() > 0 && strLine[0] == '>') {
				size_t tPos = strLine.find_first_of(' ', 0);
				if(tPos == string::npos && strLine.length() > 0)
					tPos = strLine.length() - 1;
				m_next.strAc = strLine.substr(1, tPos);
				m_next.strDe = strLine.substr(tPos+1, strLine.length()-tPos-1);
				break;
			} else {
				m_current.strSq += strLine;
			}
		}
	}
}

bool FastaReader::hasNext()
{
	return !m_current.strAc.empty();
}

ProteinItem FastaReader::currentItem() const
{
	return m_current;
}

void FastaReader::testCase()
{
	FastaReader fr("F:\\Data\\database\\1_38peptides.fasta");
	size_t tCnt = 0;

	fr.getConnection();
	for(fr.first(); fr.hasNext(); fr.next()) {
		ProteinItem stItem = fr.currentItem();
		cout<<stItem.strAc<<endl
				<<stItem.strDe<<endl
				<<stItem.strSq<<endl;
		tCnt++;
	}
	fr.endConnection();
	cout<<tCnt<<endl;
}

}
