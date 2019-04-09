
#include "../include/util.h"
#include "../include/index.h"

using namespace std;
namespace sdk {

EnvironmentInitializer *EnvironmentInitializer::m_pInstance = NULL;

void EnvironmentInitializer::init(const char *szParamPath) {
	if (m_pInstance == NULL) {
		string strParamPath = string(szParamPath);
		m_pInstance = new EnvironmentInitializer(strParamPath);
	}
}

void EnvironmentInitializer::destroy() {
	if (m_pInstance != NULL) {
		delete m_pInstance;
		m_pInstance = NULL;
	}
}

SearchParameter *EnvironmentInitializer::getSearchParameter() {
	if (m_pInstance == NULL)
		return NULL;
	else
		return m_pInstance->m_pParameter;
}


EnvironmentInitializer::EnvironmentInitializer(std::string &strParamPath) :
		m_pParameter(NULL),  m_pTrace(NULL) {
	m_pTrace = Trace::getInstance();
	m_pTrace->info("Welcome to use pSimulate v%s.", PSIMULATE_VERSION.c_str());
	m_pTrace->info("pSimulate initializing...");

	m_pTrace->debug("Reading parameters from file " + strParamPath);
	ParameterReader *pReader = ParameterReader::getInstance(strParamPath);
	pReader->loadParameter();
	m_pParameter = pReader->getParameter();

	createDirs();

	// redirect Trace to output running information
	string strLogs[] = {m_pParameter->m_strOutputDirPath, "tmps", "log.out"};
	vector<string> vLogs(strLogs, strLogs+sizeof(strLogs)/sizeof(string));
	string strLogFilePath;
	join(strLogFilePath, SLASH, vLogs);
	m_pTrace = Trace::getInstance(strLogFilePath.c_str(), LAT_CONSOLE_AND_FILE);
	m_pTrace->changeLogRank(m_pParameter->m_eLogRank);

	// notice: do not change the initializion order, since dependency exists
	m_pTrace->debug("Loading element settings...");
	ElementDict::getInstance(m_pParameter->m_strElementFilePath.c_str());

	m_pTrace->debug("Loading amino acids settings...");
	AminoAcidDict::getInstance(m_pParameter->m_strAAListFilePath.c_str());

	m_pTrace->debug("Loading enzyme settings...");
	EnzymeDict::getInstance(m_pParameter->m_strEnzymeListFilePath.c_str());

	m_pTrace->debug("Loading modification settings...");
	ModificationDict::getInstance(
			m_pParameter->m_strModificationFilePath.c_str());

	m_pTrace->debug("Loading linker settings...");
	XLinkerDict::getInstance(m_pParameter->m_strLinkerFilePath.c_str());

	m_pTrace->debug("Loading instrument settings...");
	InstrumentDict::getInstance(m_pParameter->m_strInstrumentFilePath.c_str());

	m_pTrace->debug("Loading database settings...");
	if(m_pParameter->m_bAutoReverse) {
		m_pTrace->info("Generating reverse database...");
		reverseDatabase();
	}

	m_pTrace->debug("Initialize constants...");
	ConstantStore::getInstance();

	m_pTrace->info("pSimulate is ready to simulate.\n");
}

void EnvironmentInitializer::createDirs()
{
	string strDir;

	string strTmp[] = {m_pParameter->m_strOutputDirPath, "tmps"};
	vector<string> vTmps(strTmp, strTmp+sizeof(strTmp)/sizeof(string));
	join(strDir, SLASH, vTmps);
	makeDir(strDir);

	string strReport[] = {m_pParameter->m_strOutputDirPath, "results"};
	vector<string> vReport(strReport, strReport+sizeof(strReport)/sizeof(string));
	join(strDir, SLASH, vReport);
	makeDir(strDir);

}

void EnvironmentInitializer::reverseDatabase()
{
	const size_t SUFFIX_LEN = strlen(".fasta");
	size_t tLen = m_pParameter->m_strFastaFilePath.length();
	if(tLen > SUFFIX_LEN && m_pParameter->m_strFastaFilePath.substr(tLen-SUFFIX_LEN, SUFFIX_LEN).compare(".fasta")) {
		throw "Error fasta file path.";
	}
	string strPath = m_pParameter->m_strFastaFilePath.substr(0, tLen-SUFFIX_LEN) + "_comb_rever.fasta";
	ofstream outfile(strPath.c_str());
	if(!outfile) {
		throw "Failed to create file: " + strPath;
	}

	FastaReader *pFr = new FastaReader(m_pParameter->m_strFastaFilePath.c_str());
	pFr->getConnection();
	for(pFr->first(); pFr->hasNext(); pFr->next()) {
		ProteinItem stItem = pFr->currentItem();
		// forward sequence
		outfile<<">"<<stItem.strAc<<" "<<stItem.strDe<<endl;
		outfile<<stItem.strSq<<endl;

		// backward sequence
		outfile<<">"<<m_pParameter->m_strDecoySign<<stItem.strAc<<" "<<stItem.strDe<<endl;
		reverse(stItem.strSq.begin(), stItem.strSq.end());
		outfile<<stItem.strSq<<endl;
	}
	pFr->endConnection();
	delete pFr;
	outfile.close();
	outfile.clear();
	m_pParameter->m_strFastaFilePath = strPath;
}

EnvironmentInitializer::~EnvironmentInitializer() {

	ParameterReader::destroyInstance();
	Trace::destroyInstance();
}

} /* namespace proteomics_sdk */
