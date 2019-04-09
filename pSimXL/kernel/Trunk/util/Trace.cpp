#include "../include/util.h"

using namespace std;

namespace sdk
{

void WindowsConsoleAppender::append(const char *szInfo, AppenderColorType eColor)
{
	HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);

	switch(eColor) {
	case ACT_RED:
		SetConsoleTextAttribute(hConsole, FOREGROUND_RED|FOREGROUND_INTENSITY);
		break;
	case ACT_GREEN:
		SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN|FOREGROUND_INTENSITY);
		break;
	default:
		/* default set on RGB, white */
		break;
	}
	fprintf(m_fp, "%s\n", szInfo);
	fflush(m_fp);
	SetConsoleTextAttribute(hConsole, FOREGROUND_RED|FOREGROUND_GREEN|FOREGROUND_BLUE);
}

void WindowsConsoleAppender::append(const string &strInfo, AppenderColorType eColor)
{
	append(strInfo.c_str(), eColor);
}

void LinuxConsoleAppender::append(const char *szInfo, AppenderColorType eColor)
{
	string strCode;

	switch(eColor) {
	case ACT_RED:
		strCode = "\033[1;31m";
		break;
	case ACT_GREEN:
		strCode = "\033[1;32m";
		break;
	default:
		/* default set white */
		break;
	}
	fprintf(m_fp, "%s%s\033[1;37m\n", strCode.c_str(), szInfo);
	fflush(m_fp);
}

void LinuxConsoleAppender::append(const string &strInfo, AppenderColorType eColor)
{
	append(strInfo.c_str(), eColor);
}

FileAppender::FileAppender(const char *szLogPath) : Appender(NULL)
{
	m_strLogPath = szLogPath;
	m_fp = getConnection();
}

FileAppender::~FileAppender()
{
	endConnection();
}

FILE *FileAppender::getConnection()
{
	FILE *fp = fopen(m_strLogPath.c_str(), "a");
	if(!fp) {
		ErrorInfo err("FileAppender", "getConnection", "Failed to open required log file");
		throw runtime_error(err.get());
	}
	return fp;
}

void WindowsFileAppender::append(const char *szInfo)
{
	fprintf(m_fp, "%s\r\n", szInfo);
	fflush(m_fp);
}

void WindowsFileAppender::append(const string &strInfo)
{
	append(strInfo.c_str());
}

void LinuxFileAppender::append(const char *szInfo)
{
	fprintf(m_fp, "%s\n", szInfo);
	fflush(m_fp);
}

void LinuxFileAppender::append(const string &strInfo)
{
	append(strInfo.c_str());
}

void FileAppender::endConnection()
{
	if(m_fp) {
		fclose(m_fp);
		m_fp = NULL;
	}
}

Appender *WindowsAppenderFactory::getConsoleAppender()
{
	return new WindowsConsoleAppender();
}

Appender *WindowsAppenderFactory::getFileAppender(const char *szLogPath)
{
	return new WindowsFileAppender(szLogPath);
}

Appender *LinuxAppenderFactory::getConsoleAppender()
{
	return new LinuxConsoleAppender();
}

Appender *LinuxAppenderFactory::getFileAppender(const char *szLogPath)
{
	return new LinuxFileAppender(szLogPath);
}

Trace *Trace::m_pInstance = NULL;

Trace *Trace::getInstance(LogRankType eRank, bool bReset)
{
	if(bReset)
		destroyInstance();

	if(m_pInstance == NULL)
		m_pInstance = new Trace(NULL, LAT_CONSOLE, eRank);
	return m_pInstance;
}

Trace *Trace::getInstance(const char *szLogPath, LogAppenderType eLogApp,
			LogRankType eRank)
{
	destroyInstance();

	if(m_pInstance == NULL) {
		m_pInstance = new Trace(szLogPath, eLogApp, eRank);
	}

	return m_pInstance;
}

void Trace::destroyInstance()
{
	if(m_pInstance) {
		delete m_pInstance;
		m_pInstance = NULL;
	}
}

Trace::Trace(const char *szLogPath, LogAppenderType eLogApp, LogRankType eRank) :
		m_pConsole(NULL), m_pFile(NULL), m_eLogRank(eRank),
		m_eLogApp((szLogPath == NULL) ? LAT_CONSOLE : eLogApp)
{
	AppenderFactory *pAf = NULL;
#ifdef WIN32
	pAf = new WindowsAppenderFactory();
#else
	pAf = new LinuxAppenderFactory();
#endif
	if(m_eLogApp & LAT_CONSOLE)
		m_pConsole = pAf->getConsoleAppender();
	if(m_eLogApp & LAT_FILE)
		m_pFile = pAf->getFileAppender(szLogPath);

	if(pAf) {
		delete pAf;
		pAf = NULL;
	}
}

void Trace::changeLogRank(LogRankType eRank)
{
	m_eLogRank = eRank;
}

void Trace::alert(const string &strAlert)
{
	string strStr = "[alert] " + strAlert;
	showInfo(strStr, LRT_ALERT, ACT_RED);
}

void Trace::alert(const char *strFormat, ...)
{
	va_list ap;
	va_start(ap, strFormat);
	char szBuf[STR_BUF_SIZE];
	vsprintf(szBuf, strFormat, ap);
	alert(string(szBuf));
}

void Trace::info(const string &strInfo)
{
	string strStr = "[pSimulate] " + strInfo;
	showInfo(strStr, LRT_INFO);
}

void Trace::info(const char *strFormat, ...)
{
	va_list ap;
	va_start(ap, strFormat);
	char szBuf[STR_BUF_SIZE];
	vsprintf(szBuf, strFormat, ap);
	info(string(szBuf));
}

void Trace::debug(const string &strDebug)
{
	string strStr = "[debug] " + strDebug;
	showInfo(strStr, LRT_DEBUG, ACT_GREEN);
}

void Trace::debug(const char *strFormat, ...)
{
	va_list ap;
	va_start(ap, strFormat);
	char szBuf[STR_BUF_SIZE];
	vsprintf(szBuf, strFormat, ap);
	debug(string(szBuf));
}

void Trace::breakPoint()
{
	string strStr = "[stop] Press any key to continue ...";
	showInfo(strStr, LRT_DEBUG, ACT_RED);
	getchar();
}

Trace::~Trace()
{
	if(m_pConsole) {
		delete m_pConsole;
		m_pConsole = NULL;
	}
	if(m_pFile) {
		delete m_pFile;
		m_pFile = NULL;
	}
}

void Trace::showInfo(const string &strStr, LogRankType eRank, AppenderColorType eColor)
{
	if(m_eLogRank >= eRank) {
		if(m_pConsole)
			if(ConsoleAppender *pCa = dynamic_cast<ConsoleAppender *>(m_pConsole))
				pCa->append(strStr, eColor);
	}

	if(m_pFile)
		if(FileAppender *pFa = dynamic_cast<FileAppender *>(m_pFile))
			pFa->append(strStr);
}

} /* end of namespace util */
