#include "../include/io.h"
#include "../include/util.h"

using namespace std;

namespace sdk
{

SimplePSMIO::SimplePSMIO() : m_ofOut(0){
}

SimplePSMIO::~SimplePSMIO()
{
}

void SimplePSMIO::startWrite(const string &strPath){
	try
	{
		openOutFile(strPath);
	}
	catch(exception & e)
	{
		ErrorInfo info("SimplePSMIO", "startWrite", "OpenOutFile: " + strPath + " failed.");
		throw runtime_error(info.get(e));
	}
	catch(...)
	{
		ErrorInfo info("SimplePSMIO", "startWrite", "caught an unknown exception from OpenInFile: " + strPath + " failed.");
		throw runtime_error(info.get());
	}

}

void SimplePSMIO::openOutFile(const string &strFilePath)
{
	if (( m_ofOut = fopen(strFilePath.c_str(), "wt+") )== NULL)
	{
		ErrorInfo info("SimplePSMIO", "openOutFile", "m_ofOut.fail() when open file: " + strFilePath + " failed.");
		throw runtime_error(info.get());
	}
	m_strOutPath = strFilePath;
}

void SimplePSMIO::endWrite(void)
{
	fclose(m_ofOut);
}

void SimplePSMIO::writeNext(Spectrum &spec, XLinkPeptideResult &stXLinkPep){
	fprintf(m_ofOut, "%s,", spec.m_strTitle.c_str());
	fprintf(m_ofOut, "%s,", stXLinkPep.m_stAlphaPep.m_stPep.m_strSq.c_str());
	fprintf(m_ofOut, "%s,", stXLinkPep.m_stBetaPep.m_stPep.m_strSq.c_str());
	fprintf(m_ofOut, "%d,", stXLinkPep.m_stAlphaPep.m_stLinkMod.m_nSite);
	fprintf(m_ofOut, "%d,", stXLinkPep.m_stBetaPep.m_stLinkMod.m_nSite);
	fprintf(m_ofOut, "%d,", stXLinkPep.m_ePepType);
	fprintf(m_ofOut, "%d,", stXLinkPep.m_stAlphaPep.m_stPep.m_vMods.size());
	fprintf(m_ofOut, "%d,\n", stXLinkPep.m_stBetaPep.m_stPep.m_vMods.size());
}

}



