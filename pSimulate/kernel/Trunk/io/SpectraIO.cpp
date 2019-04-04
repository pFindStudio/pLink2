#include "../include/io.h"
#include "../include/util.h"

using namespace std;

namespace sdk {


MS2OutputFactory::MS2OutputFactory() 
{
}

MS2OutputFactory::~MS2OutputFactory() 
{
}

MS2Output * MS2OutputFactory::getImporter(MS2FormatType eType)
{

	switch(eType) {
    case(MFT_MGF):
        return new MGFOutput;
	default:
        ErrorInfo err("MS2OutputFactory", "getImporter", "unkown input type.");
		throw runtime_error(err.get());
		return NULL;
	}
}


MGFOutput::MGFOutput():m_bMono(true), m_strOutPath(""), m_nOutCount(0), m_tMaxPeaksNum(MAX_PEAKS_NUM), m_ofOut(0)
{
	m_bMissCharge = false;
	m_buffString = "";
}

MGFOutput::~MGFOutput() {
	closeOutFile();
}

void MGFOutput::writeAll(const std::string &strPath, vector<Spectrum> & S)
{
	int nTot = 0;
	try
	{
		startWrite(strPath.c_str(),nTot);
	}
	catch(exception & e)
	{
		ErrorInfo info("CMGFOutput", "WriteAll", "StartLoad failed.");
		throw runtime_error(info.get(e));
	}
	catch(...)
	{
		ErrorInfo info("CMGFOutput", "WriteAll", "caught an unknown exception from StartLoad().");
		throw runtime_error(info.get());
	}
//	S.clear();
	Spectrum spec;
	int idx = 0;
	for(size_t i = 0;i < S.size();++i)
	{
//		S.pop_back(spec);
		spec = S[i];
		writeNext(spec, idx);
	}
	endWrite();
}

void MGFOutput::startWrite(const string &strPath,int & nTotal)
{
	try
	{
		openOutFile(strPath);
	}
	catch(exception & e)
	{
		ErrorInfo info("CMGFOutput", "StartWrite", "OpenOutFile: " + strPath + " failed.");
		throw runtime_error(info.get(e));
	}
	catch(...)
	{
		ErrorInfo info("CMGFOutput", "StartWrite", "caught an unknown exception from OpenInFile: " + strPath + " failed.");
		throw runtime_error(info.get());
	}

	m_vCharge.clear();
	writeFileHead();
	nTotal = 0;	//Currently set it to 0
	m_nOutCount = 0;
}

void MGFOutput::openOutFile(const string &strFilePath)
{
	if (( m_ofOut = fopen(strFilePath.c_str(), "wt+") )== NULL)
	{
		ErrorInfo info("CMGFOutput", "OpenOutFile", "m_ofOut.fail() when open file: " + strFilePath + " failed.");
		throw runtime_error(info.get());
	}
	m_strOutPath = strFilePath;
}

void MGFOutput::writeFileHead()
{

}

void MGFOutput::writeNext(Spectrum & spec, int & idx)
{
	writeMGFHead(spec);
	writeMZAndItensity(spec, true);
	fprintf(m_ofOut, "%s\n\n", "END IONS");

}

void MGFOutput::writeMGFHead(Spectrum& spec)
{
	fprintf(m_ofOut, "BEGIN IONS\n");
	fprintf(m_ofOut, "TITLE=%s\n", getTitle(spec).c_str());
	fprintf(m_ofOut, "SCANS=%d\n", spec.m_tScanNo);
	fprintf(m_ofOut, "CHARGE=%d+\n", (int)spec.m_tCharge);
	fprintf(m_ofOut, "PEPMASS=%lf\n", spec.m_lfMZ);


}

string MGFOutput::getTitle(Spectrum& spec)
{
	return spec.m_strTitle;
}

void MGFOutput::writeMZAndItensity(Spectrum& spec, bool bAdmitBlank)
{
	for (size_t i = 0; i < spec.m_vPeaks.size(); i++)
	{
		fprintf(m_ofOut, "%lf %lf\n", spec.m_vPeaks[i].m_lfMz, spec.m_vPeaks[i].m_lfIntensity);
	}

}

void MGFOutput::endWrite(void)
{
	closeOutFile();
}

void MGFOutput::closeOutFile(void)
{
	fclose(m_ofOut);
};

void MGFOutput::writeNextBatch(int nNumber ,vector<Spectrum> & S, int &idx)
{
	Spectrum spec;
	for(int i = 0; i < nNumber; ++i)
	{
		spec = S[i];
		writeNext(spec, idx);
	}
}
/* end: MGFOutput */

}
