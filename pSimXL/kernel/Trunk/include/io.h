#ifndef IO_H_
#define IO_H_
#include "sdk.h"
#include "util.h"

namespace sdk {


class SimplePSMIO {
protected:
	FILE *m_ofOut;
	std::string m_strOutPath;

public:
	SimplePSMIO();
	virtual ~SimplePSMIO();
	virtual void startWrite(const std::string &strFilePath);
	virtual void openOutFile(const std::string &strFilePath);
	virtual void endWrite(void);
	virtual void writeNext(Spectrum &spec, XLinkPeptideResult &stXLinkPep);
};


class MS2Output {
public:
	MS2Output(void) {
	}
	;
	virtual ~MS2Output() {
	}
	;

	virtual MS2FormatType getType(void) = 0;
	virtual void writeAll(const std::string &strPath,
			std::vector<Spectrum> &S) = 0;
	virtual void startWrite(const std::string &strPath, int & nTotal) = 0;
	virtual void writeNext(Spectrum &spec, int &idx) = 0;
	virtual void writeNextBatch(int nNumber, std::vector<Spectrum> & S,
			int &idx) = 0;
	virtual void endWrite(void) = 0;
};

class MS2OutputFactory {
public:
	MS2OutputFactory();
	virtual ~MS2OutputFactory();

	MS2Output * getImporter(MS2FormatType eType);
};

class MGFOutput: public MS2Output {
protected:
	bool m_bMono;
	FILE * m_ofOut;
	std::ifstream m_ifIn;
	std::string m_strOutPath;
	int m_nOutCount;
	size_t m_tMaxPeaksNum;

	////用于文件头有CHARGE信息
	std::vector<int> m_vCharge; //CHARGE的个数。以AND分开，如2+ and 3+
	std::vector<int> m_vCurrCharge; //当前的CHARGE
	bool m_bMissCharge; //看每个谱是否缺失CHARGE，缺失则以m_vCharge中的数目填补,将当前行记录到m_strMissCharge中。
	std::string m_buffString; //用于记录读到的第一对谱图数目项
	std::string m_strMissCharge;
	Spectrum m_specTemp;
	////END用于文件头有CHARGE信息

public:
	MGFOutput();
	virtual ~MGFOutput();
	virtual void writeAll(const std::string &strPath,
			std::vector<Spectrum> & S);
	virtual void startWrite(const std::string &strPath, int & nTotal);
	virtual void writeNext(Spectrum & spec, int & idx);
	virtual void writeNextBatch(int nNumber, std::vector<Spectrum> & S,
			int &idx);
	virtual void endWrite(void);

	virtual MS2FormatType getType(void) {
		return MFT_MGF;
	}
	;

protected:
	void openOutFile(const std::string &strFilePath);
	void closeOutFile(void);
	void writeMGFHead(Spectrum& spec);
	void writeMZAndItensity(Spectrum& spec, bool bAdmitBlank = false);
	void writeFileHead();
	std::string getTitle(Spectrum& spec);
};


} /* end of namespace sdk */

#endif /* IO_H_ */
