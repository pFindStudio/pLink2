
#ifndef UTIL_H_
#define UTIL_H_
#include "config.h"

namespace sdk {
/* begin: types */

typedef char byte;
typedef long long Int64;
typedef signed long long SInt64;
typedef unsigned long long UInt64;


class EnumPeptideIndex;
typedef void (EnumPeptideIndex::*EnumPeptideIndexDispatch)();

/* end: types */

/* begin: constants */

#ifdef WIN64 // 64-bit Windows Compiler
#define FORMAT_SIZE_T "%zu"
#define FORMAT_SIZE_XT "%llx"
#define FORMAT_UINT64 "%llu"
const char SLASH = '\\';
#elif WIN32 // 32-bit Windows Compiler
#define FORMAT_SIZE_T "%u"
#define FORMAT_SIZE_XT "%x"
#define FORMAT_UINT64 "%I64d"
const char SLASH = '\\';
#else // 32-bit Linux Compiler
#define FORMAT_SIZE_T "%u"
#define FORMAT_SIZE_XT "%x"
#define FORMAT_UINT64 "%I64d"
const char SLASH = '/';
#endif

#ifdef WIN64
const size_t SHIFT_BITS = 56;
const size_t SET_MASK = 0x8000000000000000;
const size_t RECOVER_MASK = 0x7fffffffffffffff;
const size_t POSITION_MASK = 0x00ffffffffffffff;
const size_t LENGTH_MASK = 0x7f00000000000000;

const size_t SHIFT_HIGH_BITS = 32;
const size_t HIGH_MASK = 0xffffffff00000000;
const size_t LOW_MASK = 0x00000000ffffffff;

const size_t ALL_BIT_ON = 0xffffffffffffffff;

#else
const size_t SHIFT_BITS = 24;
const size_t SET_MASK = 0x80000000;
const size_t RECOVER_MASK = 0x7fffffff;
const size_t POSITION_MASK = 0x00ffffff;
const size_t LENGTH_MASK = 0x7f000000;

const size_t SHIFT_HIGH_BITS = 16;
const size_t HIGH_MASK = 0xffff0000;
const size_t LOW_MASK = 0x0000ffff;

const size_t ALL_BIT_ON = 0xffffffff;
#endif

const char N_TERM_MOD = '[';
const char C_TERM_MOD = ']';
const char PROTEIN_N_CHAR = '[';
const char PROTEIN_C_CHAR = ']';
const char PEPTIDE_N_CHAR = '(';
const char PEPTIDE_C_CHAR = ')';
const char PROTEIN_N_SITE = 'A' + 26;
const char PROTEIN_C_SITE = 'A' + 27;
const char PEPTIDE_N_SITE = 'A' + 28;
const char PEPTIDE_C_SITE = 'A' + 29;
const int GAP = 'a' - 'A';
const int BOTH_ION = 0;
const int ONLY_COMMON_ION = 1;
const int ONLY_LINKER_ION = 2;
const int FIXED_PRECISION = 6;
const size_t ONE_LOAD_PROTEINS = 1500;
const size_t AVERAGE_SQ_LEN = 500;
const size_t QUERY_PEAKS_NUM = 50;
const size_t QUERY_PEAK_CATEGORY = 4;
const size_t AA_NUM = 26;
const size_t IONTYPE_NUM = 26;
const size_t MZ_INDEX_SCALE = 1000;
const size_t MZ_MULTIPLIER = 10000;
const size_t FRAGMENT_ERROR_RANGE = 201;
const size_t MAX_IONTYPE_NUM = 30;
const size_t MAX_PEPTIDE_LEN = 100;
const size_t MAX_MODIFY_NUM = 20;
const size_t MAX_LINKER_NUM = 10;
const size_t LINK_SITE_NUM = 30;
const size_t MAX_PEAKS_NUM = 1024;
const size_t FEATURE_NUM = 10;
const size_t MAX_HASH_SIZE = 10000;
const size_t MAX_INNER_ION_LENGTH = 4;
const size_t STR_BUF_SIZE = 1024;
const size_t LOG_N_SIZE = 100000;
const size_t LOG_STDNORM_SIZE = 101;
const size_t LOG_TAGLEN_ROW_SIZE = 120;
const size_t LOG_TAGLEN_COL_SIZE = 21;
const size_t LOG_PRIME_SIZE = 2 * MAX_PEPTIDE_LEN + 5;
const size_t CONTINUITY_ION_NUM = 2;
const size_t N_ION_INDEX = 0;
const size_t C_ION_INDEX = 1;
const size_t INIT_AC_NODE = 100000;
const size_t PEPTIDE_TYPE_NUM = 4;
const size_t PROTEIN_TYPE_NUM = 3;
const size_t TD_TYPE_NUM = 3;
const size_t CHARGE_RANGE = 8;
const size_t MAX_PRECURSOR_ERROR = 50;
const size_t PRECURSOR_ERROR_INTERVAL_NUM = 200;
const size_t PRECURSOR_WIN_NUM = 5;
const double PRECISION = 0.000001;
const double PART_PER_MILLION = 0.000001;
const double MILLION = 1000000.0;
const double EMPERICAL_ERROR = 20.0;
const double PROTON_MASS = 1.00727647012;
const double H_MONO_MASS = 1.007825035;
const double H_AVRG_MASS = 1.00794;
const double C_MONO_MASS = 12;
const double C_AVRG_MASS = 12.0107;
const double N_MONO_MASS = 14.003074;
const double N_AVRG_MASS = 14.0067;
const double N15_ISTP_MASS = 15.0001088982;
const double O_MONO_MASS = 15.99491463;
const double O_AVRG_MASS = 15.9994;
const double NA_MONO_MASS = 22.98976967;
const double NA_AVRG_MASS = 22.98976967;
const double K_MONO_MASS = 38.9637069;
const double K_AVRG_MASS = 39.0983;
const double CA_MONO_MASS = 39.9625912;
const double CA_AVRG_MASS = 40.0784;
const double AVERAGE_MASS = 108.0;
const double CO_MONO_MASS = (C_MONO_MASS + O_MONO_MASS); // a, x
const double NH2_MONO_MASS = (N_MONO_MASS + H_MONO_MASS * 2); // z
const double NH3_MONO_MASS = (N_MONO_MASS + H_MONO_MASS * 3); // c
const double NH3_AVRG_MASS = (N_AVRG_MASS + H_AVRG_MASS * 3);
const double H2O_MONO_MASS = (O_MONO_MASS + H_MONO_MASS * 2);
const double H2O_AVRG_MASS = (O_AVRG_MASS + H_AVRG_MASS * 2);
const double UNKNOWN_QVALUE = 100.0;
const double HCD_BASES[PRECURSOR_WIN_NUM] = { 0, 1.003, 2.006, 3.009, 4.012 };
const double INTENSITY_MULTIPLIER = 100000.0;

const std::string ALL_AA = "ACDEFGHIKLMNPQRSTVWY";

const std::string NON_AA = "B"; // "BJOUXZ"

const std::string PSIMULATE_VERSION = "1.0";

const std::string PEPTIDE_TYPE_LABEL[] = { "Common", "Mono-Linked",
		"Loop-Linked", "Cross-Linked" };
const std::string PROTEIN_TYPE_LABEL[] = { "None", "Intra-Protein",
		"Inter-Protein" };
const std::string TD_TYPE_LABEL[] = { "Decoy", "TD-DT", "Target" };

/* end: constants */

/* begin: enumerations */

enum ModificationType {
	MT_PEP_NORMAL, // on body
	MT_PEP_N, // peptide N-terminal
	MT_PEP_C, // peptide C-terminal
	MT_PRO_N, // protein N-terminal
	MT_PRO_C // protein C-terminal
};

enum EnzymeType {
	ET_SPECIFIC = 0, ET_SEMI, ET_NONE
};

enum ToleranceType {
	TT_DA, // dalton
	TT_PPM // ppm
};

enum AppenderColorType {
	ACT_WHITE, ACT_RED, ACT_GREEN
};

enum LogRankType {
	LRT_NONE = 0, LRT_ALERT = 100, LRT_INFO = 200, LRT_DEBUG = 300
};

enum LogAppenderType {
	LAT_CONSOLE = 0x01, LAT_FILE = 0x02, LAT_CONSOLE_AND_FILE = 0x03
};

enum MS2FormatType {
	MFT_MGF,
	MFT_MS2,
	MFT_PF,
	MFT_RAW,
	MFT_DTA,
	MFT_DTAS,
	MFT_PKL,
	MFT_MZML,
	MFT_SDTA
};


enum ScoreType {
	ST_XLINK_LEGACY, ST_XLINK_PERFECT, ST_XLINK_SIMPLE
};

enum IndexStepType {
	IST_COUNT, IST_FILL
};

enum FlowType {
	FT_SIMULATION
};

/* do not change the values of PeptideType, ProteinType, TDType,
 * they will be used as index of an array */
enum PeptideType {
	PT_COMMON = 0, PT_MONO = 1, PT_LOOP = 2, PT_XLINK = 3
};

enum ProteinType {
	PRT_NONE = 0, PRT_INTRA = 1, PRT_INTER = 2
};

enum TDType {
	TDT_F = 0, // all from decoy
	TDT_U = 1, // one from target, the other from decoy
	TDT_T = 2 // all from target
};


/* end: enumerations */

/* begin: structures, to save space, attention: not to contain methods */

/* struct Isotope
 *   isotope
 *   e.g. C12 12.0000000 0.988930
 */
struct Isotope {
	double lfMass;
	double lfAbundance;
};

/* struct Element
 *   element
 *   vIsotopes[0] is always not empty
 */
struct Element {
	std::string strName;
	std::vector<Isotope> vIsotopes;
};

/* struct AminoAcid
 *   amino acid
 *   e.g. from aa.ini
 *   R1=A|C(3)H(5)N(1)O(1)S(0)|
 */
struct AminoAcid {
	char cAa; /* one letter abbreviation */
	double lfMono; /* mono mass */
	std::string strFormula; /* residue formula */
};

/* struct Enzyme
 *   enzyme
 *   e.g. from enzyme.ini
 *   Trypsin=_ _ KR _ 0
 *   strName=strNCleaveSite strNExceptSite strCCleaveSite strCExceptSite
 *   Underscore(_) means blank, so the corresponding string will be empty.
 */
struct Enzyme {
	std::string strName; /* name */
	std::string strNCleaveSite; /* N-terminal cleave site */
	std::string strNExceptSite; /* N-terminal exceptional cleave site */
	std::string strCCleaveSite; /* C-terminal cleave site */
	std::string strCExceptSite; /* C-terminal exceptional cleave site */
	EnzymeType eType; /* specific, semi-specific or non-specific */
};

/* struct Modification
 *   modification
 *   e.g. from modification.ini
 *   name1367=Trioxidation[C] 1
 *   Trioxidation[C]=C NORMAL 47.984744 47.9982 0 O(3)
 *
 *   strName=Trioxidation[C]
 *   strAA=C
 *   eModType=MT_PEP_NORMAL
 *   lfMonoDiff=47.984744
 *   lfArgvDiff=47.9982
 *   vMonoNeutralLossDiff and vAvrgNeutralLossDiff is empty
 *   ignore O(3), the chemical formula
 */
struct Modification {
	std::string strName;
	std::string strAA; /* amino acid to modify */
	ModificationType eModType; /* specificity */
	double lfMonoDiff;
	double lfAvrgDiff;
	std::vector<double> vMonoNeutralLossDiff;
	std::vector<double> vAvrgNeutralLossDiff;
};

/* struct XLinker
 *   linker
 *   e.g. from xlink.ini
 *   BS3=[K [K 138.068 138.068 156.079 156.079
 *   strName=strAlphaAA strBetaAA lfMonoMassDiff lfAvrgMassDiff lfMLMonoMassDiff lfMLAvrgMassDiff
 */
struct XLinker {
	std::string strName; /* linker name */
	std::string strAlphaAA; /* amino acid sites for alpha peptide */
	std::string strBetaAA; /* amino acid sites for beta peptide */
	double lfMonoMassDiff; /* x-linker mass shift in mono */
	double lfAvrgMassDiff; /* x-linker mass shift in average */
	double lfMLMonoMassDiff; /* mono-linker mass shift in mono */
	double lfMLAvrgMassDiff; /* mono-linker mass shift in average */
};

struct Sumo {
	std::string strName; /* sumo name, unique identifier*/
	std::string strSumoSq; /* fixed sumo sequence */
	std::string strSubstrateAA; /* substrate site */
	size_t tSumoSite; /* sumo sequence linked site, e.g. DQIEAVLEQLGG(12), tSumoSite = 11*/
	double lfSumoSqMass; /* sumo sequence mass */
	double lfMonoMass; /* x-linker mass */
};

struct IonType {
	char cType; /* ion identifier, a, b, y, etc. */
	int nCharge; /* charge state */
	int nNCoexistence; /* N terminal dissociation window subscript */
	int nCCoexistence; /* C terminal dissociation window subscript */
	int nHasLinker; /* 0 -- both, 1 -- common, 2 -- only linker */
	bool bContinuity; /* consider continuity or not */
	double lfLoss; /* extra mass */
};

/* struct Instrument
 *   optimized ion types
 *   e.g. from instrument.ini
 *   each line represents one IonType
 */
struct Instrument {
	std::vector<IonType> vIonTypes; /* ion types list */
};

/* end: structures */

/* begin: classes */
class ErrorInfo {
	std::string m_strClass; /* error occurred class */
	std::string m_strMethod; /* error occurred method */
	std::string m_strDetail; /* specific error */
	std::string m_strException; /* extra error caught from outside */

public:
	ErrorInfo(const std::string &strClass = "", const std::string &strMethod =
			"", const std::string &strDetail = "");
	ErrorInfo(const std::string &strClass, const std::string &strMethod,
			const std::string &strDetail, const std::exception &e);
	~ErrorInfo() {
	}
	void setErrorClass(const std::string &strClass);
	void setErrorMethod(const std::string &strMethod);
	void setErrorDetail(const std::string &strDetail);
	void setException(const std::exception &e);
	std::string get() const;
	std::string get(const std::exception &e);
	friend std::ofstream& operator<<(std::ofstream& os, const ErrorInfo& info);
	friend std::ostream& operator<<(std::ostream& os, const ErrorInfo& info);
};

/*
 * the structure of a section in a configuration file
 */
class SectionItem {
public:
	std::string strSctn; /* section name */
	std::vector<std::string> vOpts; /* options */
	std::map<std::string, std::string> mpOptVl; /* option and its corresponding value */
};

/*
 * This class is designed to parse configuration files
 * This implementation only suits for files of small size
 * Its implementation can be changed by writing a derived class
 * and override those virtual methods below
 *
 * The format of a configuration file:
 * [section]
 * option=value, blanks besides equal sign will be omitted
 * a new line beginning with ; or # regards as a comment
 */
class ConfigParser {
protected:
	std::vector<std::string> m_vSctns; /* keep the order of sections */
	std::map<std::string, SectionItem> m_mpSctns; /* section map */

	std::string m_strFilePath; /* configuration file path */
	std::fstream m_fin; /* file reader */

	SectionItem m_stCurrent; /* current section */
	SectionItem m_stNext; /* next section */

public:
	static void testCase();

	ConfigParser();
	ConfigParser(const char *szFilePath);

	virtual ~ConfigParser();
	virtual void read(const char *szFilePath) throw (std::runtime_error);
	virtual bool hasSection(const std::string &strSctn) const;
	virtual bool hasOption(const std::string &strSctn,
			const std::string &strOpt) const;
	virtual std::vector<std::string> sections() const;
	virtual std::vector<std::string> options(const std::string &strScnt) const;
	virtual std::string get(const std::string &strScnt,
			const std::string &strOpt) const;
	virtual bool getBool(const std::string &strScnt,
			const std::string &strOpt) const;
	virtual int getInt(const std::string &strScnt,
			const std::string &strOpt) const;
	virtual double getDouble(const std::string &strScnt,
			const std::string &strOpt) const;

protected:
	virtual void parse();
	virtual void close();

	void first();
	void next();
	bool hasNext();
};

class NoSectionConfigParser : public ConfigParser {
public:
	static const std::string NO_SECTION;

	NoSectionConfigParser();
	virtual ~NoSectionConfigParser();

protected:
	NoSectionConfigParser(const char *szFilePath);
	virtual void parse();
};

class OneByOneConfigParser: public ConfigParser {
protected:
	std::string strEmpty;

public:
	OneByOneConfigParser();
	virtual ~OneByOneConfigParser();

	// overloading
	bool hasOption(const std::string &strOpt) const;
	std::vector<std::string> options() const;
	const std::string &section() const;
	std::string get(const std::string &strOpt) const;
	bool getBool(const std::string &strOpt) const;
	int getInt(const std::string &strOpt) const;
	double getDouble(const std::string &strOpt) const;

	// name hiding, no polymorphism
	void first();
	void next();
	bool hasNext();

protected:
	OneByOneConfigParser(const char *szFilePath);
	virtual void parse();
};

/* class ElementDict
 *   load and parse element.ini
 */
class ElementDict {
	static ElementDict *m_pInstance;
	std::string m_strConfigFilePath;
	std::map<std::string, Element> m_mpElementMap;

public:
	static ElementDict *getInstance(const char *szConfigFilePath);
	static ElementDict *getInstance() throw (std::runtime_error);
	static void destroyInstance();

	virtual ~ElementDict();

	const Element &getElement(const std::string &strElementName) const;
	double getMonoMass(const std::string &strElementName) const;
	double getMonoAbundance(const std::string &strElementName) const;
	Isotope getMonoIsotope(const std::string &strElementName) const;
	double getMass(const std::string &strElementName, size_t tIdx) const;
	double getAbundance(const std::string &strElementName, size_t tIdx) const;
	Isotope getIsotope(const std::string &strElementName, size_t tIdx) const;
	double calculateMolecularMass(const std::string &strFormula);

protected:
	ElementDict(const char *szConfigFilePath);
	virtual void parse();
};

/* class AminoAcidDict
 *   load and parse aa.ini
 *   each (key, value) pair stored in m_mpAAMap
 */
class AminoAcidDict {
	static AminoAcidDict *m_pInstance;
	std::string m_strConfigFilePath;
	std::map<char, AminoAcid> m_mpAAMap;

public:
	static AminoAcidDict *getInstance(const char *szConfigFilePath); /* for the first time initialization */
	static AminoAcidDict *getInstance() throw (std::runtime_error); /* for afterwards */
	static void destroyInstance();

	virtual ~AminoAcidDict();

	const AminoAcid &getAminoAcid(char cAA) const;
	double getMonoMass(char cAA) const;
	double calculatePeptideMass(const std::string &strSq);

protected:
	AminoAcidDict(const char *szConfigFilePath);
	virtual void parse();
};


/* class EnzymeDict
 *   load and parse enzyme.ini
 *   (key, value) pair stored in m_mpEnzymeMap
 */
class EnzymeDict {
	static EnzymeDict *m_pInstance;
	std::string m_strConfigFilePath;
	std::map<std::string, Enzyme> m_mpEnzymeMap;

public:
	static EnzymeDict *getInstance(const char *szConfigFilePath); /* for the first time initialization */
	static EnzymeDict *getInstance() throw (std::runtime_error); /* for afterwards */
	static void destroyInstance();

	virtual ~EnzymeDict();

	const Enzyme &getEnzyme(const std::string &strEnzymeName) const;
	std::string getNCleaveSite(const std::string &strEnzymeName) const;
	std::string getNExceptSite(const std::string &strEnzymeName) const;
	std::string getCCleaveSite(const std::string &strEnzymeName) const;
	std::string getCExceptSite(const std::string &strEnzymeName) const;
	EnzymeType getEnzymeType(const std::string &strEnzymeName) const;

protected:
	EnzymeDict(const char *szConfigFilePath);
	virtual void parse();
};

/* class ModificationDict
 *   load and parse modification.ini
 *   (key, value) pair stored in m_mpModMap
 */
class ModificationDict {
	static ModificationDict *m_pInstance;
	std::string m_strConfigFilePath;
	std::map<std::string, Modification> m_mpModMap;

public:
	static ModificationDict *getInstance(const char *szConfigFilePath);
	static ModificationDict *getInstance() throw (std::runtime_error);
	static void destroyInstance();

	virtual ~ModificationDict();

	const Modification &getModification(const std::string &strModName) const;
	const std::string &getAminoAcid(const std::string &strModName) const;
	ModificationType getModificationType(const std::string &strModName) const;
	double getMonoDiff(const std::string &strModName) const;
	double getAvrgDiff(const std::string &strModName) const;
	std::vector<double> getMonoNeutralLossDiff(
			const std::string &strModName) const;
	std::vector<double> getAvrgNeutralLossDiff(
			const std::string &strModName) const;

protected:
	ModificationDict(const char *szConfigFilePath);
	virtual void parse();
};

/* class ModificationDict
 *   load and parse modification.ini
 *   (key, value) pair stored in m_mpLinkerMap
 */
class XLinkerDict {
	static XLinkerDict *m_pInstance;
	std::string m_strConfigFilePath;
	std::map<std::string, XLinker> m_mpLinkerMap;

public:
	static XLinkerDict *getInstance(const char *szConfigFilePath);
	static XLinkerDict *getInstance() throw (std::runtime_error);
	static void destroyInstance();

	virtual ~XLinkerDict();

	const XLinker &getXLinker(const std::string &strLinkerName) const;
	const std::string &getAlphaAA(const std::string &strLinkerName) const;
	const std::string &getBetaAA(const std::string &strLinkerName) const;
	double getMonoDiff(const std::string &strLinkerName) const;
	double getAvrgDiff(const std::string &strLinkerName) const;
	double getMLMonoDiff(const std::string &strLinkerName) const;
	double getMLAvrgDiff(const std::string &strLinkerName) const;

protected:
	XLinkerDict(const char *szConfigFilePath);
	virtual void parse();
};

/* class InstrumentDict
 *   load and parse modification.ini
 *   (key, value) pair stored in m_mpInstrumentMap
 */
class InstrumentDict {
	static InstrumentDict *m_pInstance;
	std::string m_strConfigFilePath;
	std::map<std::string, Instrument> m_mpInstrumentMap;

public:
	static InstrumentDict *getInstance(const char *szConfigFilePath);
	static InstrumentDict *getInstance() throw (std::runtime_error);
	static void destroyInstance();

	virtual ~InstrumentDict();

	const Instrument &getInstrument(const std::string &strInstrumentName) const;
	const std::vector<IonType> &getIonTypes(
			const std::string &strInstrumentName) const;

protected:
	InstrumentDict(const char *szConfigFilePath);
	virtual void parse();
};

class Appender {
protected:
	FILE *m_fp;

public:
	Appender(FILE *fp) :
			m_fp(fp) {
	}
	virtual ~Appender() {
		if (m_fp) {
			if (m_fp != stdout)
				delete m_fp;
			m_fp = NULL;
		}
	}

protected:
	virtual FILE *getConnection() = 0;
	virtual void endConnection() = 0;
};

class ConsoleAppender: public Appender {
public:
	ConsoleAppender() :
			Appender(stdout) {
	}
	virtual ~ConsoleAppender() {
	}

	virtual void append(const char *szInfo,
			AppenderColorType eColor = ACT_WHITE) = 0;
	virtual void append(const std::string &strInfo, AppenderColorType eColor =
			ACT_WHITE) = 0;

protected:
	virtual FILE *getConnection() {
		return stdout;
	}
	virtual void endConnection() {
	}
};

class WindowsConsoleAppender: public ConsoleAppender {
public:
	WindowsConsoleAppender() {
	}
	virtual ~WindowsConsoleAppender() {
	}

	virtual void append(const char *szInfo,
			AppenderColorType eColor = ACT_WHITE);
	virtual void append(const std::string &strInfo, AppenderColorType eColor =
			ACT_WHITE);
};

class LinuxConsoleAppender: public ConsoleAppender {
public:
	LinuxConsoleAppender() {
	}
	virtual ~LinuxConsoleAppender() {
	}

	virtual void append(const char *szInfo,
			AppenderColorType eColor = ACT_WHITE);
	virtual void append(const std::string &strInfo, AppenderColorType eColor =
			ACT_WHITE);
};

class FileAppender: public Appender {
protected:
	std::string m_strLogPath;

public:
	FileAppender(const char *szLogPath);
	virtual ~FileAppender();

	virtual void append(const char *szInfo) = 0;
	virtual void append(const std::string &strInfo) = 0;

protected:
	virtual FILE *getConnection();
	virtual void endConnection();
};

class WindowsFileAppender: public FileAppender {
public:
	WindowsFileAppender(const char *szLogPath) :
			FileAppender(szLogPath) {
	}
	virtual ~WindowsFileAppender() {
	}

	virtual void append(const char *szInfo);
	virtual void append(const std::string &strInfo);
};

class LinuxFileAppender: public FileAppender {
public:
	LinuxFileAppender(const char *szLogPath) :
			FileAppender(szLogPath) {
	}
	virtual ~LinuxFileAppender() {
	}

	virtual void append(const char *szInfo);
	virtual void append(const std::string &strInfo);
};

class AppenderFactory {
public:
	AppenderFactory() {
	}
	virtual ~AppenderFactory() {
	}

	virtual Appender *getConsoleAppender() = 0;
	virtual Appender *getFileAppender(const char *LogPath) = 0;
};

class WindowsAppenderFactory: public AppenderFactory {
public:
	WindowsAppenderFactory() {
	}
	virtual ~WindowsAppenderFactory() {
	}

	virtual Appender *getConsoleAppender();
	virtual Appender *getFileAppender(const char *szLogPath);
};

class LinuxAppenderFactory: public AppenderFactory {
public:
	LinuxAppenderFactory() {
	}
	virtual ~LinuxAppenderFactory() {
	}

	virtual Appender *getConsoleAppender();
	virtual Appender *getFileAppender(const char *szLogPath);
};

class Trace {
	static Trace * m_pInstance;

	Appender *m_pConsole;
	Appender *m_pFile;
	LogRankType m_eLogRank;
	LogAppenderType m_eLogApp;

public:
	static Trace *getInstance(LogRankType eRank = LRT_INFO,
			bool bReset = false);
	static Trace *getInstance(const char *szLogPath, LogAppenderType eLogApp =
			LAT_FILE, LogRankType eRank = LRT_INFO);
	static void destroyInstance();

	virtual ~Trace();

	virtual void changeLogRank(LogRankType eRank);
	virtual void alert(const std::string &strAlert);
	virtual void alert(const char *strFormat, ...);
	virtual void info(const std::string &strInfo);
	virtual void info(const char *strFormat, ...);
	virtual void debug(const std::string &strDebug);
	virtual void debug(const char *strFormat, ...);
	virtual void breakPoint();

protected:
	Trace(const char *szLogPath, LogAppenderType eLogApp, LogRankType eRank =
			LRT_INFO);
	void showInfo(const std::string &strStr, LogRankType eRank,
			AppenderColorType eColor = ACT_WHITE);
};

struct PrecursorWindow {
	size_t m_tIndex;
	double m_lfMinMass;
	double m_lfMaxMass;

	PrecursorWindow();
	PrecursorWindow(size_t, double, double);

	friend bool operator<(const PrecursorWindow &pw1,
			const PrecursorWindow &pw2) {
		return pw1.m_lfMinMass < pw2.m_lfMinMass;
	}

	friend bool operator<=(const PrecursorWindow &pw1,
			const PrecursorWindow &pw2) {
		return pw1.m_lfMinMass <= pw2.m_lfMinMass;
	}

	friend bool operator>(const PrecursorWindow &pw1,
			const PrecursorWindow &pw2) {
		return pw1.m_lfMinMass > pw2.m_lfMinMass;
	}

	friend bool operator>=(const PrecursorWindow &pw1,
			const PrecursorWindow &pw2) {
		return pw1.m_lfMinMass >= pw2.m_lfMinMass;
	}
};

struct PeptideTol {
	double m_lfPeptideTol;
	double m_lfPeptideTolBase;
	ToleranceType m_ePeptideTolType;
	ToleranceType m_ePeptideTolBaseType;

	PeptideTol();
	friend bool operator<(const PeptideTol &tol1, const PeptideTol &tol2) {
		return tol1.m_lfPeptideTolBase < tol2.m_lfPeptideTolBase;
	}
};

struct DoubleItem {
	size_t m_tOrder;
	double m_lfNumber;

	DoubleItem();
	DoubleItem(size_t tOrder, double lfNumber);
	~DoubleItem() {
	}

	friend bool operator<(const DoubleItem &stItem1, const DoubleItem &stItem2);
	friend bool operator<=(const DoubleItem &stItem1,
			const DoubleItem &stItem2);
	friend bool operator>(const DoubleItem &stItem1, const DoubleItem &stItem2);
	friend bool operator>=(const DoubleItem &stItem1,
			const DoubleItem &stITem2);
	friend bool operator==(const DoubleItem &stItem1,
			const DoubleItem &stItem2);
	friend bool operator!=(const DoubleItem &stItem1,
			const DoubleItem &stItem2);
};

struct SearchParameter {

	// config
	std::string m_strElementFilePath;
	std::string m_strAAListFilePath;
	std::string m_strEnzymeListFilePath;
	std::string m_strModificationFilePath;
	std::string m_strLinkerFilePath;
	std::string m_strInstrumentFilePath;

	// database
	std::string m_strDBName;
	std::string m_strFastaFilePath;
	std::string m_strEnzymeName;
	size_t m_tMaxMissSite;
	size_t m_tMinPepLen;
	size_t m_tMaxPepLen;
	size_t m_tMinFragLen;
	size_t m_tMaxFragLen;
	double m_lfMinPepMass;
	double m_lfMaxPepMass;
	double m_lfMinFragMass;
	double m_lfMaxFragMass;
	bool m_bAutoReverse;
	size_t m_tMaxMemSize;

	// modification
	std::vector<std::string> m_vFixedMods;
	std::vector<std::string> m_vVariableMods;
	size_t m_tMaxModsNo;

	// linker
	std::vector<std::string> m_vLinkers;

	// ions
	std::string m_strRefinedInstrument;
	double m_lfFragTol;
	ToleranceType m_eFragTolType;

	double m_lfPepTol;
	ToleranceType m_ePepTolType;

	// score
	ScoreType m_eRefinedScoreType;

	// log
	LogRankType m_eLogRank;

	// flow
	FlowType m_eFlowType;
	std::vector<size_t> m_tSpectraNo;

	size_t m_tProcessorNo;
	size_t m_tProteinsBatchSize; // one load proteins
	std::string m_strDecoySign;
	std::string m_strOutputDirPath;
	int m_nLostPeakPercentage; // 丢峰的概率

	// spectrum
	std::string m_strSpecTitle;
	MS2FormatType m_eSpecType;
	std::vector<std::string> m_vSpecFilePath;

	// other useful information
	std::string m_strWorkDir;

	SearchParameter();
	virtual ~SearchParameter();

	double getError(double lfThr, double lfExp) const;
	size_t getErrorInterval(double lfError) const;
};

class ParameterReader {
	static ParameterReader *m_pInstance;

	std::string m_strFilePath;
	SearchParameter *m_pParameter;

public:
	static ParameterReader *getInstance(const std::string &strFilePath);
	static ParameterReader *getInstance();
	static void destroyInstance();

	virtual ~ParameterReader();

	virtual void loadParameter();
	virtual SearchParameter *getParameter();

protected:
	ParameterReader(const std::string &strFilePath);
};

struct ThreadState {
	size_t m_tID;
	size_t m_tCurrentSpectra;
	size_t m_tTotalSpectra;

	ThreadState();
};


class IonIndexState {
	clock_t m_tStartClock;
	time_t m_tStartTime;
	size_t m_tTotalProteins;
	size_t m_tCurrentProteins;

public:
	IonIndexState();
	virtual ~IonIndexState();

	virtual std::string getProgress();
	virtual void setClock();
	virtual void setTotalProteins(size_t tTotal);
	virtual void setCurrentProteins(size_t tCurrent);
	virtual bool isFinished() const;
	virtual clock_t getStartClock() const;
	virtual time_t getStartTime() const;
};

class ThreadMonitor {
	std::vector<bool> m_vFlowFlag;
	size_t m_tCurrentThreadID;
	bool m_bLocked;

public:
	pthread_mutex_t m_stMutex;
	pthread_cond_t m_stCond;

public:
	ThreadMonitor(size_t tThreadNo);
	~ThreadMonitor();

	bool isAllBusy() {
		for (size_t i = 0; i < m_vFlowFlag.size(); ++i) {
			if (m_vFlowFlag[i]) {
				m_tCurrentThreadID = i;
				return false;
			}
		}
		return true;
	}

	void setSignal(size_t tID) {
		m_vFlowFlag[tID] = false;
	}

	void freeSignal(size_t tID) {
		m_vFlowFlag[tID] = true;
	}

	bool isFreeSignal(size_t tID) const {
		return m_vFlowFlag[tID];
	}

	void setLock() {
		m_bLocked = true;
	}

	void freeLock() {
		m_bLocked = false;
	}

	bool isLocked() const {
		return m_bLocked;
	}

	void setCurrentID(size_t tID) {
		m_tCurrentThreadID = tID;
	}

	size_t getCurrentID() const {
		return m_tCurrentThreadID;
	}

	void waitForLock() {
		while (m_bLocked) {
#ifdef WIN32
			Sleep(100);
#else
			usleep(100);
#endif
		}
	}
};

class ConstantStore {
	static ConstantStore *m_pInstance;

	std::string m_strWorkDir;
	double m_lfAAMass[AA_NUM];
	double m_lfLogPrime[LOG_PRIME_SIZE];
	double m_lfLogNPerm[LOG_N_SIZE]; // for factorial

public:
	friend class Peptide;
	friend class XLinkPeptideResult;
	friend class MzCalculator;
	friend class OpenProbScorer;
	friend class RefinedProbScorer;
	friend class SumoFlow;
	static ConstantStore *getInstance();
	static void destroyInstance();

	virtual ~ConstantStore();

protected:
	ConstantStore();

private:
	void initAAMass();
	void initLogPrime();
	void initLogNPerm();
	bool isPrime(int nNumber);
};

class EnvironmentInitializer {
	static EnvironmentInitializer *m_pInstance;
	SearchParameter *m_pParameter;
	Trace *m_pTrace;

public:
	static void init(const char *szParamPath);
	static void destroy();
	static SearchParameter *getSearchParameter();
	virtual ~EnvironmentInitializer();

protected:
	EnvironmentInitializer(std::string &strParamPath);

	void createDirs();
	void reverseDatabase();
};

struct WindowsHash {
	int *m_pData;
	int m_nStart;
	double m_lfMaxMass;

	WindowsHash();
	~WindowsHash();

	void construct(const std::vector<PrecursorWindow> &vInfo);
	int findLowerBound(double lfMass);
	bool isBeyondBound(const std::vector<PrecursorWindow> &vInfo, double lfMass);
};

/* end: classes */

/* begin: functions */
/* trim from start */
static inline std::string &ltrim(std::string &s) {
	s.erase(s.begin(),
			std::find_if(s.begin(), s.end(),
					std::not1(std::ptr_fun<int, int>(std::isspace))));
	return s;
}

/* trim from end */
static inline std::string &rtrim(std::string &s) {
	s.erase(
			std::find_if(s.rbegin(), s.rend(),
					std::not1(std::ptr_fun<int, int>(std::isspace))).base(),
			s.end());
	return s;
}

/* trim from both ends */
static inline std::string &trim(std::string &s) {
	return ltrim(rtrim(s));
}

/* convert to upper case */
static inline std::string &toupper(std::string &s) {
	for (std::string::iterator it = s.begin(); it != s.end(); ++it) {
		if (*it >= 'a' && *it <= 'z') {
			*it = (char)(*it - GAP);
		}
	}
	return s;
}

/* convert to lower case */
static inline std::string &tolower(std::string &s) {
	for (std::string::iterator it = s.begin(); it != s.end(); ++it) {
		if (*it >= 'A' && *it <= 'Z') {
			*it = (char)(*it + GAP);
		}
	}
	return s;
}

// 函数将vStrs中的字符串用delim分隔符连接在一起,与split是相反的过程
static inline void join(std::string &strJoined, char delim,
		const std::vector<std::string> &vStrs)
{
	strJoined = "";
	for(size_t tIdx = 0; tIdx < vStrs.size(); ++tIdx) {
		strJoined += vStrs[tIdx];
		if(tIdx + 1 < vStrs.size())
			strJoined += delim;
	}
}

static inline void join(std::string &strJoined, const std::string &delim,
		const std::vector<std::string> &vStrs)
{
	strJoined = "";
	for(size_t tIdx = 0; tIdx < vStrs.size(); ++tIdx) {
		strJoined += vStrs[tIdx];
		if(tIdx + 1 < vStrs.size())
			strJoined += delim;
	}
}

// 函数将str以delim作为分隔符切分后的子串放在vec中
// 其中delim在字符串中的连续出现被看做单个的分隔符
static inline void split(std::string &str, const char delim,
		std::vector<std::string> &vec) {
	size_t last, index; // last记录子串的其实位置 index记录子串终止位置的下一个

	last = 0;

	do {
		while (str[last] == delim) // 滤掉连续出现的delim last指向新的子串起始位置
			last++;
		index = str.find_first_of(delim, last); // 找到下一个delim出现的位置记录在index中
		if (index == std::string::npos) // 找到了字符串末尾
			break;
		else {
			vec.push_back(str.substr(last, index - last));
			last = index + 1;
		}
	} while (true);
	if (last < str.length()) // 末尾出现连续的delim
		vec.push_back(str.substr(last, str.length() - last));
}

static inline size_t getCPUCoreNo() {
#ifdef WIN32
	SYSTEM_INFO info;
	GetSystemInfo(&info);
	return info.dwNumberOfProcessors;
#elif LINUX
	return 1;
//	return get_nprocs_conf();
#else
#error
#endif
}

static inline std::string getFile(const std::string &strFilePath) {
	size_t tStart = strFilePath.find_last_of(SLASH);
	size_t tEnd = strFilePath.length();
	return strFilePath.substr(tStart + 1, tEnd - tStart - 1);
}

static inline std::string getFilePath(const std::string &strFilePath) {
	size_t tStart = strFilePath.find_last_of(SLASH);
	return strFilePath.substr(0, tStart+1);
}

static inline void makeDir(const std::string &strDirPath) {
	Trace *pTrace = Trace::getInstance();
	if (access(strDirPath.c_str(), 0) == -1) {
//		pTrace->alert(strDirPath + " does not exist.");
#ifdef WIN32
		if (mkdir(strDirPath.c_str()) == -1) {
#else
			if(mkdir(strDirPath.c_str(), 777) == -1) {
#endif
			ErrorInfo err("util", "makeDir",
					"failed to create dirctory " + strDirPath);
			throw std::runtime_error(err.get());
		} else {
			pTrace->debug(getFile(strDirPath) + " has been created.");
		}
	}
}

static inline std::string getFilename(const std::string &strFilePath) {
	size_t tStart = strFilePath.find_last_of(SLASH);
	size_t tEnd = strFilePath.find_last_of('.');
	return strFilePath.substr(tStart + 1, tEnd - tStart - 1);
}

/* end: functions */

} /* end of namespace util */

#endif /* UTIL_H_ */
