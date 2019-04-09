
#include "../include/sdk.h"

using namespace std;

namespace sdk {

IonPeak::IonPeak(double lfMz, double lfIntensity, size_t tCharge) :
		m_lfMz(lfMz), m_lfIntensity(lfIntensity), m_nMz(
				MZ_MULTIPLIER * lfMz + 0.5), m_tCharge(tCharge) {
}

IonPeak::IonPeak(const IonPeak &stPeak) {
	m_lfMz = stPeak.m_lfMz;
	m_lfIntensity = stPeak.m_lfIntensity;
	m_nMz = MZ_MULTIPLIER * m_lfMz + 0.5;
	m_tCharge = stPeak.m_tCharge;
}

IonPeak &IonPeak::operator=(const IonPeak &stPeak) {
	m_lfMz = stPeak.m_lfMz;
	m_lfIntensity = stPeak.m_lfIntensity;
	m_nMz = MZ_MULTIPLIER * m_lfMz + 0.5;
	m_tCharge = stPeak.m_tCharge;
	return *this;
}

QueryPeak::QueryPeak(size_t tHighIdx, size_t tLowIdx) :
		m_tLowIdx(tLowIdx), m_tHighIdx(tHighIdx) {
}

bool Spectrum::intesityGreater(const IonPeak &stPeak1, const IonPeak &stPeak2) {
	return stPeak1.m_lfIntensity > stPeak2.m_lfIntensity;
}

bool Spectrum::intesityLesser(const IonPeak &stPeak1, const IonPeak &stPeak2) {
	return stPeak1.m_lfIntensity < stPeak2.m_lfIntensity;
}

bool Spectrum::mzGreater(const IonPeak &stPeak1, const IonPeak &stPeak2) {
	return stPeak1.m_lfMz > stPeak2.m_lfMz;
}

bool Spectrum::mzLesser(const IonPeak &stPeak1, const IonPeak &stPeak2) {
	return stPeak1.m_lfMz < stPeak2.m_lfMz;
}

Spectrum::Spectrum() :
		m_lfMH(0.0), m_lfMZ(0.0), m_lfTotalInt(0.0), m_lfSqrtMaxInt(0.0), m_tCharge(
				0), m_tScanNo(0) {
}

ostream &operator<<(ostream &stOut, const Spectrum &stSpec) {
	stOut << fixed << setprecision(FIXED_PRECISION);
	stOut << stSpec.m_strTitle << " " << stSpec.m_tCharge << " "
			<< stSpec.m_lfSqrtMaxInt << " " << stSpec.m_lfMH << " "
			<< stSpec.m_lfMZ;
	return stOut;
}

fstream &operator<<(fstream &stOut, const Spectrum &stSpec) {
	stOut << fixed << setprecision(FIXED_PRECISION);
	stOut << stSpec.m_strTitle << " " << stSpec.m_tCharge << " "
			<< stSpec.m_lfSqrtMaxInt << " " << stSpec.m_lfMH << " "
			<< stSpec.m_lfMZ;
	return stOut;
}

std::istream &operator>>(std::istream &stIn, Spectrum &stSpec) {
	stIn >> stSpec.m_strTitle >> stSpec.m_tCharge >> stSpec.m_lfSqrtMaxInt
			>> stSpec.m_lfMH >> stSpec.m_lfMZ;
	return stIn;
}

std::fstream &operator>>(std::fstream &stIn, Spectrum &stSpec) {
	stIn >> stSpec.m_strTitle >> stSpec.m_tCharge >> stSpec.m_lfSqrtMaxInt
			>> stSpec.m_lfMH >> stSpec.m_lfMZ;
	return stIn;
}

void Spectrum::copyBasicItems(const Spectrum &stSpec) {
	m_strTitle = stSpec.m_strTitle;
	m_lfMH = stSpec.m_lfMH;
	m_lfMZ = stSpec.m_lfMZ;
	m_lfSqrtMaxInt = stSpec.m_lfSqrtMaxInt;
	m_tCharge = stSpec.m_tCharge;
}

void Spectrum::createHashIndex() {
	m_vHashIndex.clear();
	if (m_vPeaks.empty())
		return;

	m_vHashIndex.reserve(MAX_HASH_SIZE);
	size_t tTemp1 = 0, tTemp2 = m_vPeaks[0].m_nMz / MZ_MULTIPLIER;

	for (size_t i = 0; i <= tTemp2; ++i) {
		m_vHashIndex.push_back(0);
	}

	for (size_t i = 1; i < m_vPeaks.size(); ++i) {
		if ((tTemp1 = size_t(1.0 * m_vPeaks[i].m_nMz / MZ_MULTIPLIER))
				> tTemp2) {
			for (size_t j = tTemp2 + 1; j < tTemp1; ++j) {
				m_vHashIndex.push_back(i - 1);
			}
			m_vHashIndex.push_back(i);
			tTemp2 = tTemp1;
		}
	}
}

void QuerySpectrum::swapOut() {
	vector<QueryPeak>().swap(this->m_vQueryBPeaks);
	vector<QueryPeak>().swap(this->m_vQueryYPeaks);
}

ModificationEntry::ModificationEntry() :
		m_nSite(-1), m_nId(-1) {
}

ModificationEntry::ModificationEntry(int nSite, int nId) :
		m_nSite(nSite), m_nId(nId) {
}

ostream &operator<<(ostream &stOut, const ModificationEntry &stEntry) {
	stOut << stEntry.m_nId << " " << stEntry.m_nSite;
	return stOut;
}

fstream &operator<<(fstream &stOut, const ModificationEntry &stEntry) {
	stOut << stEntry.m_nId << " " << stEntry.m_nSite;
	return stOut;
}

istream &operator>>(istream &stIn, ModificationEntry &stEntry) {
	stIn >> stEntry.m_nId >> stEntry.m_nSite;
	return stIn;
}

fstream &operator>>(fstream &stIn, ModificationEntry &stEntry) {
	stIn >> stEntry.m_nId >> stEntry.m_nSite;
	return stIn;
}

void ModificationEntry::clear() {
	m_nSite = m_nId = -1;
}

std::string ModificationEntry::toDebugString() {
	ostringstream oss;
	oss << "(" << m_nSite << ", " << m_nId << ")";
	return oss.str();
}

PeptideItem::PeptideItem() :
		m_tCount(0), m_tPepId(0) {
}

PeptideItem::PeptideItem(size_t tCount, size_t tPepId) :
		m_tCount(tCount), m_tPepId(tPepId) {
}

bool Peptide::massGreater(const Peptide *pPep1, const Peptide *pPep2) {
	return pPep1->m_lfMass > pPep2->m_lfMass;
}

bool Peptide::massLesser(const Peptide *pPep1, const Peptide *pPep2) {
	return pPep1->m_lfMass < pPep2->m_lfMass;
}

bool Peptide::massTagGreater(const Peptide *pPep1, const Peptide *pPep2) {
	if(pPep1->m_lfMass > pPep2->m_lfMass)
		return true;
	else if(pPep1->m_lfMass < pPep2->m_lfMass)
		return false;
	else
		return pPep1->m_lfTag > pPep2->m_lfTag;
}

bool Peptide::massTagLesser(const Peptide *pPep1, const Peptide *pPep2)
{
	if(pPep1->m_lfMass < pPep2->m_lfMass)
		return true;
	else if(pPep1->m_lfMass > pPep2->m_lfMass)
		return false;
	else
		return pPep1->m_lfTag < pPep2->m_lfTag;

}

Peptide::Peptide() :
		m_lfMass(0.0), m_ucMiss(0), m_ucEnd(0), m_lfTag(0.0) {
	m_strSq.reserve(MAX_PEPTIDE_LEN);
	m_vMods.reserve(MAX_MODIFY_NUM);
}

bool operator==(const Peptide &stPep1, const Peptide &stPep2) {
	if (fabs(stPep1.m_lfMass - stPep2.m_lfMass) > PRECISION)
		return false;

	if (stPep1.m_vMods.size() != stPep2.m_vMods.size())
		return false;

	for (size_t i = 0; i < stPep1.m_vMods.size(); ++i) {
		if (stPep1.m_vMods[i] != stPep2.m_vMods[i])
			return false;
	}

	if (stPep1.m_strSq.compare(stPep2.m_strSq))
		return false;

	return true;
}

bool operator!=(const Peptide &stPep1, const Peptide &stPep2) {
	if (fabs(stPep1.m_lfMass - stPep2.m_lfMass) > PRECISION)
		return true;

	if (stPep1.m_vMods.size() != stPep2.m_vMods.size())
		return true;

	for (size_t i = 0; i < stPep1.m_vMods.size(); ++i) {
		if (stPep1.m_vMods[i] != stPep2.m_vMods[i])
			return true;
	}

	if (stPep1.m_strSq.compare(stPep2.m_strSq))
		return true;

	return false;
}

ostream &operator<<(ostream &stOut, const Peptide &stPeptide) {
	stOut << fixed << setprecision(FIXED_PRECISION);
	if(0 == stPeptide.m_strSq.size())
		stOut << NON_AA << " ";
	else
		stOut << stPeptide.m_strSq << " ";
	stOut << stPeptide.m_lfMass << " " << stPeptide.m_vMods.size() << " ";
	for (size_t i = 0; i < stPeptide.m_vMods.size(); ++i) {
		stOut << stPeptide.m_vMods[i].m_nId << " "
				<< stPeptide.m_vMods[i].m_nSite << " ";
	}
	stOut << (int) stPeptide.m_ucEnd << " " << (int) stPeptide.m_ucMiss;
	return stOut;
}

fstream &operator<<(std::fstream &stOut, const Peptide &stPeptide) {
	stOut << fixed << setprecision(FIXED_PRECISION);
	if(0 == stPeptide.m_strSq.size())
		stOut << NON_AA << " ";
	else
		stOut << stPeptide.m_strSq << " ";
	stOut << stPeptide.m_lfMass << " " << stPeptide.m_vMods.size() << " ";
	for (size_t i = 0; i < stPeptide.m_vMods.size(); ++i) {
		stOut << stPeptide.m_vMods[i] << " ";
	}
	stOut << (int) stPeptide.m_ucEnd << " " << (int) stPeptide.m_ucMiss;
	return stOut;
}

istream &operator>>(istream &stIn, Peptide &stPeptide) {
	stIn >> stPeptide.m_strSq >> stPeptide.m_lfMass;
	if(0 == stPeptide.m_strSq.compare(NON_AA))
	{
		stPeptide.m_strSq.clear();
	}
	size_t tModSize(0);
	stIn >> tModSize;
	for (size_t i = 0; i < tModSize; ++i) {
		ModificationEntry stEntry;
		stIn >> stEntry;
		stPeptide.m_vMods.push_back(stEntry);
	}
	size_t tTmp(0);
	stIn >> tTmp;
	stPeptide.m_ucEnd = (unsigned char) tTmp;
	stIn >> tTmp;
	stPeptide.m_ucMiss = (unsigned char) tTmp;
	return stIn;
}

fstream &operator>>(fstream &stIn, Peptide &stPeptide) {
	stIn >> stPeptide.m_strSq >> stPeptide.m_lfMass;
	if(0 == stPeptide.m_strSq.compare(NON_AA))
	{
		stPeptide.m_strSq.clear();
	}
	size_t tModSize(0);
	stIn >> tModSize;
	for (size_t i = 0; i < tModSize; ++i) {
		ModificationEntry stEntry;
		stIn >> stEntry;
		stPeptide.m_vMods.push_back(stEntry);
	}
	size_t tTmp(0);
	stIn >> tTmp;
	stPeptide.m_ucEnd = (unsigned char) tTmp;
	stIn >> tTmp;
	stPeptide.m_ucMiss = (unsigned char) tTmp;
	return stIn;
}

bool Peptide::isProNTerm() const {
	return (m_ucEnd == 1) || (m_ucEnd == 3);
}

bool Peptide::isProCTerm() const {
	return (m_ucEnd == 2) || (m_ucEnd == 3);
}

void Peptide::setPeptideInfo(const char *szSq, double lfMass,
		unsigned char ucMiss, unsigned char ucEnd) {
	m_strSq = szSq;
	m_lfMass = lfMass;
	m_ucMiss = ucMiss;
	m_ucEnd = ucEnd;
}

void Peptide::setPeptideInfo(const std::string &strSq, double lfMass,
		unsigned char ucMiss, unsigned char ucEnd) {
	setPeptideInfo(strSq.c_str(), lfMass, ucMiss, ucEnd);
}

void Peptide::clear() {
	m_strSq.clear();
	m_lfMass = 0.0;
	m_ucMiss = 0;
	m_ucEnd = 0;
	m_vMods.clear();
}

string Peptide::toString() {
	ostringstream oss;
	oss << m_strSq << " " << m_lfMass << " " << m_vMods.size() << " ";
	for (size_t i = 0; i < m_vMods.size(); ++i) {
		oss << m_vMods[i].m_nId << " " << m_vMods[i].m_nSite << " ";
	}
	oss << (int) m_ucEnd << " " << (int) m_ucMiss;
	return oss.str();
}

string Peptide::toDebugString() {
	ostringstream oss;
	oss << "sequence=" << m_strSq << "\n" << "mass=" << m_lfMass << "\n"
			<< "modifications=" << m_vMods.size();
	if (m_vMods.size() > 0) {
		for (size_t i = 0; i < m_vMods.size(); ++i) {
			oss << m_vMods[i].toDebugString();
		}
	}
	oss << "\n";
	oss << "peptide_type=" << (int) m_ucEnd << "\n" << "miss_site="
			<< (int) m_ucMiss;
	return oss.str();
}

double Peptide::getGodelCode() const {
	ConstantStore *pConst = ConstantStore::getInstance();
	double lfCode = 0.0;
	size_t tIdx = 0;
	for (; tIdx < m_strSq.size(); ++tIdx) {
		lfCode += (m_strSq[tIdx] - 'A' + 1) * pConst->m_lfLogPrime[tIdx];
	}
	lfCode += m_ucEnd * pConst->m_lfLogPrime[tIdx];
	return lfCode;
}

double Peptide::calculatePeptideMass()
{
	ConstantStore *pConst = ConstantStore::getInstance();
	m_lfMass = 0.0;
	for(size_t tIdx = 0; tIdx < m_strSq.size(); ++tIdx) {
		m_lfMass += pConst->m_lfAAMass[m_strSq[tIdx]-'A'];
	}

	SearchParameter *pParameter = ParameterReader::getInstance()->getParameter();
	int nVarSize = pParameter->m_vVariableMods.size();
	ModificationDict *pModDict = ModificationDict::getInstance();
	for(size_t tModIdx = 0; tModIdx < m_vMods.size(); ++tModIdx) {
		string strMod;
		if(m_vMods[tModIdx].m_nId >= nVarSize) {
			strMod = pParameter->m_vFixedMods[m_vMods[tModIdx].m_nId-nVarSize];
		} else {
			strMod = pParameter->m_vVariableMods[m_vMods[tModIdx].m_nId];
		}
		m_lfMass += pModDict->getMonoDiff(strMod);
		vector<double> vLoss = pModDict->getMonoNeutralLossDiff(strMod);
		for(size_t tLossIdx = 0; tLossIdx < vLoss.size(); ++tLossIdx) {
			m_lfMass += vLoss[tLossIdx];
		}
	}
	return m_lfMass;
}

ostream &operator<<(ostream &stOut, const ProteinInformation &stInfo) {
	stOut << stInfo.m_tProteinID.size() << " ";
	for (size_t tIdx = 0; tIdx < stInfo.m_tProteinID.size(); ++tIdx) {
		stOut << stInfo.m_tProteinID[tIdx] << " ";
	}
	stOut << stInfo.m_tProteinSite.size() << " ";
	for (size_t tIdx = 0; tIdx < stInfo.m_tProteinSite.size(); ++tIdx) {
		stOut << stInfo.m_tProteinSite[tIdx] << " ";
	}
	return stOut;
}

fstream &operator<<(fstream &stOut, const ProteinInformation &stInfo) {
	stOut << stInfo.m_tProteinID.size() << " ";
	for (size_t tIdx = 0; tIdx < stInfo.m_tProteinID.size(); ++tIdx) {
		stOut << stInfo.m_tProteinID[tIdx] << " ";
	}
	stOut << stInfo.m_tProteinSite.size() << " ";
	for (size_t tIdx = 0; tIdx < stInfo.m_tProteinSite.size(); ++tIdx) {
		stOut << stInfo.m_tProteinSite[tIdx] << " ";
	}
	return stOut;
}

istream &operator>>(istream &stIn, ProteinInformation &stInfo) {
	size_t tSize(0);
	stIn >> tSize;
	for (size_t tIdx = 0; tIdx < tSize; ++tIdx) {
		size_t tID(0);
		stIn >> tID;
		stInfo.m_tProteinID.push_back(tID);
	}
	stIn >> tSize;
	for (size_t tIdx = 0; tIdx < tSize; ++tIdx) {
		size_t tSite(0);
		stIn >> tSite;
		stInfo.m_tProteinSite.push_back(tSite);
	}
	return stIn;
}

fstream &operator>>(fstream &stIn, ProteinInformation &stInfo) {
	size_t tSize(0);
	stIn >> tSize;
	for (size_t tIdx = 0; tIdx < tSize; ++tIdx) {
		size_t tID(0);
		stIn >> tID;
		stInfo.m_tProteinID.push_back(tID);
	}
	stIn >> tSize;
	for (size_t tIdx = 0; tIdx < tSize; ++tIdx) {
		size_t tSite(0);
		stIn >> tSite;
		stInfo.m_tProteinSite.push_back(tSite);
	}
	return stIn;
}

string ProteinInformation::toString() {
	ostringstream oss;
	oss << m_tProteinID.size() << " ";
	for (size_t tIdx = 0; tIdx < m_tProteinID.size(); ++tIdx) {
		oss << m_tProteinID[tIdx] << " ";
	}
	oss << m_tProteinSite.size() << " ";
	for (size_t tIdx = 0; tIdx < m_tProteinSite.size(); ++tIdx) {
		oss << m_tProteinSite[tIdx] << " ";
	}
	return oss.str();
}

bool PeptideResult::scoreGreater(const PeptideResult &stResult1,
		const PeptideResult &stResult2) {
	return stResult1.m_lfScore > stResult2.m_lfScore;
}

bool PeptideResult::scoreLesser(const PeptideResult &stResult1,
		const PeptideResult &stResult2) {
	return stResult1.m_lfScore < stResult2.m_lfScore;
}

PeptideResult::PeptideResult() :
		m_tWinId(0), m_lfOpenMass(0.0), m_lfScore(0.0), m_pProteinInfo(NULL) {
}

PeptideResult::PeptideResult(const PeptideResult &stResult) :
		m_stPep(stResult.m_stPep), m_tWinId(stResult.m_tWinId), m_lfOpenMass(
				stResult.m_lfOpenMass), m_lfScore(stResult.m_lfScore), m_stLinkMod(
				stResult.m_stLinkMod), m_pProteinInfo(NULL) {
	if (stResult.m_pProteinInfo != NULL)
		m_pProteinInfo = new ProteinInformation(*stResult.m_pProteinInfo);
}

PeptideResult::~PeptideResult() {
	if (m_pProteinInfo != NULL) {
		delete m_pProteinInfo;
		m_pProteinInfo = NULL;
	}
}

PeptideResult &PeptideResult::operator=(const PeptideResult &stResult) {
	m_stPep = stResult.m_stPep;
	m_tWinId = stResult.m_tWinId;
	m_lfOpenMass = stResult.m_lfOpenMass;
	m_lfScore = stResult.m_lfScore;
	m_stLinkMod = stResult.m_stLinkMod;

	ProteinInformation *pInfo = m_pProteinInfo;
	if (stResult.m_pProteinInfo != NULL)
		m_pProteinInfo = new ProteinInformation(*stResult.m_pProteinInfo);
	else
		m_pProteinInfo = NULL;
	if (pInfo != NULL) {
		delete pInfo;
	}

	return *this;
}

bool operator==(const PeptideResult &stResult1,
		const PeptideResult &stResult2) {
	if (stResult1.m_stLinkMod != stResult2.m_stLinkMod)
		return false;

	if (stResult1.m_stPep != stResult2.m_stPep)
		return false;

	return true;
}

bool operator!=(const PeptideResult &stResult1,
		const PeptideResult &stResult2) {
	if (stResult1.m_stLinkMod != stResult2.m_stLinkMod)
		return true;

	if (stResult1.m_stPep != stResult2.m_stPep)
		return true;

	return false;
}

ostream &operator<<(ostream &stOut, const PeptideResult &stResult) {
	stOut << fixed << setprecision(FIXED_PRECISION);
	stOut << stResult.m_stPep << " " << stResult.m_lfOpenMass << " "
			<< stResult.m_lfScore << " " << stResult.m_stLinkMod << " "
			<< stResult.m_tWinId << " ";

	if (stResult.m_pProteinInfo == NULL) {
		stOut << "0";
	} else {
		stOut << "1" << " " << *stResult.m_pProteinInfo;
	}
	return stOut;
}

fstream &operator<<(fstream &stOut, const PeptideResult &stResult) {
	stOut << fixed << setprecision(FIXED_PRECISION);
	stOut << stResult.m_stPep << " " << stResult.m_lfOpenMass << " "
			<< stResult.m_lfScore << " " << stResult.m_stLinkMod << " "
			<< stResult.m_tWinId << " ";

	if (stResult.m_pProteinInfo == NULL) {
		stOut << "0";
	} else {
		stOut << "1" << " " << *stResult.m_pProteinInfo;
	}
	return stOut;
}

istream &operator>>(istream &stIn, PeptideResult &stResult) {
	stIn >> stResult.m_stPep >> stResult.m_lfOpenMass >> stResult.m_lfScore
			>> stResult.m_stLinkMod >> stResult.m_tWinId;

	int nFlag(0);
	stIn >> nFlag;
	if (nFlag == 1) {
		stResult.m_pProteinInfo = new ProteinInformation;
		stIn >> *stResult.m_pProteinInfo;
	}
	return stIn;
}

fstream &operator>>(fstream &stIn, PeptideResult &stResult) {
	stIn >> stResult.m_stPep >> stResult.m_lfOpenMass >> stResult.m_lfScore
			>> stResult.m_stLinkMod >> stResult.m_tWinId;

	int nFlag(0);
	stIn >> nFlag;
	if (nFlag == 1) {
		stResult.m_pProteinInfo = new ProteinInformation;
		stIn >> *stResult.m_pProteinInfo;
	}
	return stIn;
}

bool PeptideResult::isTarget(const char *szDecoySign) const {
	string strFastaPath = ParameterReader::getInstance()->getParameter()->m_strFastaFilePath;
	AccessionIndexer *pACIndexer = AccessionIndexer::getAccessionIndex(strFastaPath);
	const vector<size_t> &vID = m_pProteinInfo->m_tProteinID;
	for (size_t t = 0; t < vID.size(); ++t) {
		if (pACIndexer->getACByID(vID[t]).substr(0, strlen(szDecoySign)).compare(szDecoySign)) {
			return true;
		}
	}
	return false;
}

void PeptideResult::clear() {
	m_stPep.clear();
	m_tWinId = 0;
	m_lfOpenMass = 0.0;
	m_lfScore = 0.0;
	m_stLinkMod.clear();
}

string PeptideResult::toString() {
	ostringstream oss;
	oss << m_stPep << " " << m_lfOpenMass << " " << m_lfScore << " "
			<< m_stLinkMod.m_nId << " " << m_stLinkMod.m_nSite << " ";
	if (m_pProteinInfo != NULL) {
		oss << m_tWinId << " " << *m_pProteinInfo;
	} else {
		oss << m_tWinId;
	}
	return oss.str();
}

string PeptideResult::toDebugString() {
	ostringstream oss;
	oss << m_stPep.toDebugString() << "\n" << "open_mass=" << m_lfOpenMass
			<< "\n" << "linker_site=" << m_stLinkMod.m_nSite << "\n"
			<< "linker_id=" << m_stLinkMod.m_nId;
	return oss.str();
}

ostream &operator<<(ostream &stOut, const RatioFeature &stInfo) {
	stOut << fixed << setprecision(FIXED_PRECISION);
	stOut << stInfo.m_nLength << " " << stInfo.m_nTagLen << " "
			<< stInfo.m_lfMatchedPeakRatio << " " << stInfo.m_lfMatchedIonRatio;
	return stOut;
}

fstream &operator<<(fstream &stOut, const RatioFeature &stInfo) {
	stOut << fixed << setprecision(FIXED_PRECISION);
	stOut << stInfo.m_nLength << " " << stInfo.m_nTagLen << " "
			<< stInfo.m_lfMatchedPeakRatio << " " << stInfo.m_lfMatchedIonRatio;
	return stOut;
}

istream &operator>>(istream &stIn, RatioFeature &stInfo) {
	stIn >> stInfo.m_nLength >> stInfo.m_nTagLen >> stInfo.m_lfMatchedPeakRatio
			>> stInfo.m_lfMatchedIonRatio;
	return stIn;
}

fstream &operator>>(fstream &stIn, RatioFeature &stInfo) {
	stIn >> stInfo.m_nLength >> stInfo.m_nTagLen >> stInfo.m_lfMatchedPeakRatio
			>> stInfo.m_lfMatchedIonRatio;
	return stIn;
}

void RatioFeature::clear() {
	m_nLength = m_nTagLen = 0;
	m_lfMatchedPeakRatio = m_lfMatchedIonRatio = 0.0;
}

ostream &operator<<(ostream &stOut, const XLinkFeature &stInfo) {
	stOut << fixed << setprecision(FIXED_PRECISION);
	stOut << stInfo.m_stAlpha << " " << stInfo.m_stBeta;
	return stOut;
}

fstream &operator<<(fstream &stOut, const XLinkFeature &stInfo) {
	stOut << fixed << setprecision(FIXED_PRECISION);
	stOut << stInfo.m_stAlpha << " " << stInfo.m_stBeta;
	return stOut;
}

istream &operator>>(istream &stIn, XLinkFeature &stInfo) {
	stIn >> stInfo.m_stAlpha >> stInfo.m_stBeta;
	return stIn;
}

fstream &operator>>(fstream &stIn, XLinkFeature &stInfo) {
	stIn >> stInfo.m_stAlpha >> stInfo.m_stBeta;
	return stIn;
}

void XLinkFeature::clear() {
	m_stAlpha.clear();
	m_stBeta.clear();
}

FeatureInformation::FeatureInformation() :
		m_lfErrorMean(0.0), m_lfErrorDeviation(0.0) {
}

ostream &operator<<(ostream &stOut, const FeatureInformation &stInfo) {
	stOut << fixed << setprecision(FIXED_PRECISION);
	stOut << stInfo.m_stRatio << " " << stInfo.m_lfErrorMean << " "
			<< stInfo.m_lfErrorDeviation << " " << stInfo.m_stXLink;
	return stOut;
}

fstream &operator<<(fstream &stOut, const FeatureInformation &stInfo) {
	stOut << fixed << setprecision(FIXED_PRECISION);
	stOut << stInfo.m_stRatio << " " << stInfo.m_lfErrorMean << " "
			<< stInfo.m_lfErrorDeviation << " " << stInfo.m_stXLink;
	return stOut;
}
istream &operator>>(istream &stIn, FeatureInformation &stInfo) {
	stIn >> stInfo.m_stRatio >> stInfo.m_lfErrorMean
			>> stInfo.m_lfErrorDeviation >> stInfo.m_stXLink;
	return stIn;
}

fstream &operator>>(fstream &stIn, FeatureInformation &stInfo) {
	stIn >> stInfo.m_stRatio >> stInfo.m_lfErrorMean
			>> stInfo.m_lfErrorDeviation >> stInfo.m_stXLink;
	return stIn;
}

string FeatureInformation::toString() {
	return "";
}

string FeatureInformation::toDebugString() {
	return "";
}

void FeatureInformation::clear() {
	m_stRatio.clear();
	m_lfErrorMean = m_lfErrorDeviation = 0.0;
	m_stXLink.clear();
}

ReportInformation::ReportInformation() :
		m_lfMatchedPercentage(0.0), m_lfTotalIntensity(0.0), m_tMatchedNumber(
				0), m_tTotalNumber(0) {
}

ostream &operator<<(ostream &stOut, const ReportInformation &stInfo) {
	stOut << fixed << setprecision(FIXED_PRECISION);
	stOut << stInfo.m_lfMatchedPercentage << " " << stInfo.m_lfTotalIntensity
			<< " " << stInfo.m_tMatchedNumber << " " << stInfo.m_tTotalNumber
			<< " ";
	stOut << stInfo.m_bNContinuity.size() << " ";
	for (size_t tIdx = 0; tIdx < stInfo.m_bNContinuity.size(); ++tIdx) {
		stOut << stInfo.m_bNContinuity[tIdx] << " ";
	}
	stOut << stInfo.m_bCContinuity.size() << " ";
	for (size_t tIdx = 0; tIdx < stInfo.m_bCContinuity.size(); ++tIdx) {
		stOut << stInfo.m_bCContinuity[tIdx] << " ";
	}
	stOut << stInfo.m_nNCoexistence.size() << " ";
	for (size_t tIdx = 0; tIdx < stInfo.m_nNCoexistence.size(); ++tIdx) {
		stOut << stInfo.m_nNCoexistence[tIdx] << " ";
	}
	stOut << stInfo.m_nCCoexistence.size() << " ";
	for (size_t tIdx = 0; tIdx < stInfo.m_nCCoexistence.size(); ++tIdx) {
		stOut << stInfo.m_nCCoexistence[tIdx] << " ";
	}
	stOut << stInfo.m_vIonNumber.size() << " ";
	for (size_t tIdx = 0; tIdx < stInfo.m_vIonNumber.size(); ++tIdx) {
		stOut << stInfo.m_vIonNumber[tIdx] << " ";
	}
	if (!stInfo.m_vIonMatchedNumber.empty())
		stOut << stInfo.m_vIonMatchedNumber.size() << " ";
	else
		stOut << stInfo.m_vIonMatchedNumber.size();
	for (size_t tIdx = 0; tIdx < stInfo.m_vIonMatchedNumber.size(); ++tIdx) {
		stOut << stInfo.m_vIonMatchedNumber[tIdx];
		if (tIdx + 1 < stInfo.m_vIonMatchedNumber.size())
			stOut << " ";
	}
	return stOut;
}

fstream &operator<<(fstream &stOut, const ReportInformation &stInfo) {
	stOut << fixed << setprecision(FIXED_PRECISION);
	stOut << stInfo.m_lfMatchedPercentage << " " << stInfo.m_lfTotalIntensity
			<< " " << stInfo.m_tMatchedNumber << " " << stInfo.m_tTotalNumber
			<< " ";
	stOut << stInfo.m_bNContinuity.size() << " ";
	for (size_t tIdx = 0; tIdx < stInfo.m_bNContinuity.size(); ++tIdx) {
		stOut << stInfo.m_bNContinuity[tIdx] << " ";
	}
	stOut << stInfo.m_bCContinuity.size() << " ";
	for (size_t tIdx = 0; tIdx < stInfo.m_bCContinuity.size(); ++tIdx) {
		stOut << stInfo.m_bCContinuity[tIdx] << " ";
	}
	stOut << stInfo.m_nNCoexistence.size() << " ";
	for (size_t tIdx = 0; tIdx < stInfo.m_nNCoexistence.size(); ++tIdx) {
		stOut << stInfo.m_nNCoexistence[tIdx] << " ";
	}
	stOut << stInfo.m_nCCoexistence.size() << " ";
	for (size_t tIdx = 0; tIdx < stInfo.m_nCCoexistence.size(); ++tIdx) {
		stOut << stInfo.m_nCCoexistence[tIdx] << " ";
	}
	stOut << stInfo.m_vIonNumber.size() << " ";
	for (size_t tIdx = 0; tIdx < stInfo.m_vIonNumber.size(); ++tIdx) {
		stOut << stInfo.m_vIonNumber[tIdx] << " ";
	}
	if (!stInfo.m_vIonMatchedNumber.empty())
		stOut << stInfo.m_vIonMatchedNumber.size() << " ";
	else
		stOut << stInfo.m_vIonMatchedNumber.size();
	for (size_t tIdx = 0; tIdx < stInfo.m_vIonMatchedNumber.size(); ++tIdx) {
		stOut << stInfo.m_vIonMatchedNumber[tIdx];
		if (tIdx + 1 < stInfo.m_vIonMatchedNumber.size())
			stOut << " ";
	}
	return stOut;
}

istream &operator>>(istream &stIn, ReportInformation &stInfo) {
	stIn >> stInfo.m_lfMatchedPercentage >> stInfo.m_lfTotalIntensity
			>> stInfo.m_tMatchedNumber >> stInfo.m_tTotalNumber;
	size_t tSize(0);
	stIn >> tSize;
	stInfo.m_bNContinuity.resize(tSize, false);
	for (size_t tIdx = 0; tIdx < tSize; ++tIdx) {
		bool bContinuity(false);
		stIn >> bContinuity;
		stInfo.m_bNContinuity[tIdx] = bContinuity;
	}
	stIn >> tSize;
	stInfo.m_bCContinuity.resize(tSize, false);
	for (size_t tIdx = 0; tIdx < tSize; ++tIdx) {
		bool bContinuity(false);
		stIn >> bContinuity;
		stInfo.m_bCContinuity[tIdx] = bContinuity;
	}
	stIn >> tSize;
	stInfo.m_nNCoexistence.resize(tSize, false);
	for (size_t tIdx = 0; tIdx < tSize; ++tIdx) {
		stIn >> stInfo.m_nNCoexistence[tIdx];
	}
	stIn >> tSize;
	stInfo.m_nCCoexistence.resize(tSize, false);
	for (size_t tIdx = 0; tIdx < tSize; ++tIdx) {
		stIn >> stInfo.m_nCCoexistence[tIdx];
	}
	stIn >> tSize;
	stInfo.m_vIonNumber.resize(tSize, false);
	for (size_t tIdx = 0; tIdx < tSize; ++tIdx) {
		stIn >> stInfo.m_vIonNumber[tIdx];
	}
	stIn >> tSize;
	stInfo.m_vIonMatchedNumber.resize(tSize, false);
	for (size_t tIdx = 0; tIdx < tSize; ++tIdx) {
		stIn >> stInfo.m_vIonMatchedNumber[tIdx];
	}
	return stIn;
}

fstream &operator>>(fstream &stIn, ReportInformation &stInfo) {
	stIn >> stInfo.m_lfMatchedPercentage >> stInfo.m_lfTotalIntensity
			>> stInfo.m_tMatchedNumber >> stInfo.m_tTotalNumber;
	size_t tSize(0);
	stIn >> tSize;
	stInfo.m_bNContinuity.resize(tSize, false);
	for (size_t tIdx = 0; tIdx < tSize; ++tIdx) {
		bool bContinuity(false);
		stIn >> bContinuity;
		stInfo.m_bNContinuity[tIdx] = bContinuity;
	}
	stIn >> tSize;
	stInfo.m_bCContinuity.resize(tSize, false);
	for (size_t tIdx = 0; tIdx < tSize; ++tIdx) {
		bool bContinuity(false);
		stIn >> bContinuity;
		stInfo.m_bCContinuity[tIdx] = bContinuity;
	}
	stIn >> tSize;
	stInfo.m_nNCoexistence.resize(tSize, false);
	for (size_t tIdx = 0; tIdx < tSize; ++tIdx) {
		stIn >> stInfo.m_nNCoexistence[tIdx];
	}
	stIn >> tSize;
	stInfo.m_nCCoexistence.resize(tSize, false);
	for (size_t tIdx = 0; tIdx < tSize; ++tIdx) {
		stIn >> stInfo.m_nCCoexistence[tIdx];
	}
	stIn >> tSize;
	stInfo.m_vIonNumber.resize(tSize, false);
	for (size_t tIdx = 0; tIdx < tSize; ++tIdx) {
		stIn >> stInfo.m_vIonNumber[tIdx];
	}
	stIn >> tSize;
	stInfo.m_vIonMatchedNumber.resize(tSize, false);
	for (size_t tIdx = 0; tIdx < tSize; ++tIdx) {
		stIn >> stInfo.m_vIonMatchedNumber[tIdx];
	}
	return stIn;
}

string ReportInformation::toString() {
	stringstream oss;
	oss << m_lfMatchedPercentage << " " << m_lfTotalIntensity << " "
			<< m_tMatchedNumber << " " << m_tTotalNumber << " ";
	oss << m_bNContinuity.size() << " ";
	for (size_t tIdx = 0; tIdx < m_bNContinuity.size(); ++tIdx) {
		oss << m_bNContinuity[tIdx] << " ";
	}
	oss << m_bCContinuity.size() << " ";
	for (size_t tIdx = 0; tIdx < m_bCContinuity.size(); ++tIdx) {
		oss << m_bCContinuity[tIdx] << " ";
	}
	oss << m_nNCoexistence.size() << " ";
	for (size_t tIdx = 0; tIdx < m_nNCoexistence.size(); ++tIdx) {
		oss << m_nNCoexistence[tIdx] << " ";
	}
	oss << m_nCCoexistence.size() << " ";
	for (size_t tIdx = 0; tIdx < m_nCCoexistence.size(); ++tIdx) {
		oss << m_nCCoexistence[tIdx] << " ";
	}
	oss << m_vIonNumber.size() << " ";
	for (size_t tIdx = 0; tIdx < m_vIonNumber.size(); ++tIdx) {
		oss << m_vIonNumber[tIdx] << " ";
	}
	if (!m_vIonMatchedNumber.empty())
		oss << m_vIonMatchedNumber.size() << " ";
	else
		oss << m_vIonMatchedNumber.size();
	for (size_t tIdx = 0; tIdx < m_vIonMatchedNumber.size(); ++tIdx) {
		oss << m_vIonMatchedNumber[tIdx];
		if (tIdx + 1 < m_vIonMatchedNumber.size())
			oss << " ";
	}
	return oss.str();
}

string ReportInformation::toDebugString() {
	ostringstream oss;
	oss << "matched_intensity_percentage=" << m_lfMatchedPercentage << "\n"
			<< "total_intensity=" << m_lfTotalIntensity << "\n"
			<< "matched_number=" << m_tMatchedNumber << "\n" << "total_number="
			<< m_tTotalNumber << "\n";
	oss << "n_continuity=";
	for (size_t tIdx = 0; tIdx < m_bNContinuity.size(); ++tIdx) {
		oss << m_bNContinuity[tIdx] << " ";
	}
	oss << "\nc_continuity=";
	for (size_t tIdx = 0; tIdx < m_bCContinuity.size(); ++tIdx) {
		oss << m_bCContinuity[tIdx] << " ";
	}
	oss << "\nn_coexistence=";
	for (size_t tIdx = 0; tIdx < m_nNCoexistence.size(); ++tIdx) {
		oss << m_nNCoexistence[tIdx] << " ";
	}
	oss << "\nc_coexistence=";
	for (size_t tIdx = 0; tIdx < m_nCCoexistence.size(); ++tIdx) {
		oss << m_nCCoexistence[tIdx] << " ";
	}
	oss << "\nion_number=";
	for (size_t tIdx = 0; tIdx < m_vIonNumber.size(); ++tIdx) {
		oss << m_vIonNumber[tIdx] << " ";
	}
	oss << "\nmatched_number=";
	for (size_t tIdx = 0; tIdx < m_vIonMatchedNumber.size(); ++tIdx) {
		oss << m_vIonMatchedNumber[tIdx] << " ";
	}
	return oss.str();
}

XLinkPeptideItem::XLinkPeptideItem(string strSq, size_t tSpecIdx,
		size_t tRankIdx, bool bAlpha) :
		m_strSq(strSq), m_tSpecIdx(tSpecIdx), m_tRankIdx(tRankIdx), m_bAlpha(
				bAlpha) {

}

XLinkPeptideResult::XLinkPeptideResult() :
		m_ePepType(PT_COMMON), m_eProType(PRT_NONE), m_eTDType(TDT_F), m_lfScore(
				0.0), m_lfEvalue(0.0), m_lfQvalue(0.0), m_lfCalcMH(0.0) {
}

XLinkPeptideResult::~XLinkPeptideResult() {
}

ostream &operator<<(ostream &stOut, const XLinkPeptideResult &stResult) {
	stOut.setf(ios::fixed, ios::floatfield);
	stOut.precision(FIXED_PRECISION);
	stOut << stResult.m_ePepType << " " << stResult.m_eProType << " "
			<< stResult.m_eTDType << " " << stResult.m_lfCalcMH << " "
			<< stResult.m_lfScore << " ";

	stOut.setf(ios::scientific, ios::floatfield);
	stOut.precision(FIXED_PRECISION);
	stOut << stResult.m_lfEvalue << " " << stResult.m_lfQvalue << " ";

	stOut.setf(ios::fixed, ios::floatfield);
	stOut.precision(FIXED_PRECISION);
	stOut << stResult.m_stAlphaPep << " " << stResult.m_stBetaPep << " "
			<< stResult.m_stFeatureInfo;
	return stOut;
}

fstream &operator<<(fstream &stOut, const XLinkPeptideResult &stResult) {
	stOut.setf(ios::fixed, ios::floatfield);
	stOut.precision(FIXED_PRECISION);
	stOut << stResult.m_ePepType << " " << stResult.m_eProType << " "
			<< stResult.m_eTDType << " " << stResult.m_lfCalcMH << " "
			<< stResult.m_lfScore << " ";

	stOut.setf(ios::scientific, ios::floatfield);
	stOut.precision(FIXED_PRECISION);
	stOut << stResult.m_lfEvalue << " " << stResult.m_lfQvalue << " ";

	stOut.setf(ios::fixed, ios::floatfield);
	stOut.precision(FIXED_PRECISION);
	stOut << stResult.m_stAlphaPep << " " << stResult.m_stBetaPep << " "
			<< stResult.m_stFeatureInfo;
	return stOut;
}

istream &operator>>(istream &stIn, XLinkPeptideResult &stResult) {
	int nType(0);
	stIn >> nType;
	stResult.m_ePepType = (PeptideType) nType;
	stIn >> nType;
	stResult.m_eProType = (ProteinType) nType;
	stIn >> nType;
	stResult.m_eTDType = (TDType) nType;
	stIn >> stResult.m_lfCalcMH >> stResult.m_lfScore >> stResult.m_lfEvalue
			>> stResult.m_lfQvalue >> stResult.m_stAlphaPep
			>> stResult.m_stBetaPep >> stResult.m_stFeatureInfo;
	return stIn;
}

fstream &operator>>(fstream &stIn, XLinkPeptideResult &stResult) {
	int nType(0);
	stIn >> nType;
	stResult.m_ePepType = (PeptideType) nType;
	stIn >> nType;
	stResult.m_eProType = (ProteinType) nType;
	stIn >> nType;
	stResult.m_eTDType = (TDType) nType;
	stIn >> stResult.m_lfCalcMH >> stResult.m_lfScore >> stResult.m_lfEvalue
			>> stResult.m_lfQvalue >> stResult.m_stAlphaPep
			>> stResult.m_stBetaPep >> stResult.m_stFeatureInfo;
	return stIn;
}

bool operator<(const XLinkPeptideResult &stResult1,
		const XLinkPeptideResult &stResult2) {
	return stResult1.m_lfScore < stResult2.m_lfScore;
}

bool operator>(const XLinkPeptideResult &stResult1,
		const XLinkPeptideResult &stResult2) {
	return stResult1.m_lfScore > stResult2.m_lfScore;
}

bool operator==(const XLinkPeptideResult &stResult1,
		const XLinkPeptideResult &stResult2) {
	if (fabs(stResult1.m_lfScore - stResult2.m_lfScore) > PRECISION)
		return false;

	if (stResult1.m_ePepType != stResult2.m_ePepType)
		return false;

	if (fabs(stResult1.m_lfCalcMH - stResult2.m_lfCalcMH) > PRECISION)
		return false;

	if (stResult1.m_stAlphaPep != stResult2.m_stAlphaPep)
		return false;

	if (stResult1.m_stBetaPep != stResult2.m_stBetaPep)
		return false;

	return true;
}

bool operator!=(const XLinkPeptideResult &stResult1,
		const XLinkPeptideResult &stResult2) {
	if (fabs(stResult1.m_lfScore - stResult2.m_lfScore) > PRECISION)
		return true;

	if (stResult1.m_ePepType != stResult2.m_ePepType)
		return true;

	if (fabs(stResult1.m_lfCalcMH - stResult2.m_lfCalcMH) > PRECISION)
		return true;

	if (stResult1.m_stAlphaPep != stResult2.m_stAlphaPep)
		return true;

	if (stResult1.m_stBetaPep != stResult2.m_stBetaPep)
		return true;

	return false;
}

void XLinkPeptideResult::setProteinTDType(const char *szDecoySign)
{
	// initialize protein type, target and decoy type
	if(m_ePepType == PT_XLINK) {
		int nAlpha = 0, nBeta = 0;
		if(m_stAlphaPep.m_pProteinInfo == NULL) {
			ErrorInfo err("XLinkPeptideResult", "setProteinTDType", "no protein information.");
			throw runtime_error(err.get());
		}
		if(m_stAlphaPep.isTarget(szDecoySign))
			nAlpha = 1;

		if(m_stBetaPep.m_pProteinInfo == NULL) {
			ErrorInfo err("XLinkPeptideResult", "setProteinTDType", "no protein information, xlink.");
			throw runtime_error(err.get());
		}

		if(m_stBetaPep.isTarget(szDecoySign))
			nBeta = 1;

		if(nAlpha+nBeta == 2) {
			m_eTDType = TDT_T;
		} else if(nAlpha+nBeta == 1) {
			m_eTDType = TDT_U;
		} else {
			m_eTDType = TDT_F;
		}

		m_eProType = getProteinType(szDecoySign);
	} else {
		if(m_stAlphaPep.m_pProteinInfo == NULL) {
			ErrorInfo err("XLinkPeptideResult", "setProteinTDType", "no protein information, non-xlink.");
			throw runtime_error(err.get());
		}

		if(m_stAlphaPep.isTarget(szDecoySign)) {
			m_eTDType = TDT_T;
		} else {
			m_eTDType = TDT_F;
		}

		m_eProType = PRT_NONE;
	}
}

ProteinType XLinkPeptideResult::getProteinType(const char *szDecoySign) const
{
	if(m_ePepType != PT_XLINK)
		return PRT_NONE;
	string strFastaPath = ParameterReader::getInstance()->getParameter()->m_strFastaFilePath;
	AccessionIndexer *pACIndexer = AccessionIndexer::getAccessionIndex(strFastaPath);

	vector<size_t> vAlphaIDs, vBetaIDs;


	const vector<size_t> &vID = m_stAlphaPep.m_pProteinInfo->m_tProteinID;
	for(size_t i = 0; i < vID.size(); ++i) {
		if(vID[i]%2 == 0)
			vAlphaIDs.push_back(vID[i]);
		else
			vAlphaIDs.push_back(vID[i]-1);
	}

	const vector<size_t> &vID2 = m_stBetaPep.m_pProteinInfo->m_tProteinID;
	for(size_t i = 0; i < vID2.size(); ++i) {
		if(vID[i]%2 == 0)
			vBetaIDs.push_back(vID2[i]);
		else
			vBetaIDs.push_back(vID2[i]-1);
	}

	sort(vAlphaIDs.begin(), vAlphaIDs.end());
	sort(vBetaIDs.begin(), vBetaIDs.end());
	size_t i = 0, j = 0;
	while(i < vAlphaIDs.size() && j < vBetaIDs.size() ) {
		if(vAlphaIDs[i] == vBetaIDs[j]) {
			return PRT_INTRA;
		} else if(vAlphaIDs[i] < vBetaIDs[j]) {
			++i;
		} else {
			++j;
		}
	}
	return PRT_INTER;
}

int XLinkPeptideResult::getModificationId(std::string &strMod)
{
	SearchParameter *pParameter = ParameterReader::getInstance()->getParameter();
	int nModId = -1;
	for(size_t t = 0; t < pParameter->m_vVariableMods.size(); ++t) {
		if(!pParameter->m_vVariableMods[t].compare(strMod)) {
			nModId = t;
			break;
		}
	}
	for(size_t t = 0; t < pParameter->m_vFixedMods.size(); ++t) {
		if(!pParameter->m_vFixedMods[t].compare(strMod)) {
			nModId = t + pParameter->m_vVariableMods.size();
			break;
		}
	}
	return nModId;
}

int XLinkPeptideResult::unifiedPosition(int nIdx, bool bAlpha) {
	SearchParameter *pParameter = ParameterReader::getInstance()->getParameter();
	ModificationDict *pModDict = ModificationDict::getInstance();
	int nVarSize = pParameter->m_vVariableMods.size();
	if (bAlpha) { // alpha
		const string &strSq = m_stAlphaPep.m_stPep.m_strSq;
		const ModificationEntry &stEntry = m_stAlphaPep.m_stPep.m_vMods[nIdx];
		ModificationType eType = MT_PEP_NORMAL;
		if (stEntry.m_nId >= nVarSize) {
			eType = pModDict->getModificationType(
					pParameter->m_vFixedMods[stEntry.m_nId - nVarSize]);
		} else {
			eType = pModDict->getModificationType(
					pParameter->m_vVariableMods[stEntry.m_nId]);
		}

		if (stEntry.m_nSite > 0 && stEntry.m_nSite + 1 < (int) strSq.length()) {
			return stEntry.m_nSite + 1;
		} else if (stEntry.m_nSite == 0) {
			return (eType == MT_PEP_NORMAL) ? 1 : 0;
		} else if (stEntry.m_nSite + 1 == (int) strSq.length()) {
			return (eType == MT_PEP_NORMAL) ?
					stEntry.m_nSite + 1 : stEntry.m_nSite + 2;
		} else {
			ErrorInfo err("FinalResultIO", "unifiedPosition",
					"error modification position");
			throw runtime_error(err.get());
			return -1;
		}
	} else { // beta
		int nBase = m_stAlphaPep.m_stPep.m_strSq.length() + 4;
		const string &strSq = m_stBetaPep.m_stPep.m_strSq;
		const ModificationEntry &stEntry = m_stBetaPep.m_stPep.m_vMods[nIdx];
		ModificationType eType = MT_PEP_NORMAL;
		if (stEntry.m_nId >= nVarSize) {
			eType = pModDict->getModificationType(
					pParameter->m_vFixedMods[stEntry.m_nId - nVarSize]);
		} else {
			eType = pModDict->getModificationType(
					pParameter->m_vVariableMods[stEntry.m_nId]);
		}

		if (stEntry.m_nSite > 0 && stEntry.m_nSite + 1 < (int) strSq.length()) {
			return nBase + stEntry.m_nSite;
		} else if (stEntry.m_nSite == 0) {
			return (eType == MT_PEP_NORMAL) ? nBase : nBase - 1;
		} else if (stEntry.m_nSite + 1 == (int) strSq.length()) {
			return (eType == MT_PEP_NORMAL) ?
					nBase + stEntry.m_nSite : nBase + strSq.length();
		} else {
			ErrorInfo err("XLinkPeptideResult", "unifiedPosition",
					"error modification position");
			throw runtime_error(err.get());
			return -1;
		}
	}
}


void XLinkPeptideResult::getModifications(vector<string> &vMods) {
	const SearchParameter *pParameter = ParameterReader::getInstance()->getParameter();
	vMods.clear();

	const vector<string> &vFixMods = pParameter->m_vFixedMods;
	const vector<string> &vVarMods = pParameter->m_vVariableMods;
	int nVarSize = vVarMods.size();

	switch (m_ePepType) {
	case PT_COMMON:
	case PT_MONO:
	case PT_LOOP: {
		// modifications
		vector<ModificationEntry> &stMods =
				m_stAlphaPep.m_stPep.m_vMods;
		for (size_t tModIdx = 0; tModIdx < stMods.size(); ++tModIdx) {
			ostringstream oss;
			if (stMods[tModIdx].m_nId >= nVarSize) {
				oss << vFixMods[stMods[tModIdx].m_nId - nVarSize];
			} else {
				oss << vVarMods[stMods[tModIdx].m_nId];
			}
			oss << "("
					<< unifiedPosition(tModIdx, true) << ")";
			vMods.push_back(oss.str());
		}
		break;
	}
	case PT_XLINK: {
		const vector<ModificationEntry> &stMods1 =
				m_stAlphaPep.m_stPep.m_vMods;
		const vector<ModificationEntry> &stMods2 =
				m_stBetaPep.m_stPep.m_vMods;
		for (size_t tModIdx = 0; tModIdx < stMods1.size(); ++tModIdx) {
			ostringstream oss;
			if (stMods1[tModIdx].m_nId >= nVarSize) {
				oss << vFixMods[stMods1[tModIdx].m_nId - nVarSize];
			} else {
				oss << vVarMods[stMods1[tModIdx].m_nId];
			}
			oss << "("
					<< unifiedPosition(tModIdx, true) << ")";
			vMods.push_back(oss.str());
		}
		for (size_t tModIdx = 0; tModIdx < stMods2.size(); ++tModIdx) {
			ostringstream oss;
			if (stMods2[tModIdx].m_nId >= nVarSize) {
				oss << vFixMods[stMods2[tModIdx].m_nId - nVarSize];
			} else {
				oss << vVarMods[stMods2[tModIdx].m_nId];
			}
			oss << "("
					<< unifiedPosition(tModIdx, false) << ")";
			vMods.push_back(oss.str());
		}
		break;
	}
	}
}

void XLinkPeptideResult::getModifications(string &strMods) {
	strMods = "";
	vector<string> vMods;

	getModifications(vMods);
	for (size_t tIdx = 0; tIdx < vMods.size(); ++tIdx) {
		strMods += vMods[tIdx];
		if (tIdx + 1 < vMods.size())
			strMods += ";";
	}
	if (vMods.size() == 0)
		strMods = "null";
}

/** return a vector with all modifications in this peptide pair
 *  if both peptide do not contain any modifications, the value is 0
 *  for var modifications: 1, 2, ... var_num
 *  for fixed modifications: var_num+1, var_num+2, ... var_num+fix_num
 */
void XLinkPeptideResult::getModifications(set<size_t> &vModTypes) {
	if (m_ePepType == PT_XLINK) {
		if (m_stAlphaPep.m_stPep.m_vMods.empty()
				&& m_stBetaPep.m_stPep.m_vMods.empty()) {
			vModTypes.insert(0);
		} else {
			for (size_t tIdx = 0; tIdx < m_stAlphaPep.m_stPep.m_vMods.size();
					++tIdx) {
				vModTypes.insert(m_stAlphaPep.m_stPep.m_vMods[tIdx].m_nId + 1);
			}
			for (size_t tIdx = 0; tIdx < m_stBetaPep.m_stPep.m_vMods.size();
					++tIdx) {
				vModTypes.insert(m_stBetaPep.m_stPep.m_vMods[tIdx].m_nId + 1);
			}
		}
	} else {
		if (m_stAlphaPep.m_stPep.m_vMods.empty()) {
			vModTypes.insert(0);
		} else {
			for (size_t tIdx = 0; tIdx < m_stAlphaPep.m_stPep.m_vMods.size();
					++tIdx) {
				vModTypes.insert(m_stAlphaPep.m_stPep.m_vMods[tIdx].m_nId + 1);
			}
		}
	}
}

void XLinkPeptideResult::setPeptideSq(std::string &strSq)
{
	vector<string> vec;
	switch(m_ePepType) {
	case PT_XLINK:
	{
		split(strSq, '-', vec);
		size_t tLeft = vec[0].find('(');
		size_t tRight = vec[0].find(')');
		if(tLeft == string::npos || tRight == string::npos) {
			ErrorInfo err("XLinkPeptideResult", "parsePeptideSq", strSq.c_str());
			throw err.get();
			return;
		}
		m_stAlphaPep.m_stPep.m_strSq = vec[0].substr(0, tLeft);
		m_stAlphaPep.m_stLinkMod.m_nSite = atoi(vec[0].substr(tLeft+1, tRight-tLeft-1).c_str()) - 1;
		tLeft = vec[1].find('(');
		tRight = vec[1].find(')');
		if(tLeft == string::npos || tRight == string::npos) {
			ErrorInfo err("XLinkPeptideResult", "parsePeptideSq", strSq.c_str());
			throw err.get();
			return;
		}
		m_stBetaPep.m_stPep.m_strSq = vec[1].substr(0, tLeft);
		m_stBetaPep.m_stLinkMod.m_nSite = atoi(vec[1].substr(tLeft+1, tRight-tLeft-1).c_str()) - 1;
		break;
	}
	case PT_LOOP:
	{
		split(strSq, '-', vec);
		size_t tLeft = vec[0].find('(');
		size_t tRight = vec[0].find(')');
		if(tLeft == string::npos || tRight == string::npos) {
			ErrorInfo err("XLinkPeptideResult", "parsePeptideSq", strSq.c_str());
			throw err.get();
			return;
		}
		m_stAlphaPep.m_stPep.m_strSq = vec[0].substr(0, tLeft);
		m_stAlphaPep.m_stLinkMod.m_nSite = atoi(vec[0].substr(tLeft+1, tRight-tLeft-1).c_str()) - 1;
		tLeft = vec[0].find_first_of('(', tRight);
		tRight = vec[0].find_first_of(')', tLeft);
		m_stBetaPep.m_stLinkMod.m_nSite = atoi(vec[0].substr(tLeft+1, tRight-tLeft-1).c_str()) - 1;
		break;
	}
	case PT_MONO:
	{
		split(strSq, '-', vec);
		size_t tLeft = vec[0].find('(');
		size_t tRight = vec[0].find(')');
		if(tLeft == string::npos || tRight == string::npos) {
			ErrorInfo err("XLinkPeptideResult", "parsePeptideSq", strSq.c_str());
			throw err.get();
			return;
		}
		m_stAlphaPep.m_stPep.m_strSq = vec[0].substr(0, tLeft);
		m_stAlphaPep.m_stLinkMod.m_nSite = atoi(vec[0].substr(tLeft+1, tRight-tLeft-1).c_str()) - 1;
		break;
	}
	case PT_COMMON:
		m_stAlphaPep.m_stPep.m_strSq = vec[0];
		break;
	}
}

void XLinkPeptideResult::setPeptideLinkedId(string &strLinker)
{
	SearchParameter *pParameter = ParameterReader::getInstance()->getParameter();
	int tLinkerId = 0;
	for(size_t t = 0; t < pParameter->m_vLinkers.size(); ++t) {
		if(!pParameter->m_vLinkers[t].compare(strLinker)) {
			tLinkerId = t;
			break;
		}
	}
	if((int)m_ePepType >= 1) {
		m_stAlphaPep.m_stLinkMod.m_nId = tLinkerId;
	}

	if((int)m_ePepType >= 2) {
		m_stBetaPep.m_stLinkMod.m_nId = tLinkerId;
	}
}

void XLinkPeptideResult::setPeptideType(std::string &strType)
{
	if(!strType.compare("Cross-Linked")) {
		m_ePepType = PT_XLINK;
	} else if(!strType.compare("Loop-Linked")) {
		m_ePepType = PT_LOOP;
	} else if(!strType.compare("Mono-Linked")) {
		m_ePepType = PT_MONO;
	} else {
		m_ePepType = PT_COMMON;
	}
}

void XLinkPeptideResult::setModifications(string &strMods)
{
	if(!strMods.compare("null")) {
		return;
	}
	vector<string> vec;
	split(strMods, ';', vec);
	for(size_t tModId = 0; tModId < vec.size(); ++tModId) {
		ModificationEntry stEntry;
		size_t tLeft = vec[tModId].find('(');
		size_t tRight = vec[tModId].find(')');
		if(tLeft == string::npos || tRight == string::npos) {
			ErrorInfo err("XLinkPeptideResult", "setModifications", strMods);
			throw err.get();
			return;
		}
		string strMod = vec[tModId].substr(0, tLeft);
		stEntry.m_nId = getModificationId(strMod);
		int nSite = atoi(vec[tModId].substr(tLeft+1, tRight-tLeft-1).c_str());
		if(nSite <= (int)m_stAlphaPep.m_stPep.m_strSq.length()+1) {
			if(nSite == 0) {
				stEntry.m_nSite = nSite;
			} else if(nSite == (int)m_stAlphaPep.m_stPep.m_strSq.length()+1){
				stEntry.m_nSite = nSite - 2;
			} else {
				stEntry.m_nSite = nSite - 1;
			}
			m_stAlphaPep.m_stPep.m_vMods.push_back(stEntry);
		} else {
			if(nSite == (int)m_stAlphaPep.m_stPep.m_strSq.length()+3) {
				stEntry.m_nSite = 0;
			} else if(nSite == (int)(m_stAlphaPep.m_stPep.m_strSq.length()+3+m_stBetaPep.m_stPep.m_strSq.length()+1)) {
				stEntry.m_nSite = nSite - m_stAlphaPep.m_stPep.m_strSq.length() - 5;
			} else {
				stEntry.m_nSite = nSite - m_stAlphaPep.m_stPep.m_strSq.length() - 4;
			}
			m_stBetaPep.m_stPep.m_vMods.push_back(stEntry);
		}
	}
}

void XLinkPeptideResult::clear() {
	m_stAlphaPep.clear();
	m_stBetaPep.clear();
	m_ePepType = PT_COMMON;
	m_eProType = PRT_NONE;
	m_eTDType = TDT_F;
	m_lfScore = 0.0;
	m_lfEvalue = 0.0;
	m_lfQvalue = 0.0;
	m_lfCalcMH = 0.0;
	m_stFeatureInfo.clear();
}

string XLinkPeptideResult::toString() {
	ostringstream oss;
	oss << m_ePepType << " " << m_lfCalcMH << " " << m_lfScore << " "
			<< m_stAlphaPep << " " << m_stBetaPep << " " << m_stFeatureInfo;
	return oss.str();
}

string XLinkPeptideResult::toDebugString() {
	ostringstream oss;
	oss << "peptide_type=" << m_ePepType << "\n" << "calculate_MH="
			<< m_lfCalcMH << "\n" << "score=" << m_lfScore << "\n"
			<< "alpha_pep: \n" << m_stAlphaPep.toDebugString() << "\n"
			<< "beta_pep: \n" << m_stBetaPep.toDebugString() << "\n"
			<< "features: \n" << m_stFeatureInfo.toDebugString() << "\n";
	return oss.str();
}

double XLinkPeptideResult::getGodelCode() const {
	if (m_ePepType == PT_XLINK) {
		ConstantStore *pConst = ConstantStore::getInstance();
		double lfCode = m_stAlphaPep.m_stPep.getGodelCode();
		size_t tAlphaLen = m_stAlphaPep.m_stPep.m_strSq.length();
		lfCode += ('Z' + 1 - 'A') * pConst->m_lfLogPrime[tAlphaLen]; // connector
		for (size_t tIdx = 0; tIdx < m_stBetaPep.m_stPep.m_strSq.length();
				++tIdx) {
			lfCode += (m_stBetaPep.m_stPep.m_strSq[tIdx] - 'A')
					* pConst->m_lfLogPrime[tIdx + tAlphaLen + 1];
		}
		return lfCode;
	} else {
		return m_stAlphaPep.m_stPep.getGodelCode();
	}
}

double XLinkPeptideResult::calculateXLinkPeptideMass()
{
	m_stAlphaPep.m_stPep.calculatePeptideMass();
	if(m_ePepType == PT_XLINK)
		m_stBetaPep.m_stPep.calculatePeptideMass();
	SearchParameter *pParameter = ParameterReader::getInstance()->getParameter();
	XLinkerDict *pDict = XLinkerDict::getInstance();
	string strLinker = pParameter->m_vLinkers[m_stAlphaPep.m_stLinkMod.m_nId];

	double lfMass = m_stAlphaPep.m_stPep.m_lfMass + H2O_MONO_MASS;
	if(m_ePepType == PT_XLINK) {
		lfMass += m_stBetaPep.m_stPep.m_lfMass + H2O_MONO_MASS;
		lfMass += pDict->getMonoDiff(strLinker);
		m_stAlphaPep.m_lfOpenMass = m_stBetaPep.m_stPep.m_lfMass + H2O_MONO_MASS + pDict->getMonoDiff(strLinker);
		m_stBetaPep.m_lfOpenMass = m_stAlphaPep.m_stPep.m_lfMass + H2O_MONO_MASS + pDict->getMonoDiff(strLinker);
	} else if(m_ePepType == PT_LOOP) {
		lfMass += pDict->getMonoDiff(strLinker);
		m_stAlphaPep.m_lfOpenMass = pDict->getMonoDiff(strLinker);
	} else if(m_ePepType == PT_MONO) {
		lfMass += pDict->getMLMonoDiff(strLinker);
		m_stAlphaPep.m_lfOpenMass = pDict->getMLMonoDiff(strLinker);
	}
	m_lfCalcMH = lfMass + PROTON_MASS;
	return lfMass;
}

OpenMatchResult::OpenMatchResult(size_t tSize) {
	m_vPepItemsQ.reserve(tSize);
	m_vPepsQ.reserve(tSize);
}

OpenMatchResult::~OpenMatchResult() {

}

void OpenMatchResult::buildMinHeap(size_t tCount,
		const PeptideResult &stResult) {
	PeptideItem stItem(tCount, m_vPepItemsQ.size());
	m_vPepItemsQ.push_back(stItem);
	m_vPepsQ.push_back(stResult);

	for (int i = ((int) m_vPepItemsQ.size() - 2) >> 1; i >= 0; --i)
		heapify(i);
}

void OpenMatchResult::heapify(size_t tCount, const PeptideResult &stResult) {
	PeptideResult &topResult = m_vPepsQ[m_vPepItemsQ[0].m_tPepId];
	if (topResult.m_lfScore < stResult.m_lfScore) {
		m_vPepItemsQ[0].m_tCount = tCount;
		m_vPepsQ[m_vPepItemsQ[0].m_tPepId] = stResult;
		heapify(0);
	}
}

void OpenMatchResult::heapify(int nIdx) {
	int nSize = (int) m_vPepItemsQ.size();
	int nPos = (nIdx << 1) + 1;

	while (nPos < nSize) {
		if (nPos + 1 < nSize && m_vPepItemsQ[nPos + 1] < m_vPepItemsQ[nPos])
			++nPos;
		if (m_vPepItemsQ[nIdx] < m_vPepItemsQ[nPos])
			break;
		swap(m_vPepItemsQ[nIdx], m_vPepItemsQ[nPos]);
		nIdx = nPos;
		nPos = (nIdx << 1) + 1;
	}
}

void OpenMatchResult::removeRedundant() {
	try {
		map<string, vector<PeptideItem> > mpNonRedundant;
		map<string, vector<PeptideResult> > mpResults;
		for (size_t i = 0; i < m_vPepItemsQ.size(); ++i) {
			string strPepSq = m_vPepsQ[m_vPepItemsQ[i].m_tPepId].m_stPep.m_strSq;
			if (mpNonRedundant.find(strPepSq) == mpNonRedundant.end()) {
				mpNonRedundant[strPepSq].push_back(m_vPepItemsQ[i]);
				mpResults[strPepSq].push_back(
						m_vPepsQ[m_vPepItemsQ[i].m_tPepId]);
			} else {
				size_t j = 0;
				for (j = 0; j < mpResults[strPepSq].size(); ++j) {
					if (mpResults[strPepSq][j].m_stLinkMod.m_nId
							== m_vPepsQ[m_vPepItemsQ[i].m_tPepId].m_stLinkMod.m_nId)
						break;
				}
				if (j >= mpResults[strPepSq].size()) {
					mpNonRedundant[strPepSq].push_back(m_vPepItemsQ[i]);
					mpResults[strPepSq].push_back(
							m_vPepsQ[m_vPepItemsQ[i].m_tPepId]);
				}
			}
		}

		m_vPepItemsQ.clear();
		m_vPepsQ.clear();

		size_t tIdx = 0;
		map<string, vector<PeptideItem> >::iterator it = mpNonRedundant.begin();
		for (; it != mpNonRedundant.end(); ++it) {
			for (size_t j = 0; j < mpNonRedundant[(*it).first].size(); ++j) {
				m_vPepsQ.push_back(mpResults[(*it).first][j]);
				(it->second)[j].m_tPepId = tIdx;
				m_vPepItemsQ.push_back((*it).second[j]);
				++tIdx;
			}
		}

		for (int i = ((int) m_vPepItemsQ.size() - 2) >> 1; i >= 0; --i)
			heapify(i);
	} catch (exception &e) {
		ErrorInfo err("OpenMatchResult", "removeRedundant",
				"caught an exception.", e);
		throw runtime_error(err.get());
	} catch (...) {
		ErrorInfo err("OpenMatchResult", "removeRedundant",
				"caught an unknown exception.");
		throw runtime_error(err.get());
	}
}

void OpenMatchResult::swapOut() {
	vector<PeptideItem>().swap(m_vPepItemsQ);
	vector<PeptideResult>().swap(m_vPepsQ);
}

/* begin: struct GroupCounter */
GroupCounter::GroupCounter() :
		m_lfSystemaicError(0.0) {
	clear();
}

void GroupCounter::clear() {
	fill(m_vPepCounts, m_vPepCounts + PEPTIDE_TYPE_NUM, 0);
	fill(m_vProCounts, m_vProCounts + PROTEIN_TYPE_NUM, 0);
	fill(m_vChargeCounts, m_vChargeCounts + CHARGE_RANGE, 0);
	fill(m_vErrorCounts, m_vErrorCounts + PRECURSOR_ERROR_INTERVAL_NUM, 0);
	fill(m_vModCounts, m_vModCounts + MAX_MODIFY_NUM, 0);
	m_lfSystemaicError = 0.0;
}
/* end: struct GroupCounter */

/* begin: struct XLinkPSMFeature */
XLinkPSMFeature::XLinkPSMFeature() :
		m_lfScore(0.0), m_lfEvalue(0.0) {
}

bool operator<(const XLinkPSMFeature &stItem1, const XLinkPSMFeature &stItem2) {
	return stItem1.m_lfScore < stItem2.m_lfScore;
}

bool operator<=(const XLinkPSMFeature &stItem1,
		const XLinkPSMFeature &stItem2) {
	return stItem1.m_lfScore <= stItem2.m_lfScore;
}

bool operator>(const XLinkPSMFeature &stItem1, const XLinkPSMFeature &stItem2) {
	return stItem1.m_lfScore > stItem2.m_lfScore;
}

bool operator>=(const XLinkPSMFeature &stItem1,
		const XLinkPSMFeature &stItem2) {
	return stItem1.m_lfScore >= stItem2.m_lfScore;
}
/* end: struct XLinkPSMFeature */

/* begin: struct PSMSortItem */

bool operator<(const PSMSortItem &stItem1, const PSMSortItem &stItem2) {
	if (stItem1.m_lfQvalue == stItem2.m_lfQvalue) {
		return stItem1.m_lfScore < stItem2.m_lfScore;
	} else {
		return stItem1.m_lfQvalue > stItem2.m_lfQvalue;
	}
}

bool operator<=(const PSMSortItem &stItem1, const PSMSortItem &stItem2) {
	if (stItem1.m_lfQvalue == stItem2.m_lfQvalue) {
		return stItem1.m_lfScore <= stItem2.m_lfScore;
	} else {
		return stItem1.m_lfQvalue > stItem2.m_lfQvalue;
	}
}

bool operator>(const PSMSortItem &stItem1, const PSMSortItem &stItem2) {
	if (stItem1.m_lfQvalue == stItem2.m_lfQvalue) {
		return stItem1.m_lfScore > stItem2.m_lfScore;
	} else {
		return stItem1.m_lfQvalue < stItem2.m_lfQvalue;
	}
}

bool operator>=(const PSMSortItem &stItem1, const PSMSortItem &stItem2) {
	if (stItem1.m_lfQvalue == stItem2.m_lfQvalue) {
		return stItem1.m_lfScore >= stItem2.m_lfScore;
	} else {
		return stItem1.m_lfQvalue < stItem2.m_lfQvalue;
	}
}

/* end: struct PSMSortItem */

ReportSummary::ReportSummary() :
		m_tSpectraTotal(0), m_tSpecTotal(0) {
	std::fill(m_vSpectraNum, m_vSpectraNum + PEPTIDE_TYPE_NUM, 0);
	std::fill(m_vSpecNum, m_vSpecNum + PEPTIDE_TYPE_NUM, 0);
	std::fill(m_vPepNum, m_vPepNum + PEPTIDE_TYPE_NUM, 0);
	std::fill(m_vSitesNum, m_vSitesNum + PEPTIDE_TYPE_NUM, 0);
	std::fill(m_vGroupNum, m_vGroupNum + PEPTIDE_TYPE_NUM, 0);
}

} /* namespace sdk */
