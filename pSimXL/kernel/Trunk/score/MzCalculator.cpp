
#include "../include/scorers.h"

using namespace std;

namespace sdk
{

MzTriple::MzTriple() :
		nMz(0), lfMz(0.0), yIonTypeOrder(UCHAR_MAX),
		yPepPosOrder1(UCHAR_MAX), yPepPosOrder2(UCHAR_MAX),
		yNTerm1(FALSE), yNTerm2(FALSE),
		yAddToTag(FALSE), yContainLinker(FALSE)
{
}

ostream &operator<<(ostream &stOut, const MzTriple &stTriple)
{
	stOut<<stTriple.nMz<<" "
		 <<(int)stTriple.yIonTypeOrder<<" "
		 <<(int)stTriple.yPepPosOrder1<<" "
		 <<(int)stTriple.yPepPosOrder2<<" "
		 <<(int)stTriple.yNTerm1<<" "
		 <<(int)stTriple.yNTerm2<<" "
		 <<(int)stTriple.yAddToTag<<" "
		 <<(int)stTriple.yContainLinker;
	return stOut;
}

void MzTriple::clear()
{
	nMz = 0;
	yIonTypeOrder = yPepPosOrder1 = yPepPosOrder2 = UCHAR_MAX;
	yNTerm1 = yNTerm2 = yAddToTag = yContainLinker = FALSE;
}

string MzTriple::toString()
{
	ostringstream oss;
	oss<<nMz<<" "
	   <<(int)yIonTypeOrder<<" "
	   <<(int)yPepPosOrder1<<" "
	   <<(int)yPepPosOrder2<<" "
	   <<(int)yNTerm1<<" "
	   <<(int)yNTerm2<<" "
	   <<(int)yAddToTag<<" "
	   <<(int)yContainLinker;
	return oss.str();
}

PeptideAAMass::PeptideAAMass() :
		m_lfPepMass(0.0), m_nPepLen(0)
{
	m_nPepAAMass.reserve(MAX_PEPTIDE_LEN + 5);
	m_bPepAAMass.reserve(MAX_PEPTIDE_LEN + 5);
	m_yPepAAMass.reserve(MAX_PEPTIDE_LEN + 5);
}

PeptideAAMass::~PeptideAAMass()
{
	clear();
}

void PeptideAAMass::resize(size_t tSize)
{
	m_nPepAAMass.resize(tSize, 0);
	m_bPepAAMass.resize(tSize, 0);
	m_yPepAAMass.resize(tSize, 0);
}

void PeptideAAMass::clear()
{
	m_lfPepMass = 0;
	m_nPepLen = 0;
	m_nPepAAMass.clear();
	m_bPepAAMass.clear();
	m_yPepAAMass.clear();
}

MzCalculator::MzCalculator() :
		m_pAAMass(NULL), m_pParameter(NULL), m_pPeptide(NULL), m_pPeptideResult(NULL),
		m_pXLinkPeptideResult(NULL), m_pAlphaPep(NULL), m_pBetaPep(NULL)
{
	try {
		m_pAAMass = ConstantStore::getInstance()->m_lfAAMass;
		m_pParameter = ParameterReader::getInstance()->getParameter();

		ModificationDict *pDict = ModificationDict::getInstance();
		for(size_t i = 0; i < m_pParameter->m_vVariableMods.size(); ++i) {
			const Modification &mod = pDict->getModification(m_pParameter->m_vVariableMods[i]);
			m_vMods.push_back(mod);
		}
		for(size_t i = 0; i < m_pParameter->m_vFixedMods.size(); ++i) {
			const Modification &mod = pDict->getModification(m_pParameter->m_vFixedMods[i]);
			m_vMods.push_back(mod);
		}
	} catch(exception &e) {
		ErrorInfo err("MzCalculator", "MzCalculator", "initialize MzCaculator failed!");
		throw runtime_error(err.get());
	} catch(...) {
		ErrorInfo err("MzCalculator", "MzCalculator", "caught an unknown exception.");
		throw runtime_error(err.get());
	}
}

MzCalculator::~MzCalculator()
{
	if(m_pAAMass) {
		m_pAAMass = NULL;
	}

	if(m_pParameter) {
		m_pParameter = NULL;
	}

	close();
}

void MzCalculator::init(Peptide *pPeptide)
{
	try {
		close();

		m_pPeptide = pPeptide;
		m_pAlphaPep = new PeptideAAMass;
		m_pAlphaPep->m_nPepLen = m_pPeptide->m_strSq.length();
		m_pAlphaPep->resize(m_pAlphaPep->m_nPepLen);

		// basic mass
		for(int i = 0; i < m_pAlphaPep->m_nPepLen; ++i) {
			m_pAlphaPep->m_nPepAAMass[i] = m_pAAMass[m_pPeptide->m_strSq[i]-'A'];
		}

		// modification
		vector<ModificationEntry> &vMods = m_pPeptideResult->m_stPep.m_vMods;
		for(size_t i = 0; i < vMods.size(); ++i) {
			m_pAlphaPep->m_nPepAAMass[vMods[i].m_nSite] += m_vMods[vMods[i].m_nId].lfMonoDiff;
		}

		// compute b and y ions masses
		m_pAlphaPep->m_bPepAAMass[0] = m_pAlphaPep->m_nPepAAMass[0] + PROTON_MASS;
		m_pAlphaPep->m_yPepAAMass[0] = m_pAlphaPep->m_nPepAAMass[m_pAlphaPep->m_nPepLen-1] + PROTON_MASS + H2O_MONO_MASS;
		for(int i = 1; i < m_pAlphaPep->m_nPepLen-1; ++i) {
			m_pAlphaPep->m_bPepAAMass[i] = m_pAlphaPep->m_bPepAAMass[i-1] + m_pAlphaPep->m_nPepAAMass[i];
			m_pAlphaPep->m_yPepAAMass[i] = m_pAlphaPep->m_yPepAAMass[i-1] + m_pAlphaPep->m_nPepAAMass[m_pAlphaPep->m_nPepLen-1-i];
		}
	} catch(exception &e) {
		ErrorInfo err("MzCalculator", "init", "initialize Peptide failed!");
		throw runtime_error(err.get());
	} catch(...) {
		ErrorInfo err("MzCalculator", "init", "caught an unknown exception.");
		throw runtime_error(err.get());
	}
}

void MzCalculator::init(PeptideResult *pPeptide)
{
	try {
		close();

		m_pPeptideResult = pPeptide;
		m_pAlphaPep = new PeptideAAMass;
		m_pAlphaPep->m_nPepLen = m_pPeptideResult->m_stPep.m_strSq.length();
		m_pAlphaPep->resize(m_pAlphaPep->m_nPepLen);

		// basic mass
		for(int i = 0; i < m_pAlphaPep->m_nPepLen; ++i) {
			m_pAlphaPep->m_nPepAAMass[i] = m_pAAMass[m_pPeptideResult->m_stPep.m_strSq[i]-'A'];
		}

		// modification
		vector<ModificationEntry> &vMods = m_pPeptideResult->m_stPep.m_vMods;
		for(size_t i = 0; i < vMods.size(); ++i) {
			m_pAlphaPep->m_nPepAAMass[vMods[i].m_nSite] += m_vMods[vMods[i].m_nId].lfMonoDiff;
		}

		ModificationEntry &stLinker = m_pPeptideResult->m_stLinkMod;
		// open mass
		if(stLinker.m_nSite != -1)
			m_pAlphaPep->m_nPepAAMass[stLinker.m_nSite] += m_pPeptideResult->m_lfOpenMass;

		// compute b and y ions masses
		m_pAlphaPep->m_bPepAAMass[0] = m_pAlphaPep->m_nPepAAMass[0] + PROTON_MASS;
		m_pAlphaPep->m_yPepAAMass[0] = m_pAlphaPep->m_nPepAAMass[m_pAlphaPep->m_nPepLen-1] + PROTON_MASS + H2O_MONO_MASS;
		for(int i = 1; i < m_pAlphaPep->m_nPepLen-1; ++i) {
			m_pAlphaPep->m_bPepAAMass[i] = m_pAlphaPep->m_bPepAAMass[i-1] + m_pAlphaPep->m_nPepAAMass[i];
			m_pAlphaPep->m_yPepAAMass[i] = m_pAlphaPep->m_yPepAAMass[i-1] + m_pAlphaPep->m_nPepAAMass[m_pAlphaPep->m_nPepLen-1-i];
		}
	} catch(exception &e) {
		ErrorInfo err("MzCalculator", "init", "initialize PeptideResult failed!");
		throw runtime_error(err.get());
	} catch(...) {
		ErrorInfo err("MzCalculator", "init", "caught an unknown exception.");
		throw runtime_error(err.get());
	}
}

void MzCalculator::init(XLinkPeptideResult *pPeptide)
{
	try {
		close();

		m_pXLinkPeptideResult = pPeptide;
		m_pAlphaPep = new PeptideAAMass;

		// peptide length
		m_pAlphaPep->m_nPepLen = m_pXLinkPeptideResult->m_stAlphaPep.m_stPep.m_strSq.length();
		m_pAlphaPep->resize(m_pAlphaPep->m_nPepLen);

		// peptide mass
		m_pAlphaPep->m_lfPepMass = m_pXLinkPeptideResult->m_stAlphaPep.m_stPep.m_lfMass + H2O_MONO_MASS;

		// basic mass
		for(int i = 0; i < m_pAlphaPep->m_nPepLen; ++i) {
			m_pAlphaPep->m_nPepAAMass[i] = m_pAAMass[m_pXLinkPeptideResult->m_stAlphaPep.m_stPep.m_strSq[i]-'A'];
		}

		// modifications
		vector<ModificationEntry> &vMods = m_pXLinkPeptideResult->m_stAlphaPep.m_stPep.m_vMods;
		for(size_t i = 0; i < vMods.size(); ++i) {
			m_pAlphaPep->m_nPepAAMass[vMods[i].m_nSite] += m_vMods[vMods[i].m_nId].lfMonoDiff;
		}

		if(m_pXLinkPeptideResult->m_ePepType == PT_XLINK) {
			m_pBetaPep = new PeptideAAMass;

			// peptide length
			m_pBetaPep->m_nPepLen = m_pXLinkPeptideResult->m_stBetaPep.m_stPep.m_strSq.length();
			m_pBetaPep->resize(m_pBetaPep->m_nPepLen);

			// peptide mass
			m_pBetaPep->m_lfPepMass = m_pXLinkPeptideResult->m_stBetaPep.m_stPep.m_lfMass + H2O_MONO_MASS;

			// basic mass
			for(int i = 0; i < m_pBetaPep->m_nPepLen; ++i) {
				m_pBetaPep->m_nPepAAMass[i] = m_pAAMass[m_pXLinkPeptideResult->m_stBetaPep.m_stPep.m_strSq[i]-'A'];
			}

			// modifications
			vector<ModificationEntry> &vMods = m_pXLinkPeptideResult->m_stBetaPep.m_stPep.m_vMods;
			for(size_t i = 0; i < vMods.size(); ++i) {
				m_pBetaPep->m_nPepAAMass[vMods[i].m_nSite] += m_vMods[vMods[i].m_nId].lfMonoDiff;
			}

			// b ions of alpha peptide
			double lfTmpMass(PROTON_MASS);
			for(int nSite = 0; nSite < m_pAlphaPep->m_nPepLen; ++nSite) {
				lfTmpMass += m_pAlphaPep->m_nPepAAMass[nSite];
				if(nSite == m_pXLinkPeptideResult->m_stAlphaPep.m_stLinkMod.m_nSite) {
					lfTmpMass += m_pXLinkPeptideResult->m_stAlphaPep.m_lfOpenMass;
				}
				m_pAlphaPep->m_bPepAAMass[nSite] = lfTmpMass;
			}

			// y ions of alpha peptide
			lfTmpMass = H2O_MONO_MASS + PROTON_MASS;
			for(int nSite = m_pAlphaPep->m_nPepLen-1; nSite >= 0; --nSite) {
				lfTmpMass += m_pAlphaPep->m_nPepAAMass[nSite];
				if(nSite == m_pXLinkPeptideResult->m_stAlphaPep.m_stLinkMod.m_nSite) {
					lfTmpMass += m_pXLinkPeptideResult->m_stAlphaPep.m_lfOpenMass;
				}
				m_pAlphaPep->m_yPepAAMass[m_pAlphaPep->m_nPepLen-nSite-1] = lfTmpMass;
			}

			// b ions of beta peptide
			lfTmpMass = PROTON_MASS;
			for(int nSite = 0; nSite < m_pBetaPep->m_nPepLen; ++nSite) {
				lfTmpMass += m_pBetaPep->m_nPepAAMass[nSite];
				if(nSite == m_pXLinkPeptideResult->m_stBetaPep.m_stLinkMod.m_nSite) {
					lfTmpMass += m_pXLinkPeptideResult->m_stBetaPep.m_lfOpenMass;
				}
				m_pBetaPep->m_bPepAAMass[nSite] = lfTmpMass;
			}

			// y ions of beta peptide
			lfTmpMass = H2O_MONO_MASS + PROTON_MASS;
			for(int nSite = m_pBetaPep->m_nPepLen-1; nSite >= 0; --nSite) {
				lfTmpMass += m_pBetaPep->m_nPepAAMass[nSite];
				if(nSite == m_pXLinkPeptideResult->m_stBetaPep.m_stLinkMod.m_nSite) {
					lfTmpMass += m_pXLinkPeptideResult->m_stBetaPep.m_lfOpenMass;
				}
				m_pBetaPep->m_yPepAAMass[m_pBetaPep->m_nPepLen-nSite-1] = lfTmpMass;
			}
		} else if(m_pXLinkPeptideResult->m_ePepType == PT_LOOP) {
			// b ions
			double lfTmpMass(PROTON_MASS);
			int nSite(0);
			while(nSite < m_pXLinkPeptideResult->m_stAlphaPep.m_stLinkMod.m_nSite) {
				lfTmpMass += m_pAlphaPep->m_nPepAAMass[nSite];
				m_pAlphaPep->m_bPepAAMass[nSite] = lfTmpMass;
				++nSite;
			}
			while(nSite < m_pXLinkPeptideResult->m_stBetaPep.m_stLinkMod.m_nSite) {
				lfTmpMass += m_pAlphaPep->m_nPepAAMass[nSite];
				m_pAlphaPep->m_bPepAAMass[nSite] = PROTON_MASS + m_pAlphaPep->m_lfPepMass + m_pXLinkPeptideResult->m_stAlphaPep.m_lfOpenMass;
				++nSite;
			}
			lfTmpMass += m_pXLinkPeptideResult->m_stAlphaPep.m_lfOpenMass;
			while(nSite < m_pAlphaPep->m_nPepLen) {
				lfTmpMass += m_pAlphaPep->m_nPepAAMass[nSite];
				m_pAlphaPep->m_bPepAAMass[nSite] = lfTmpMass;
				++nSite;
			}

			// y ions
			lfTmpMass = H2O_MONO_MASS + PROTON_MASS;
			nSite = m_pAlphaPep->m_nPepLen - 1;
			while(nSite > m_pXLinkPeptideResult->m_stBetaPep.m_stLinkMod.m_nSite) {
				lfTmpMass += m_pAlphaPep->m_nPepAAMass[nSite];
				m_pAlphaPep->m_yPepAAMass[m_pAlphaPep->m_nPepLen-nSite-1] = lfTmpMass;
				--nSite;
			}
			while(nSite > m_pXLinkPeptideResult->m_stAlphaPep.m_stLinkMod.m_nSite) {
				lfTmpMass += m_pAlphaPep->m_nPepAAMass[nSite];
				m_pAlphaPep->m_yPepAAMass[m_pAlphaPep->m_nPepLen-nSite-1] = PROTON_MASS + m_pAlphaPep->m_lfPepMass + m_pXLinkPeptideResult->m_stAlphaPep.m_lfOpenMass;
				--nSite;
			}
			lfTmpMass += m_pXLinkPeptideResult->m_stAlphaPep.m_lfOpenMass;
			while(nSite >= 0) {
				lfTmpMass += m_pAlphaPep->m_nPepAAMass[nSite];
				m_pAlphaPep->m_yPepAAMass[m_pAlphaPep->m_nPepLen-nSite-1] = lfTmpMass;
				--nSite;
			}
		} else { // common and mono cases
			ModificationEntry &stLinker = m_pXLinkPeptideResult->m_stAlphaPep.m_stLinkMod;
			// open mass
			if(stLinker.m_nSite != -1)
				m_pAlphaPep->m_nPepAAMass[stLinker.m_nSite] += m_pXLinkPeptideResult->m_stAlphaPep.m_lfOpenMass;

			// compute b and y ions masses
			m_pAlphaPep->m_bPepAAMass[0] = m_pAlphaPep->m_nPepAAMass[0] + PROTON_MASS;
			m_pAlphaPep->m_yPepAAMass[0] = m_pAlphaPep->m_nPepAAMass[m_pAlphaPep->m_nPepLen-1] + PROTON_MASS + H2O_MONO_MASS;
			for(int i = 1; i < m_pAlphaPep->m_nPepLen-1; ++i) {
				m_pAlphaPep->m_bPepAAMass[i] = m_pAlphaPep->m_bPepAAMass[i-1] + m_pAlphaPep->m_nPepAAMass[i];
				m_pAlphaPep->m_yPepAAMass[i] = m_pAlphaPep->m_yPepAAMass[i-1] + m_pAlphaPep->m_nPepAAMass[m_pAlphaPep->m_nPepLen-1-i];
			}
		}
	} catch(exception &e) {
		ErrorInfo err("MzCalculator", "init", "initialize XLinkPeptideResult failed!");
		throw runtime_error(err.get());
	} catch(...) {
		ErrorInfo err("MzCalculator", "init", "caught an unknown exception.");
		throw runtime_error(err.get());
	}
}

void MzCalculator::attachNTermIons(vector<MzTriple> &vIonMz, const IonType &stType, byte yOrder)
{
	if(m_pPeptide) {
	} else if(m_pPeptideResult) {
		// n terminal ions
		for(int nIdx = 0; nIdx < m_pAlphaPep->m_nPepLen-1; ++nIdx) {
			MzTriple triple;
			triple.lfMz = (m_pAlphaPep->m_bPepAAMass[nIdx] + (stType.nCharge - 1) * PROTON_MASS - stType.lfLoss) / stType.nCharge;
			triple.nMz =  triple.lfMz * MZ_MULTIPLIER + 0.5;
			triple.yIonTypeOrder = yOrder;
			triple.yPepPosOrder1 = nIdx;
			triple.yAddToTag = TRUE;
			vIonMz.push_back(triple);
		}

		// for neutral loss
		for(size_t j = 0; j < m_pPeptideResult->m_stPep.m_vMods.size(); ++j) {
			int nId = m_pPeptideResult->m_stPep.m_vMods[j].m_nId;
			int nSite = m_pPeptideResult->m_stPep.m_vMods[j].m_nSite;
			for(size_t k = 0; k < m_vMods[nId].vMonoNeutralLossDiff.size(); ++k) {
				for(int nIdx = nSite; nIdx < m_pAlphaPep->m_nPepLen-1; ++nIdx) {
					MzTriple triple;
					double lfMz = m_pAlphaPep->m_bPepAAMass[nIdx] - m_vMods[nId].vMonoNeutralLossDiff[k];
					if(lfMz < 0)
						continue;
					triple.lfMz = (lfMz + (stType.nCharge - 1) * PROTON_MASS - stType.lfLoss) / stType.nCharge;
					triple.nMz = triple.lfMz * MZ_MULTIPLIER + 0.5;
					triple.yIonTypeOrder = yOrder;
					triple.yPepPosOrder1 = nIdx;
					triple.yAddToTag = FALSE;
					vIonMz.push_back(triple);
				}
			}
		}
	} else if(m_pXLinkPeptideResult) {
		PeptideType eType = m_pXLinkPeptideResult->m_ePepType;
		int nAlphaSite(m_pXLinkPeptideResult->m_stAlphaPep.m_stLinkMod.m_nSite);
		int nBetaSite(m_pXLinkPeptideResult->m_stBetaPep.m_stLinkMod.m_nSite);

		// N terminal ions of alpha peptide
		for(int nIdx = 0; nIdx < m_pAlphaPep->m_nPepLen-1; ++nIdx) {
			bool bContainLinker(false);
			if (eType == PT_LOOP) {
				if (nIdx >= nAlphaSite && nIdx < nBetaSite)
					continue;
			} else if(eType == PT_XLINK) {
				if (nIdx >= nAlphaSite) { // contain linker
					if (stType.nHasLinker == ONLY_COMMON_ION) // only consider common ion
						continue;
					bContainLinker = true;
				} else { // not contain linker
					if (stType.nHasLinker == ONLY_LINKER_ION) // only consider x-link ion
						continue;
					bContainLinker = false;
				}
			}

			MzTriple triple;
			triple.lfMz = (m_pAlphaPep->m_bPepAAMass[nIdx] + (stType.nCharge - 1) * PROTON_MASS - stType.lfLoss) / stType.nCharge;
			triple.nMz = triple.lfMz * MZ_MULTIPLIER + 0.5;
			triple.yIonTypeOrder = yOrder;
			triple.yPepPosOrder1 = nIdx;
			triple.yNTerm1 = TRUE;
			triple.yAddToTag = stType.bContinuity;
			triple.yContainLinker = bContainLinker;
			vIonMz.push_back(triple);
		}

		for(size_t j = 0; j < m_pXLinkPeptideResult->m_stAlphaPep.m_stPep.m_vMods.size(); ++j) {
			// for neutral loss
			int nId = m_pXLinkPeptideResult->m_stAlphaPep.m_stPep.m_vMods[j].m_nId;
			int nSite = m_pXLinkPeptideResult->m_stAlphaPep.m_stPep.m_vMods[j].m_nSite;
			for(size_t k = 0; k < m_vMods[nId].vMonoNeutralLossDiff.size(); ++k) {
				for(int nIdx = nSite; nIdx < m_pAlphaPep->m_nPepLen-1; ++nIdx) {
					bool bContainLinker(false);
					if (eType == PT_LOOP) {
						if (nIdx >= nAlphaSite && nIdx < nBetaSite)
							continue;
					} else if(eType == PT_XLINK) {
						if (nIdx >= nAlphaSite) { // contain linker
							if (stType.nHasLinker == ONLY_COMMON_ION) // only consider common ion
								continue;
							bContainLinker = true;
						} else { // not contain linker
							if (stType.nHasLinker == ONLY_LINKER_ION) // only consider x-link ion
								continue;
							bContainLinker = false;
						}
					}

					MzTriple triple;
					double lfMz = m_pAlphaPep->m_bPepAAMass[nIdx] - m_vMods[nId].vMonoNeutralLossDiff[k];
					if(lfMz < 0)
						continue;
					triple.lfMz =  (lfMz + (stType.nCharge - 1) * PROTON_MASS - stType.lfLoss) / stType.nCharge;
					triple.nMz = triple.lfMz * MZ_MULTIPLIER + 0.5;
					triple.yIonTypeOrder = yOrder;
					triple.yPepPosOrder1 = nIdx;
					triple.yNTerm1 = TRUE;
					triple.yAddToTag = FALSE;
					triple.yContainLinker = bContainLinker;
					vIonMz.push_back(triple);
				}
			}
		}

		// n terminal ions of beta peptide, if PT_XLINK is set on
		if(eType == PT_XLINK) {
			for(int nIdx = 0; nIdx < m_pBetaPep->m_nPepLen-1; ++nIdx) {
				bool bContainLinker(false);
				if(nIdx >= nBetaSite) {
					if(stType.nHasLinker == ONLY_COMMON_ION)
						continue;
					bContainLinker = true;
				} else {
					if(stType.nHasLinker == ONLY_LINKER_ION)
						continue;
					bContainLinker = false;
				}

				MzTriple triple;
				triple.lfMz = (m_pBetaPep->m_bPepAAMass[nIdx] + (stType.nCharge - 1) * PROTON_MASS - stType.lfLoss) / stType.nCharge;
				triple.nMz = triple.lfMz * MZ_MULTIPLIER + 0.5;
				triple.yIonTypeOrder = yOrder;
				triple.yPepPosOrder1 = nIdx + m_pAlphaPep->m_nPepLen;
				triple.yNTerm1 = TRUE;
				triple.yAddToTag = stType.bContinuity;
				triple.yContainLinker = bContainLinker;
				vIonMz.push_back(triple);
			}

			for(size_t j = 0; j < m_pXLinkPeptideResult->m_stBetaPep.m_stPep.m_vMods.size(); ++j) {
				// for neutral loss
				int nId = m_pXLinkPeptideResult->m_stBetaPep.m_stPep.m_vMods[j].m_nId;
				int nSite = m_pXLinkPeptideResult->m_stBetaPep.m_stPep.m_vMods[j].m_nSite;
				for(size_t k = 0; k < m_vMods[nId].vMonoNeutralLossDiff.size(); ++k) {
					for(int nIdx = nSite; nIdx < m_pBetaPep->m_nPepLen-1; ++nIdx) {
						bool bContainLinker(false);
						if(nIdx >= nBetaSite) {
							if(stType.nHasLinker == ONLY_COMMON_ION)
								continue;
							bContainLinker = true;
						} else {
							if(stType.nHasLinker == ONLY_LINKER_ION)
								continue;
							bContainLinker = false;
						}

						MzTriple triple;
						double lfMz = m_pBetaPep->m_bPepAAMass[nIdx] - m_vMods[nId].vMonoNeutralLossDiff[k];
						if(lfMz < 0)
							continue;
						triple.lfMz = (lfMz + (stType.nCharge - 1) * PROTON_MASS - stType.lfLoss) / stType.nCharge;
						triple.nMz =  triple.lfMz * MZ_MULTIPLIER + 0.5;
						triple.yIonTypeOrder = yOrder;
						triple.yPepPosOrder1 = nIdx + m_pAlphaPep->m_nPepLen;
						triple.yNTerm1 = TRUE;
						triple.yAddToTag = FALSE;
						triple.yContainLinker = bContainLinker;
						vIonMz.push_back(triple);
					}
				}
			}
		}
	}
}

void MzCalculator::attachCTermIons(std::vector<MzTriple> &vIonMz, const IonType &stType, byte yOrder)
{
	if(m_pPeptide) {

	} else if(m_pPeptideResult) {
		for(int nIdx = 0; nIdx < m_pAlphaPep->m_nPepLen-1; ++nIdx) {
			MzTriple triple;
			triple.lfMz = (m_pAlphaPep->m_yPepAAMass[nIdx] + (stType.nCharge - 1) * PROTON_MASS - stType.lfLoss) / stType.nCharge;
			triple.nMz = triple.lfMz * MZ_MULTIPLIER + 0.5;
			triple.yIonTypeOrder = yOrder;
			triple.yAddToTag = TRUE;
			triple.yPepPosOrder1 = nIdx;
			vIonMz.push_back(triple);
		}

		// for neutral loss
		for(size_t j = 0; j < m_pPeptideResult->m_stPep.m_vMods.size(); ++j) {
			int nId = m_pPeptideResult->m_stPep.m_vMods[j].m_nId;
			int nSite = m_pPeptideResult->m_stPep.m_vMods[j].m_nSite;
			for(size_t k = 0; k < m_vMods[nId].vMonoNeutralLossDiff.size(); ++k) {
				for(int nIdx = nSite; nIdx < m_pAlphaPep->m_nPepLen-1; ++nIdx) {
					MzTriple triple;
					double lfMz = m_pAlphaPep->m_yPepAAMass[nIdx] - m_vMods[nId].vMonoNeutralLossDiff[k];
					if(lfMz < 0)
						continue;
					triple.lfMz = (lfMz + (stType.nCharge - 1) * PROTON_MASS - stType.lfLoss) / stType.nCharge;
					triple.nMz = triple.lfMz * MZ_MULTIPLIER + 0.5;
					triple.yIonTypeOrder = yOrder;
					triple.yPepPosOrder1 = nIdx;
					triple.yAddToTag = FALSE;
					vIonMz.push_back(triple);
				}
			}
		}
	} else if(m_pXLinkPeptideResult) {
		PeptideType eType = m_pXLinkPeptideResult->m_ePepType;
		int nAlphaSite(m_pXLinkPeptideResult->m_stAlphaPep.m_stLinkMod.m_nSite);
		int nBetaSite(m_pXLinkPeptideResult->m_stBetaPep.m_stLinkMod.m_nSite);

		// C terminal ions of alpha peptide
		for(int nIdx = 0; nIdx < m_pAlphaPep->m_nPepLen-1; ++nIdx) {
			bool bContainLinker(false);
			if (eType == PT_LOOP) {
				if (m_pAlphaPep->m_nPepLen-1-nIdx > nAlphaSite && m_pAlphaPep->m_nPepLen-1-nIdx <= nBetaSite)
					continue;
			} else if(eType == PT_XLINK) {
				if (m_pAlphaPep->m_nPepLen-1-nIdx <= nAlphaSite) { // contain linker
					if (stType.nHasLinker == ONLY_COMMON_ION) // only consider common ion
						continue;
					bContainLinker = true;
				} else { // not contain linker
					if (stType.nHasLinker == ONLY_LINKER_ION) // only consider x-link ion
						continue;
					bContainLinker = false;
				}
			}

			MzTriple triple;
			triple.lfMz = (m_pAlphaPep->m_yPepAAMass[nIdx] + (stType.nCharge - 1) * PROTON_MASS - stType.lfLoss) / stType.nCharge;
			triple.nMz = triple.lfMz * MZ_MULTIPLIER + 0.5;
			triple.yIonTypeOrder = yOrder;
			triple.yPepPosOrder1 = nIdx;
			triple.yNTerm1 = FALSE;
			triple.yAddToTag = stType.bContinuity;
			triple.yContainLinker = bContainLinker;
			vIonMz.push_back(triple);
		}

		for(size_t j = 0; j < m_pXLinkPeptideResult->m_stAlphaPep.m_stPep.m_vMods.size(); ++j) {
			// for neutral loss
			int nId = m_pXLinkPeptideResult->m_stAlphaPep.m_stPep.m_vMods[j].m_nId;
			int nSite = m_pXLinkPeptideResult->m_stAlphaPep.m_stPep.m_vMods[j].m_nSite;
			for(size_t k = 0; k < m_vMods[nId].vMonoNeutralLossDiff.size(); ++k) {
				for(int nIdx = nSite; nIdx < m_pAlphaPep->m_nPepLen-1; ++nIdx) {
					bool bContainLinker(false);
					if (eType == PT_LOOP) {
						if (m_pAlphaPep->m_nPepLen-1-nIdx > nAlphaSite && m_pAlphaPep->m_nPepLen-1-nIdx <= nBetaSite)
							continue;
					} else if(eType == PT_XLINK) {
						if (m_pAlphaPep->m_nPepLen-1-nIdx <= nAlphaSite) { // contain linker
							if (stType.nHasLinker == ONLY_COMMON_ION) // only consider common ion
								continue;
							bContainLinker = true;
						} else { // not contain linker
							if (stType.nHasLinker == ONLY_LINKER_ION) // only consider x-link ion
								continue;
							bContainLinker = false;
						}
					}

					MzTriple triple;
					double lfMz = m_pAlphaPep->m_yPepAAMass[nIdx] - m_vMods[nId].vMonoNeutralLossDiff[k];
					if(lfMz < 0)
						continue;
					triple.lfMz = (lfMz + (stType.nCharge - 1) * PROTON_MASS - stType.lfLoss) / stType.nCharge;
					triple.nMz = triple.lfMz * MZ_MULTIPLIER + 0.5;
					triple.yIonTypeOrder = yOrder;
					triple.yPepPosOrder1 = nIdx;
					triple.yNTerm1 = FALSE;
					triple.yAddToTag = FALSE;
					triple.yContainLinker = bContainLinker;
					vIonMz.push_back(triple);
				}
			}
		}

		// C terminal ions of beta peptide, if PT_XLINK is set on
		if(eType == PT_XLINK) {
			for(int nIdx = 0; nIdx < m_pBetaPep->m_nPepLen-1; ++nIdx) {
				bool bContainLinker(false);
				if(m_pBetaPep->m_nPepLen-1-nIdx <= nBetaSite) {
					if(stType.nHasLinker == ONLY_COMMON_ION)
						continue;
					bContainLinker = true;
				} else {
					if(stType.nHasLinker == ONLY_LINKER_ION)
						continue;
					bContainLinker = false;
				}

				MzTriple triple;
				triple.lfMz = (m_pBetaPep->m_yPepAAMass[nIdx] + (stType.nCharge - 1) * PROTON_MASS - stType.lfLoss) / stType.nCharge;
				triple.nMz = triple.lfMz * MZ_MULTIPLIER + 0.5;
				triple.yIonTypeOrder = yOrder;
				triple.yPepPosOrder1 = nIdx + m_pAlphaPep->m_nPepLen;
				triple.yNTerm1 = FALSE;
				triple.yAddToTag = stType.bContinuity;
				triple.yContainLinker = bContainLinker;
				vIonMz.push_back(triple);
			}

			for(size_t j = 0; j < m_pXLinkPeptideResult->m_stBetaPep.m_stPep.m_vMods.size(); ++j) {
				// for neutral loss
				int nId = m_pXLinkPeptideResult->m_stBetaPep.m_stPep.m_vMods[j].m_nId;
				int nSite = m_pXLinkPeptideResult->m_stBetaPep.m_stPep.m_vMods[j].m_nSite;
				for(size_t k = 0; k < m_vMods[nId].vMonoNeutralLossDiff.size(); ++k) {
					for(int nIdx = nSite; nIdx < m_pBetaPep->m_nPepLen-1; ++nIdx) {
						bool bContainLinker(false);
						if(m_pBetaPep->m_nPepLen-1-nIdx <= nBetaSite) {
							if(stType.nHasLinker == ONLY_COMMON_ION)
								continue;
							bContainLinker = true;
						} else {
							if(stType.nHasLinker == ONLY_LINKER_ION)
								continue;
							bContainLinker = false;
						}

						MzTriple triple;
						double lfMz = m_pBetaPep->m_yPepAAMass[nIdx] - m_vMods[nId].vMonoNeutralLossDiff[k];
						if(lfMz < 0)
							continue;
						triple.lfMz = (lfMz + (stType.nCharge - 1) * PROTON_MASS - stType.lfLoss) / stType.nCharge;
						triple.nMz = triple.lfMz * MZ_MULTIPLIER + 0.5;
						triple.yIonTypeOrder = yOrder;
						triple.yPepPosOrder1 = nIdx + m_pAlphaPep->m_nPepLen;
						triple.yNTerm1 = FALSE;
						triple.yAddToTag = FALSE;
						triple.yContainLinker = bContainLinker;
						vIonMz.push_back(triple);
					}
				}
			}
		}
	} else {
		// do nothing
	}
}

void MzCalculator::close()
{
	if(m_pPeptide)
		m_pPeptide = NULL;

	if(m_pPeptideResult)
		m_pPeptideResult = NULL;

	if(m_pXLinkPeptideResult)
		m_pXLinkPeptideResult = NULL;

	if(m_pAlphaPep) {
		delete m_pAlphaPep;
		m_pAlphaPep = NULL;
	}

	if(m_pBetaPep) {
		delete m_pBetaPep;
		m_pBetaPep = NULL;
	}
}

}
