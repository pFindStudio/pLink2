
#ifndef DECORATORS_H_
#define DECORATORS_H_
#include "util.h"
#include "sdk.h"

namespace sdk
{

template<typename Invoke, typename Func>
class ModificationDecorator {
protected:
	const SearchParameter *m_pParameter;
	std::map<char, std::vector<size_t> > m_mpFixedMods;
	std::map<char, std::vector<size_t> > m_mpVarMods;
	double m_lfMaxMod;
	double m_lfLeastNegativeMod;

public:
	ModificationDecorator(const SearchParameter *pParameter) :
		m_pParameter(pParameter)
	{
		initModifications();
	}

	virtual ~ModificationDecorator() {}

	void addModifications(Peptide &stPep, Invoke *pFlow, Func func)
	{
		// fixed modification
		size_t tFixedModNum = 0;

		if (stPep.m_strSq.length() <= 0)
			return;

		ModificationDict *pDict = ModificationDict::getInstance();
		const std::vector<std::string> &vVarMods = m_pParameter->m_vVariableMods;
		const std::vector<std::string> &vFixedMods = m_pParameter->m_vFixedMods;
		// fixed modification on N-terminal of protein
		for (size_t i = 0; i < m_mpFixedMods[N_TERM_MOD].size(); ++i) {
			size_t tIdx = m_mpFixedMods[N_TERM_MOD][i] - vVarMods.size();
			const Modification &mod = pDict->getModification(vFixedMods[tIdx]);

			if(mod.eModType == MT_PRO_N && !stPep.isProNTerm())
				continue;

			// take modification site into consideration
			if (mod.strAA.find(stPep.m_strSq[0]) == std::string::npos)
				continue;

			// notice when scoring psm, fixed modifications are attached after variable modifications
			ModificationEntry stEntry(0, m_mpFixedMods[N_TERM_MOD][i]);
			stPep.m_vMods.push_back(stEntry);
			stPep.m_lfMass += mod.lfMonoDiff;
			++tFixedModNum;
		}

		// fixed modification on C-terminal of protein
		for (size_t i = 0; i < m_mpFixedMods[C_TERM_MOD].size(); ++i) {
			size_t tIdx = m_mpFixedMods[C_TERM_MOD][i] - vVarMods.size();
			const Modification &mod = pDict->getModification(vFixedMods[tIdx]);

			if(mod.eModType == MT_PRO_C && !stPep.isProCTerm())
				continue;

			// take modification site into consideration
			if (mod.strAA.find(stPep.m_strSq[stPep.m_strSq.length()-1]) == std::string::npos)
				continue;

			ModificationEntry stEntry(stPep.m_strSq.length()-1, m_mpFixedMods[C_TERM_MOD][i]);
			stPep.m_vMods.push_back(stEntry);
			stPep.m_lfMass += mod.lfMonoDiff;
			++tFixedModNum;
		}

		// fixed modification on the body of peptide
		for (size_t i = 0; i < stPep.m_strSq.length(); ++i) {
			char cSite = stPep.m_strSq[i];
			for (size_t j = 0; j < m_mpFixedMods[cSite].size(); ++j) {
				size_t tIdx = m_mpFixedMods[cSite][j] - vVarMods.size();
				const Modification &mod = pDict->getModification(vFixedMods[tIdx]);

				ModificationEntry stEntry(i, m_mpFixedMods[cSite][j]);
				stPep.m_vMods.push_back(stEntry);
				stPep.m_lfMass += mod.lfMonoDiff;
				++tFixedModNum;
			}
		}

		// variable modification
		addVariableMods(stPep, -1, tFixedModNum, pFlow, func);
	}
	double getMaximumMod() const
	{
		return m_lfMaxMod;
	}
	double getLeastNegativeMod() const
	{
		return m_lfLeastNegativeMod;
	}

private:
	void initModifications()
	{
		m_mpFixedMods.clear();
		m_mpVarMods.clear();
		m_lfLeastNegativeMod = 0.0;
		m_lfMaxMod = 0.0;

		ModificationDict *pDict = ModificationDict::getInstance();
		const std::vector<std::string> &vVarMods = m_pParameter->m_vVariableMods;
		for(size_t t = 0; t < vVarMods.size(); ++t) {
			const Modification &mod = pDict->getModification(vVarMods[t]);
			if(mod.eModType == MT_PEP_N || mod.eModType == MT_PRO_N) {
				m_mpVarMods[N_TERM_MOD].push_back(t);
			} else if(mod.eModType == MT_PEP_C || mod.eModType == MT_PRO_N) {
				m_mpVarMods[C_TERM_MOD].push_back(t);
			} else {
				for(size_t i = 0; i < mod.strAA.length(); ++i) {
					m_mpVarMods[mod.strAA[i]].push_back(t);
				}
			}

			// initialize least negative modification
			if(mod.lfMonoDiff < m_lfLeastNegativeMod)
				m_lfLeastNegativeMod = mod.lfMonoDiff;

			// initialize maximum modification
			if(mod.lfMonoDiff > m_lfMaxMod)
				m_lfMaxMod = mod.lfMonoDiff;
		}
		const std::vector<std::string> &vFixedMods = m_pParameter->m_vFixedMods;
		for(size_t t = 0; t < vFixedMods.size(); ++t) {
			const Modification &mod = pDict->getModification(vFixedMods[t]);
			if(mod.eModType == MT_PEP_N || mod.eModType == MT_PRO_N) {
				m_mpFixedMods[N_TERM_MOD].push_back(t + vVarMods.size());
			} else if(mod.eModType == MT_PEP_C || mod.eModType == MT_PRO_N) {
				m_mpFixedMods[C_TERM_MOD].push_back(t + vVarMods.size());
			} else {
				for(size_t i = 0; i < mod.strAA.length(); ++i) {
					m_mpFixedMods[mod.strAA[i]].push_back(t + vVarMods.size());
				}
			}

			// initialize least negative modification
			if(mod.lfMonoDiff < m_lfLeastNegativeMod)
				m_lfLeastNegativeMod = mod.lfMonoDiff;

			// initialize maximum modification
			if(mod.lfMonoDiff > m_lfMaxMod)
				m_lfMaxMod = mod.lfMonoDiff;
		}
	}
	void addVariableMods(Peptide &stPep, int nCurPos, size_t tFixedModNum, Invoke *pFlow, Func func)
	{
		if (stPep.m_strSq.length() <= 0)
			return;

		if ((stPep.m_vMods.size() - tFixedModNum) > m_pParameter->m_tMaxModsNo)
			return;

		ModificationDict *pDict = ModificationDict::getInstance();
		const std::vector<std::string> &vVarMods = m_pParameter->m_vVariableMods;

		if (nCurPos == -1) {
			// variable modification on N-terminal of protein
			// if without
			addVariableMods(stPep, nCurPos+1, tFixedModNum, pFlow, func);
			// if with
			if ((stPep.m_vMods.size() - tFixedModNum) < m_pParameter->m_tMaxModsNo) {
				for (size_t i = 0; i < m_mpVarMods[N_TERM_MOD].size(); ++i) {
					size_t tIdx = m_mpVarMods[N_TERM_MOD][i];
					const Modification &mod = pDict->getModification(vVarMods[tIdx]);

					if (mod.eModType == MT_PRO_N && !stPep.isProNTerm())
						continue;

					// take modification site into consideration
					if (mod.strAA.find(stPep.m_strSq[0]) == std::string::npos)
						continue;

					ModificationEntry stEntry(nCurPos+1, m_mpVarMods[N_TERM_MOD][i]);
					stPep.m_lfMass += mod.lfMonoDiff;
					stPep.m_vMods.push_back(stEntry);
					addVariableMods(stPep, nCurPos+1, tFixedModNum, pFlow, func);
					stPep.m_vMods.pop_back();
					stPep.m_lfMass -= mod.lfMonoDiff;
				}
			}
		} else if (nCurPos == (int)stPep.m_strSq.length()) {
			// variable modification on C-terminal of protein
			// if without
			(pFlow->*func)();
			// if with
			if ((stPep.m_vMods.size() - tFixedModNum) < m_pParameter->m_tMaxModsNo) {
				for (size_t i = 0; i < m_mpVarMods[MT_PRO_C].size(); ++i) {
					size_t tIdx = m_mpVarMods[C_TERM_MOD][i];
					const Modification &mod = pDict->getModification(vVarMods[tIdx]);

					if (mod.eModType == MT_PRO_C && !stPep.isProCTerm())
						continue;

					// take modification site into consideration
					if (mod.strAA.find(stPep.m_strSq[0]) == std::string::npos)
						continue;

					ModificationEntry stEntry(nCurPos-1, m_mpVarMods[C_TERM_MOD][i]);
					stPep.m_lfMass += mod.lfMonoDiff;
					stPep.m_vMods.push_back(stEntry);
					(pFlow->*func)();
					stPep.m_vMods.pop_back();
					stPep.m_lfMass -= mod.lfMonoDiff;
				}
			}
		} else {
			while (nCurPos < (int)stPep.m_strSq.length()
					&& m_mpVarMods.find(stPep.m_strSq[nCurPos]) == m_mpVarMods.end()) {
				++nCurPos;
			}

			if (nCurPos == (int)stPep.m_strSq.length())
				addVariableMods(stPep, nCurPos, tFixedModNum, pFlow, func);
			else {
				// var mod on body of peptide
				// if without
				addVariableMods(stPep, nCurPos+1, tFixedModNum, pFlow, func);
				// if with
				char cSite = stPep.m_strSq[nCurPos];
				if ((stPep.m_vMods.size() - tFixedModNum) < m_pParameter->m_tMaxModsNo) {
					for (size_t i = 0; i < m_mpVarMods[cSite].size(); ++i) {
						size_t tIdx = m_mpVarMods[cSite][i];
						const Modification &mod = pDict->getModification(vVarMods[tIdx]);

						ModificationEntry stEntry(nCurPos, m_mpVarMods[cSite][i]);
						stPep.m_lfMass += mod.lfMonoDiff;
						stPep.m_vMods.push_back(stEntry);
						addVariableMods(stPep, nCurPos+1, tFixedModNum, pFlow, func);
						stPep.m_vMods.pop_back();
						stPep.m_lfMass -= mod.lfMonoDiff;
					}
				}
			}
		}
	}
};

class XLinkDecorator {
protected:
	const SearchParameter *m_pParameter;
	bool m_bMonoLinked[MAX_LINKER_NUM][LINK_SITE_NUM];
	bool m_bLinked[MAX_LINKER_NUM][LINK_SITE_NUM][LINK_SITE_NUM];
	double m_lfLeastNegativeLinker;

public:
	XLinkDecorator(const SearchParameter *m_pParameter);
	virtual ~XLinkDecorator() {}

	bool isXLinked(const Peptide &stPep);
	bool isXLinkSiteLegal(const Peptide &stPep, size_t tIdx, size_t tLinkerId);
	bool isXLinkSiteLegal(const Peptide &stPep1, size_t tIdx1,
			const Peptide &stPep2, size_t tIdx2, size_t tLinkerId);
	double getLeastNegativeLinker() const;
	double getLinkerMass(bool bMono, size_t tLinkerId);

private:
	void initLinkers();
	void setOnLinkers();
	size_t getCandXLinkSiteIndex(const Peptide &stPep, size_t tIdx, size_t tLinkerId);
};

}

#endif /* DECORATORS_H_ */
