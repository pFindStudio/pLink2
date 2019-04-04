#include "../include/flow.h"

using namespace std;

namespace sdk 
{

Flow::Flow() :
		m_pParameter(ParameterReader::getInstance()->getParameter())
{
	int spectra_total_num=0;
	for(int i=0;i<m_pParameter->m_tSpectraNo.size();++i){
		spectra_total_num+=m_pParameter->m_tSpectraNo[i];
	}
	m_vSpectra.resize(spectra_total_num);
}

Flow::~Flow()
{
}

void Flow::getPrecursorBorder(PrecursorWindow &stWindow, size_t tSpecIndex)
{
	double lfNeutralMass = m_vSpectra[tSpecIndex].m_lfMH - PROTON_MASS - H2O_MONO_MASS;
	double lfTol(0.0);
	if(m_pParameter->m_ePepTolType == TT_PPM) {
		lfTol = m_vSpectra[tSpecIndex].m_lfMZ * m_vSpectra[tSpecIndex].m_tCharge
				* m_pParameter->m_lfPepTol * PART_PER_MILLION;
	} else {
		lfTol = m_pParameter->m_lfPepTol;
	}
	stWindow.m_lfMinMass = lfNeutralMass - lfTol;
	stWindow.m_lfMaxMass = lfNeutralMass + lfTol;
}



bool Flow::isModLinkSiteConflict(const Peptide &stPep, size_t tIdx)
{
	for(size_t i = 0; i < stPep.m_vMods.size(); ++i) {
		if(stPep.m_vMods[i].m_nSite == (int)tIdx)
			return true;
	}
	return false;
}

FlowFactory::FlowFactory() :
        m_pParameter(ParameterReader::getInstance()->getParameter())
{
}

FlowFactory::~FlowFactory()
{
}

Flow *FlowFactory::getFlow(FlowType eType) const
{
    switch(eType) {
    case FT_SIMULATION:
     	return new BenchMarkFlow();
    default:
        ErrorInfo err("FlowFactory", "getFlow", "unkown flow type.");
        throw runtime_error(err.get());
        return NULL;
    }
}

} /* namespace sdk */

