#include "../include/util.h"

using namespace std;

namespace sdk
{

WindowsHash::WindowsHash()
: m_pData(NULL), m_nStart(0), m_lfMaxMass(0.0)
{
}

WindowsHash::~WindowsHash()
{
	if(m_pData != NULL) {
		delete [] m_pData;
		m_pData = NULL;
	}
}

void WindowsHash::construct(const vector<PrecursorWindow> &vInfo)
{
	if(vInfo.empty())
		return;

	m_nStart = (int)vInfo[0].m_lfMaxMass;
	size_t tSize = (int)vInfo[vInfo.size() - 1].m_lfMaxMass - m_nStart + 10;
	m_pData = new int[tSize];
	memset(m_pData, 0, sizeof(int) * tSize);
	m_lfMaxMass = 0.0;
	int nPreIdx = -1;
	for(int i = 0; i < (int)vInfo.size(); ++i) {
		int nCurIdx = (int)vInfo[i].m_lfMaxMass - m_nStart;
		if(nCurIdx > nPreIdx) {
			for(int j = nCurIdx; j > nPreIdx; j--)
				m_pData[j] = i;
			nPreIdx = nCurIdx;
		}
		if(vInfo[i].m_lfMaxMass > m_lfMaxMass)
			m_lfMaxMass = vInfo[i].m_lfMaxMass;
	}
}

bool WindowsHash::isBeyondBound(const vector<PrecursorWindow> &vInfo, double lfMass)
{
	if(lfMass<vInfo[0].m_lfMinMass || lfMass>m_lfMaxMass)
		return true;
	else
		return false;
}

int WindowsHash::findLowerBound(double lfMass)
{
	if((int)lfMass - m_nStart >= 0)
		return m_pData[(int)lfMass - m_nStart];
	else
		return 0;
}

}
