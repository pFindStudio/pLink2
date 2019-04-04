#include "../include/util.h"

using namespace std;

namespace sdk
{

ConstantStore *ConstantStore::m_pInstance = NULL;

ConstantStore *ConstantStore::getInstance()
{
	if(m_pInstance == NULL) {
		m_pInstance = new ConstantStore();
	}
	return m_pInstance;
}

void ConstantStore::destroyInstance()
{
	if(m_pInstance) {
		delete m_pInstance;
		m_pInstance = NULL;
	}
}

ConstantStore::~ConstantStore()
{
}

ConstantStore::ConstantStore()
{
	try {
		char szBuffer[STR_BUF_SIZE];
		getcwd(szBuffer, STR_BUF_SIZE);
		m_strWorkDir = szBuffer;

		initAAMass();
		initLogPrime();
		initLogNPerm();
	} catch(exception &e) {
		ErrorInfo err("ConstantStore", "ConstantStore", "caught an exception in the constructor of ConstantStore.");
		throw runtime_error(err.get());
	}
}

void ConstantStore::initAAMass()
{
	AminoAcidDict *pDict = AminoAcidDict::getInstance();

	for(size_t i = 0; i < AA_NUM; i++) {
		m_lfAAMass[i] = pDict->getMonoMass('A' + i);
	}
}

void ConstantStore::initLogPrime()
{
	int nCurInteger = 2;
	for(size_t i = 0; i < LOG_PRIME_SIZE; ++i) {
		while(!isPrime(nCurInteger))
			++nCurInteger;
		m_lfLogPrime[i] = log(nCurInteger);
		++nCurInteger;
	}
}

bool ConstantStore::isPrime(int nNumber)
{
	for(int n = 2; n <= sqrt(nNumber); ++n)
		if(nNumber % n == 0)
			return false;
	return true;
}

void ConstantStore::initLogNPerm()
{
	double lfSum = 0.0;
	m_lfLogNPerm[0] = 0.0;
	for(size_t i = 1; i < LOG_N_SIZE; ++i) {
		lfSum += log(i);
		m_lfLogNPerm[i] = lfSum;
	}
}

}
