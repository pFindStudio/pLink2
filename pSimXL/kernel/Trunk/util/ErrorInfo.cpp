#include "../include/util.h"

using namespace std;

namespace sdk
{

ErrorInfo::ErrorInfo(const string &strClass, const string &strMethod,
		const string &strDetail)
{
	setErrorClass(strClass);
	setErrorMethod(strMethod);
	setErrorDetail(strDetail);
}

ErrorInfo::ErrorInfo(const string &strClass, const string &strMethod,
		const string &strDetail, const exception &e)
{
	setErrorClass(strClass);
	setErrorMethod(strMethod);
	setErrorDetail(strDetail);
	setException(e);
}

void ErrorInfo::setErrorClass(const string &strClass)
{
	m_strClass = strClass;
}

void ErrorInfo::setErrorMethod(const string &strMethod)
{
	m_strMethod = strMethod;
}

void ErrorInfo::setErrorDetail(const string &strDetail)
{
	m_strDetail = strDetail;
}

void ErrorInfo::setException(const exception &e)
{
	m_strException = e.what();
}

string ErrorInfo::get() const
{
	string strError = m_strException;
	strError += "\t at " + m_strClass + "::" + m_strMethod + "() "
			+ m_strDetail + "\n";
	return strError;
}

string ErrorInfo::get(const exception &e)
{
	setException(e);
	return get();
}

ostream &operator<<(ostream &os, const ErrorInfo &info)
{
	time_t stCurTime;

	time(&stCurTime);
	os<<endl<<"=========================="<<endl
			<<ctime(&stCurTime)<<endl
			<<info.get()<<endl;
	return os;
}

ofstream &operator<<(ofstream &os, const ErrorInfo &info)
{
	time_t stCurTime;

	time(&stCurTime);
	os<<endl<<"=========================="<<endl
			<<ctime(&stCurTime)<<endl
			<<info.get()<<endl;
	return os;
}

}
