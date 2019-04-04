#include "../include/io.h"
#include "../include/index.h"
#include "../include/searcher.h"

using namespace std;

namespace sdk
{

Searcher::Searcher()
{
}

Searcher::~Searcher()
{
}

void Searcher::search()
{
	const SearchParameter *pParameter = ParameterReader::getInstance()->getParameter();

	FlowFactory stFlowFactory;
	Flow *flow = stFlowFactory.getFlow(pParameter->m_eFlowType);

	flow->init();
	flow->run();

	delete flow;
}



}


