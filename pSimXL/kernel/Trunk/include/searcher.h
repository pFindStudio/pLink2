#ifndef SEARCHER_H_
#define SEARCHER_H_
#include "util.h"
#include "sdk.h"
#include "io.h"
#include "index.h"
#include "scorers.h"
#include "flow.h"

namespace sdk
{

class Searcher {

public:
	Searcher();
	virtual ~Searcher();

	virtual void search();
};

}


#endif /* IONINDEXSEARCHER_H_ */
