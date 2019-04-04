#include "../include/searcher.h"

using namespace std;
using namespace sdk;

int main(int argc, char *argv[]) {
	try {
		if (argc == 2) {
			EnvironmentInitializer::init(argv[1]);
			SearchParameter *pParameter = EnvironmentInitializer::getSearchParameter();

			Searcher stSearcher;
			stSearcher.search();

			EnvironmentInitializer::destroy();
		} else {
			cerr << "Usage:" << endl;
			cerr << "pSimulate.exe config file" << endl;
		}
	} catch (exception &e) {
		ErrorInfo err("pSimulate", "main", "caught an exception.", e);
		cerr << err.get() << endl;
	} catch (...) {
		ErrorInfo err("pSimulate", "main", "caught an unknown exception.");
		cerr << err.get() << endl;
	}
	return 0;
}

