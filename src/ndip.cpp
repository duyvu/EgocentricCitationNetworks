//============================================================================
// Name        : ndip.cpp
// Author      : Duy Vu
// Version     :
// Copyright   : Free to Use
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "icmlArXivHepTh.h"

int main(int argc, char *argv[]) {

	trainTiedEgoCentricNetworkModel(atoi(argv[1]), atof(argv[2]),
			atof(argv[3]), argv[4]);

	//testTiedEgoCentricNetworkModel(atoi(argv[1]), atof(argv[2]), atof(argv[3]),
	//		argv[4]);

	return 0;
}
