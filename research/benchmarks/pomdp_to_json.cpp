#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>

#include <iostream>

#include "Pomdp.h"
#include "zmdpMainConfig.h"
#include "zmdpMainConfig.cc" // embed default config file

using namespace std;

void usage(void) {
  cerr <<
    "usage: pomdp_to_json <foo.POMDP>\n"
		"  This writes the POMDP data in a form that can be read eg. by python via\n"
		"  ast.literal_eval(<output>).\n"
    "  -h or --help   Display this help\n";
  exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]) {
  int argi;
  char *pomdpFileName = 0;
  for (argi=1; argi < argc; argi++) {
    if (0 == strcmp("-h",argv[argi]) || 0 == strcmp("--help",argv[argi])) {
      usage();
    } else if (0 == pomdpFileName) {
      pomdpFileName = argv[argi];
    } else {
      cerr << "too many arguments" << endl;
      usage();
    }
  }
  if ( 0 == pomdpFileName ) {
    usage();
  }

  zmdp::ZMDPConfig* config = new zmdp::ZMDPConfig();
  config->readFromString("<defaultConfig>", defaultConfig.data);

  zmdp::Pomdp p(std::string(pomdpFileName), config);

	cout << "{\"Type\":\"POMDP\""
			 << ",\"dim\":" << p.getBeliefSize()
			 << ",\"nU\":"  << p.getNumActions()
			 << ",\"nO\":"   << p.getNumObservations()
			 << ",\"discount\":" << p.getDiscount();
	cout << ",\"R(s,a)\":";

  for (int s=0; s < p.getBeliefSize(); s++) {
		if (s == 0) printf("["); else printf(",");
		for (int a=0; a < p.getNumActions(); a++) {
      if (a == 0) printf("["); else printf(",");
			printf("%g", p.R(s,a));
    }
		printf("]");
  }
	printf("]");
	
	cout << ",\"b(s)\":";
  for (int s=0; s < p.getBeliefSize(); s++) {
		if (s == 0) printf("["); else printf(",");
		printf("%g", p.getInitialBelief()(s));
	}
	printf("]");

	cout << ",\"T(a,s,sp)\":";
  for (int a=0; a < p.getNumActions(); a++) {
		if (a == 0) printf("["); else printf(",");
    for (int s=0; s < p.getBeliefSize(); s++) {
			if (s == 0) printf("["); else printf(",");
      for (int sp=0; sp < p.getBeliefSize(); sp++) {
				if (sp == 0) printf("["); else printf(",");
				printf("%g", p.T[a](s,sp));
      }
			printf("]");
		}
		printf("]");
	}
	printf("]");

	cout << ",\"O(a,sp,o)\":";
  for (int a=0; a < p.getNumActions(); a++) {
		if (a == 0) printf("["); else printf(",");
    for (int sp=0; sp < p.getBeliefSize(); sp++) {
			if (sp == 0) printf("["); else printf(",");     
			for (int o=0; o < p.getNumObservations(); o++) {
				if (o == 0) printf("["); else printf(",");
				printf("%g", p.O[a](sp,o));
      }
			printf("]");
    }
		printf("]");
  }
	printf("]");

	cout << "}" << endl;
}
