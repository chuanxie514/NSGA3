

#include "global.h"
#include "nsga3.h"

int main() {

	size_t num_exp = 3;	// 总实验次数
	numVariables = 12; 	// 变量的数量，此处为12
	strcpy_s(strTestInstance, "DTLZ2");

	ifstream indata("TestNSGA3.txt");
	if (!indata) {
		cerr << "Error: file could not be opened" << endl;
		exit(1);
	}

	char temp[1024];
	while (!indata.eof()) {
		indata >> temp >> numObjectives;
		indata >> temp >> max_gen;
		indata >> temp >> p_boundary;
		indata >> temp >> p_inside;

		for (size_t run = 1; run <= num_exp; ++run) {
			printf(" Running experiment %lu for %s problem with %d objectives\n",
				run, strTestInstance, numObjectives);

			seed = (seed + 111) % 1235;
			rnd_uni_init = -(long)seed;

			TNSGA3 NSGA3;
			NSGA3.run(max_gen, run);
		}
	}
	indata.close();
	return 1;
}
