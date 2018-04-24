#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "dna.h"

using namespace std;

int main(int argc, char *argv[]) {
	if (argc < 2) {
		cerr << "Wrong number of arguments" << endl;
		exit(EXIT_FAILURE);
	}
	int cols = 80;
	if (argc == 3) {
		cols = atoi(argv[2]);
	}

	string filename = argv[1];
	ifstream infile(filename);
	DNA::DNA d;
	if (!infile.fail()) {
		d = DNA::DNA(infile);
	}

	cout << d.revcomp().toFasta(cols);

	DNA::DNA d2 = d.revcomp().revcomp();
	if (d2 == d) {cout << "Equality works" << endl;}
	d = DNA::DNA(">test", "GATTTACT");

	string query = "GTA";
	int idx = d.find(query, 0);
	cout << "found query at " << idx << endl;
	cout << d.getSequence() << endl;
	cout << d.revcomp().getSequence() << endl;
}