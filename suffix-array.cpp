#include <fstream>
#include <iostream>

#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/sequence.h>

using std::cout;

// compile with:
// g++ suffix-array.cpp -o suffix-array -O3

int main(int argc, char** argv) {
	if (argc != 2) {
		cout << "suffix-aray requires 1 argument (name of fasta file).";
		return 1;
	}
	std::fstream in(argv[1], std::ios::binary | std::ios::in);
	// seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(in);
	// read file one record at a time
	seqan::String<char> id;
	seqan::String<char> seq;
	while (!in.eof()) {
	// 	// if (readRecord(id, seq, reader) != 0) {
	// 	// 	return 1;  // error reading record from file
	// 	// }
		if (in.peek() == '>') {
			readMeta(in, id, seqan::Fasta());			
			read(in, seq, seqan::Fasta());
			cout << ">" << id << "\n";
			seqan::String<unsigned> sa;
			// build a suffix array using the Skew7 algorithm
			resize(sa, length(seq));
			createSuffixArray(sa, seq, seqan::Skew7());
			for (unsigned i=0; i<length(sa); i++) {
				cout << sa[i] << "\n";
			}
		}
	}
	in.close();
	return 0;
}
