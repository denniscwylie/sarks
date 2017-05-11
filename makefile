CPP = g++
CPPFLAGS = -O3
# INC = -I/work/03319/dwylie/seqan/seqan-library-1.4.2/include

all: suffix-array windginiimp

suffix-array: suffix-array.cpp
	$(CPP) ${INC} suffix-array.cpp -o suffix-array $(CPPFLAGS)

windginiimp: windginiimp.cpp
	$(CPP) windginiimp.cpp -o windginiimp $(CPPFLAGS)

install: suffix-array windginiimp
	$(eval PREFIX ?= ~)
	install -m 0755 suffix-array ${PREFIX}/bin/
	install -m 0755 windginiimp ${PREFIX}/bin/
	rm suffix-array
	rm windginiimp

clean:
	rm suffix-array
	rm windginiimp

.PHONY: all clean install
