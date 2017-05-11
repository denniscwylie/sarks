CPP = g++
CPPFLAGS = -O3

all: suffix-array windginiimp

suffix-array: suffix-array.cpp
	$(CPP) suffix-array.cpp -o suffix-array $(CPPFLAGS)

windginiimp: windginiimp.cpp
	$(CPP) windginiimp.cpp -o windginiimp $(CPPFLAGS)

install: suffix-array windginiimp
	install -m 0755 suffix-array ${PREFIX}/bin
	install -m 0755 windginiimp ${PREFIX}/bin
	rm suffix-array
	rm windginiimp

clean:
	rm suffix-array
	rm windginiimp

.PHONY: all clean install
