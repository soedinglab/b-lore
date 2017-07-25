CC := g++
LIBFLAGS := -O2 -shared -fPIC

all: margloglik.so zstates.so

margloglik.so: margloglik.cpp
	$(CC) $< $(LIBFLAGS) -o $@

zstates.so: zstates_py.cpp
	$(CC) $< $(LIBFLAGS) -o $@

clean:
	@echo " Cleaning..."
	@rm -rf zstates.so margloglik.so
