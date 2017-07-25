CC := g++
LIBFLAGS := -O2 -shared -fPIC
LIBDIR := blore/lib

all: $(LIBDIR)/margloglik.so $(LIBDIR)/zstates.so

$(LIBDIR)/margloglik.so: $(LIBDIR)/margloglik.cpp
	$(CC) $< $(LIBFLAGS) -o $@

$(LIBDIR)/zstates.so: $(LIBDIR)/zstates_py.cpp
	$(CC) $< $(LIBFLAGS) -o $@

clean:
	@echo " Cleaning..."
	@rm -rf $(LIBDIR)/zstates.so $(LIBDIR)/margloglik.so
