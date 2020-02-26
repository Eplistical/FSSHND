CXX = g++
SUBDIRS := $(wildcard */.)

all: libmisc.so

libmisc.so: $(SUBDIRS)
	$(CXX) -shared -o $@ */*.o
	

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	rm */*.o

.PHONY: all clean $(SUBDIRS)

