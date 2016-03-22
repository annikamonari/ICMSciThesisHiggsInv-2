CC := g++ # Main compiler
ROOTSYS := /Applications/root_v5.34.34/
DYLD_LIBRARY_PATH := $ROOTSYS/lib

SRCDIR := src
BUILDDIR := build
TARGET := bin/main # Main executable of project-- type 'make' and this gets built.
 
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT)) # dynamic, includes any file with .cpp extension

# OBJECTS - build object list based on available sources
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o)) 
CFLAGS := -g -Wall $(inherited) `root-config --cflags`
LDFLAGS := $(inherited) `root-config --libs --glibs  --cflags` -lTMVA -lMinuit -lXMLIO 
INC := -I include # ensures all header files are acessible

# use echos to make sure commands are running

$(TARGET): $(OBJECTS)
	@echo "Linking..."
	$(CC) $(LDFLAGS) $^ -o $(TARGET)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@echo "Compiling..."
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo "Cleaning..." 
	$(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: clean