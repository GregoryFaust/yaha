# Determine the yaha build number
BUILDNUM := 83
# Decide if we are compiling for user mode, or internal lab usage.
# This almost exclusively determines the options available to the user when running YAHA.
# But also controls timing information.
USERMODE := TRUE

# Set up the sub-directory names.
SDIR := src
ODIR := obj
BDIR := bin

CC	     ?= gcc
CXX      ?= g++
CFLAGS   ?= -Wall -O3
CXXFLAGS ?= -Wall -O3
CPPFLAGS ?=
LDFLAGS  ?=

CFLAGS	 += -MMD -MP -std=gnu99
CXXFLAGS += -MMD -MP
LDFLAGS  += -pthread

# Set up flags depending on mode
ifdef USERMODE
PROG	 := yaha
CPPFLAGS += -D COMPILE_USER_MODE -D BUILDNUM=$(BUILDNUM)
else
PROG	 := yaha$(BUILDNUM)
CPPFLAGS += -D BUILDNUM=$(BUILDNUM)
CFLAGS   += -g
CXXFLAGS += -g
endif

# The list of object files.
OBJfiles := Main.o AlignArgs.o AlignHelpers.o AlignExtFrag.o AlignOutput.o BaseSeq.o Compress.o \
	FileHelpers.o GraphPath.o Index.o Query.o QueryMatch.o QueryState.o SW.o Math.o
# Convert to use the obj dir
OBJS := $(patsubst %,$(ODIR)/%,$(OBJfiles))

# Make everything, all of which is defined below
YAHA: $(ODIR) $(BDIR) $(BDIR)/$(PROG)

# Make the directories
$(ODIR):
	mkdir -p $(ODIR)
$(BDIR):
	mkdir -p $(BDIR)

# Link the program
$(BDIR)/$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS)

# Include the dependencies.
# The actual objects will be built by the below generic rules based on these dependencies.
# The - suppresses warning messages about missing files when building from scratch.
-include $(OBJS:.o=.d)

# Make the object files.
# The built in rules will miss the subdirectories.
$(ODIR)/%.o: $(SDIR)/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -Rf $(ODIR)
	rm -Rf $(BDIR)
