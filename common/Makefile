include ../common.mk

all: bin make.dep libcommon

include ../make.rules

LIBOBJS = 
BINOBJS = $(addprefix bin/, $(LIBOBJS))

libcommon: bin/libcommon.$(LIB_EXT)

bin/libcommon.$(LIB_EXT): $(BINOBJS)
	$(LIB_TOOL) $(LIB_CREATE) $@ $(BINOBJS) 
