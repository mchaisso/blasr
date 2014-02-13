#
#
#
SHELL = /bin/bash
GCC = g++
SRC_DIR = src
BUILD_DIR = build
TARGET = jabon

CCOPT = -g
LIBS = -lboost_program_options -lpthread -L$(BUILD_DIR) -lsoap
DEFINES = -DDB_SHORT -DREAD_36 -DMAXGAP=3 -DMAXHITS=10000
INCLUDES = -Isrc/libsoap
ifneq ($(shell uname -s),Darwin)
    STATIC   = -static
endif

OBJS = \
  $(SRC_DIR)/jabon.o 			\
  $(SRC_DIR)/JabonMapper.o 		\
  $(SRC_DIR)/SoapShortHit.o 	\
  $(SRC_DIR)/SoapHitChain.o  	\
  $(SRC_DIR)/SoapMappingHit.o \
  $(SRC_DIR)/Rectangle.o 		\
  $(SRC_DIR)/RectangleTree.o 		\
  $(SRC_DIR)/GusfieldChainer.o 	

all: init libsoap $(TARGET)

init:
	@( [[ ! -e $(BUILD_DIR) ]] && mkdir $(BUILD_DIR) ) || true

$(TARGET): $(OBJS)
	$(GCC) $(DEFINES) $(INCLUDES) $(CCOPT) $(OBJS) $(STATIC) -o $(BUILD_DIR)/$@ $(LIBS)

%.o: %.cpp
	$(GCC) $(DEFINES) $(INCLUDES) $(CCOPT) -c $< -o $@

libsoap:
	make -C src/libsoap
		
INSTALL_DIR ?= .
install:
	/usr/bin/install -m 555 $(BUILD_DIR)/$(TARGET) $(INSTALL_DIR)/
	

clean:
	rm -rf $(BUILD_DIR)
	rm -f $(OBJS)
	make -C src/libsoap clean
