
include $(CONFIG)

ifdef IS_LINUX
	SUFFIX = .so
	CFLAGS += -fPIC
	ifeq ($(BIT64), 0)
		PTHREAD = -lpthread
	endif
endif

ifdef IS_WINDOWS
	SUFFIX = .dll
	LFLAGS += -Wl,--enable-auto-import -Wl,--enable-runtime-pseudo-reloc
	ifeq ($(BIT64), 0)
		PTHREAD = -lpthreadGCE2
	endif
endif

SRC    = AccessionIndexer.cpp AminoAcidDict.cpp ConfigParser.cpp ConstantStore.cpp Decorators.cpp Dictionaries.cpp \
ElementDict.cpp ErrorInfo.cpp FastaReader.cpp ModificationDict.cpp ParameterReader.cpp sdk.cpp \
Trace.cpp WindowsHash.cpp \
NoSectionConfigParser.cpp OneByOneConfigParser.cpp EnvironmentInitializer.cpp Bias.cpp
OBJS   = $(SRC:%.cpp=%.o)
LIBS   = $(PTHREAD)
TARGET = libutil$(SUFFIX)

all: $(TARGET)
ifeq ($(UNITTEST), 1)
	cd testNoSectionConfigParser && $(MAKE) BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
	cd testElementDict && $(MAKE) BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
	cd testAADict && $(MAKE) BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
	cd testSumoDict && $(MAKE) BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
	cd testACIndexer && $(MAKE) BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
endif

clean:
	-rm -f $(OBJS) $(TARGET)
ifeq ($(UNITTEST), 1)
	cd testNoSectionConfigParser && $(MAKE) clean BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
	cd testElementDict && $(MAKE) clean BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
	cd testAADict && $(MAKE) clean BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
	cd testSumoDict && $(MAKE) clean BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
	cd testACIndexer && $(MAKE) clean BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
endif

install:
	cp $(TARGET) $(INSTALL_DIR)/$(TARGET)
ifeq ($(UNITTEST), 1)
	cd testNoSectionConfigParser && $(MAKE) install BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
	cd testElementDict && $(MAKE) install BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
	cd testAADict && $(MAKE) install BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
	cd testSumoDict && $(MAKE) install BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
	cd testACIndexer && $(MAKE) install BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
endif
	
$(OBJS): $(SRC)
	$(CXX) -c $(CFLAGS) $(SRC)
	 	
$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(LFLAGS) -shared $(OBJS) -L$(BUILD_DIR) $(LIBS) 
	cp $(TARGET) $(BUILD_DIR)/$(TARGET)
