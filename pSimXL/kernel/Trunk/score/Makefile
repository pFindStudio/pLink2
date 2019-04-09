
include $(CONFIG)

ifdef IS_LINUX
	SUFFIX = .so
	CFLAGS += -fPIC
endif

ifdef IS_WINDOWS
	SUFFIX = .dll
	LFLAGS += -Wl,--enable-auto-import -Wl,--enable-runtime-pseudo-reloc
endif

SRC    = MzCalculator.cpp RefinedScorer.cpp BMIonGenerator.cpp BMIonPerfectGenerator.cpp BMIonSimpleGenerator.cpp
OBJS   = $(SRC:.cpp=.o)
LIBS   = -lutil
TARGET = libscore$(SUFFIX)

all: $(TARGET)
ifeq ($(UNITTEST), 1)
	cd showIonTypes && $(MAKE) BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
endif

clean:
	rm -f $(OBJS) $(TARGET)
ifeq ($(UNITTEST), 1)
	cd showIonTypes && $(MAKE) clean BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
endif

install:
	cp $(TARGET) $(INSTALL_DIR)/$(TARGET)
ifeq ($(UNITTEST), 1)
	cd showIonTypes && $(MAKE) install BIT64=$(BIT64) UNITTEST=$(UNITTEST) BUILD_DIR=$(BUILD_DIR) INSTALL_DIR=$(INSTALL_DIR) CONFIG=$(CONFIG)
endif
	
$(OBJS): $(SRC)
	$(CXX) $(CFLAGS) -c $(SRC)
	 	
$(TARGET): $(OBJS)
	$(CXX) $(LFLAGS) -L$(BUILD_DIR) -shared -o $(TARGET) $(OBJS) $(LIBS) 
	cp $(TARGET) $(BUILD_DIR)/$(TARGET)
