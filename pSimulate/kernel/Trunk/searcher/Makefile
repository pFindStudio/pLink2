include $(CONFIG)

ifdef IS_LINUX
	SUFFIX = 
	CFLAGS += -fPIC
	ifeq ($(BIT64), 0)
		PTHREAD = -lpthread
	endif
endif

ifdef IS_WINDOWS
	SUFFIX = .exe
	LFLAGS += -Wl,--enable-auto-import -Wl,--enable-runtime-pseudo-reloc
	ifeq ($(BIT64), 0)
		PTHREAD = -lpthreadGCE2
	endif
endif

SRC = main.cpp 
LIBS   = -lutil -lio -lindex -lflow -lsearch $(PTHREAD)
TARGET = pSimulate$(SUFFIX)

OBJS   = $(SRC:.cpp=.o)

all:	$(TARGET)
ifeq ($(UNITTEST), 1)
endif

clean:
	rm -f $(OBJS) $(TARGET)
ifeq ($(UNITTEST), 1)
endif

install:
	cp $(TARGET) $(INSTALL_DIR)/$(TARGET)
ifeq ($(UNITTEST), 1)
endif
	
$(OBJS): $(SRC)
	$(CXX) $(CXXFLAGS) -c $(SRC)

$(TARGET):	$(OBJS)
	$(CXX) $(CXXFLAGS) $(LINK_FLAGS) -L$(BUILD_DIR) -o $(TARGET) $(OBJS) $(LIBS)
	cp $(TARGET) $(BUILD_DIR)/$(TARGET)
	
