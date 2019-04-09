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

SRC = Searcher.cpp  
LIBS   = -lutil -lio -lindex -lflow $(PTHREAD)
TARGET = libsearch$(SUFFIX)

OBJS   = $(SRC:.cpp=.o)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)

install:
	cp $(TARGET) $(INSTALL_DIR)/$(TARGET)
	
$(OBJS): $(SRC)
	$(CXX) $(CFLAGS) -c $(SRC)
	 	
$(TARGET): $(OBJS)
	$(CXX) $(LFLAGS) -L$(BUILD_DIR) -shared -o $(TARGET) $(OBJS) $(LIBS) 
	cp $(TARGET) $(BUILD_DIR)/$(TARGET)
	
