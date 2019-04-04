include $(CONFIG)

ifdef IS_LINUX
	SUFFIX = .so
	CFLAGS += -fPIC
endif

ifdef IS_WINDOWS
	SUFFIX = .dll
	LFLAGS += -Wl,--enable-auto-import -Wl,--enable-runtime-pseudo-reloc
endif

SRC    = FlowFactory.cpp BenchMarkFlow.cpp
OBJS   = $(SRC:.cpp=.o)
LIBS   = -lutil -lindex -lscore -lio
TARGET = libflow$(SUFFIX)

all: $(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)

install:
	cp $(TARGET) $(INSTALL_DIR)/$(TARGET)
	
$(OBJS): $(SRC)
	$(CXX) $(CFLAGS) -c $(SRC)
	 	
$(TARGET): $(OBJS)
	$(CXX) $(LFLAGS) -L$(BUILD_DIR) -shared -o $(TARGET) $(OBJS) $(LIBS) 
	cp $(TARGET) $(BUILD_DIR)/$(TARGET)
