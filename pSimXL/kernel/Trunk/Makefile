include $(CONFIG)

MAKE := make

all:
	mkdir -p $(BUILD_DIR)
	mkdir -p $(INSTALL_DIR)
	cp -r $(BUILD_LIST) $(BUILD_DIR)
	cp -r $(RUN_LIST) $(INSTALL_DIR)
	cd util && $(MAKE)
	cd io && $(MAKE)
	cd index && $(MAKE)
	cd score && $(MAKE)
	cd flow && $(MAKE)
	cd search && $(MAKE)
	cd searcher && $(MAKE)
	
install:
	cd util && $(MAKE) install
	cd io && $(MAKE) install
	cd index && $(MAKE) install
	cd score && $(MAKE) install
	cd flow && $(MAKE) install
	cd search && $(MAKE) install
	cd searcher && $(MAKE) install

clean:
	cd util && $(MAKE) clean
	cd io && $(MAKE) clean
	cd index && $(MAKE) clean
	cd score && $(MAKE) clean
	cd flow && $(MAKE) clean
	cd search && $(MAKE) clean
	cd searcher && $(MAKE) clean
	rm -fr $(BUILD_DIR)
	rm -fr $(INSTALL_DIR)
