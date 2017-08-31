LIBRARY_DIR=path/to/itensor
HEADERS=itebd.h

#################################################################

include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk


cond-mat.0605597: %: %.cpp $(HEADERS)
	$(CCCOM) $< -o $@ $(CCFLAGS) $(LIBFLAGS)

%: %.cpp $(HEADERS)
	$(CCCOM) $< -o $@ $(CCFLAGS) $(LIBFLAGS)

.PHONY: clean
clean:
	find . -maxdepth 1 -type f -name '*.[0-9]*' ! -name '*.cpp' -exec rm {} +
