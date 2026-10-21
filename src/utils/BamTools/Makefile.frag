src/utils/BamTools/include/BamAlignment.mapping.hpp: src/utils/BamTools/mapping/BamAlignment.py src/utils/BamTools/mapping/BamAlignment.map

src/utils/BamTools/include/%.mapping.hpp: src/utils/BamTools/mapping/%.py src/utils/BamTools/mapping/%.map
	$(PYTHON) $^ > $@

# bedtools-crossplatform: stream-callback hFILE backend (was in the bundled htslib).
BUILT_OBJECTS += obj/HtsCallback.o

obj/HtsCallback.o: src/utils/BamTools/HtsCallback.cpp obj/HtsCallback.d
	$(CXX_COMPILE)
