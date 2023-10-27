TARGET1=efficiency
TARGET2=check_deltaXYdis
TARGET3=read_dxy
TARGET4=plot_sigmaPar
TARGET5=test_FnuDivideAlign
TARGET6=test_FnuQualityCheck
# TARGET7=calc_dxy
TARGET8=align_and_measure_momentum

FEDRALIBS := -lEIO -lEdb -lEbase -lEdr -lScan -lAlignment -lEmath -lEphys -lvt -lDataConversion
CUDA_ROOT=/usr/local/cuda
MY_TOOL=/home/kokui/LEPP/FASERnu/Tools

all: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET8)

$(TARGET1): $(TARGET1).cpp
	g++ $^ -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $@

$(TARGET2): $(TARGET2).cpp
	g++ $^ -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $@

$(TARGET3): $(TARGET3).cpp
	g++ $^ -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $@

$(TARGET4): $(TARGET4).cpp
	g++ $^ -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $@

$(TARGET5): $(TARGET5).cpp src/FnuDivideAlign.cu
	nvcc $^ -Iinclude -I`root-config --incdir` -I$(FEDRA_ROOT)/include -I$(CUDA_ROOT)/include -I$(CUDA_ROOT)/samples/common/inc -L`root-config --libdir` -L$(FEDRA_ROOT)/lib -lCore -lEve -lMathCore -lRint -lThread -lTree -lRIO -lASImage -lGpad -lHist -lGraf -lGraf3d -lcudart -lPhysics -lEdb -lEIO -lEbase -lEdr -lvt -lEmath -lAlignment -lEphys -lDataConversion -o $@ -w

$(TARGET6): $(TARGET6).cpp src/FnuQualityCheck.cpp include/FnuQualityCheck.h
	g++ $^ -Iinclude -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $@

# can't compile with ROOT6
# $(TARGET7): $(TARGET7).cu
# 	nvcc $^ -I`root-config --incdir` -I$(FEDRA_ROOT)/include -I$(CUDA_ROOT)/include -I$(CUDA_ROOT)/samples/common/inc -L`root-config --libdir` -L$(FEDRA_ROOT)/lib -lCore -lEve -lMathCore -lRint -lThread -lTree -lRIO -lASImage -lGpad -lHist -lGraf -lGraf3d -lcudart -lPhysics -lEdb -lEIO -lEbase -lEdr -lvt -lEmath -lAlignment -lEphys -lDataConversion -o $@ -w

$(TARGET8) : align_and_measure_momentum.cpp FnuMomCoord.o FnuDivideAlign.o FnuQualityCheck.o
	nvcc -o $@ $^ -w -I$(MY_TOOL)/FnuMomCoord/include -I$(MY_TOOL)/DataQualityCheck/include -I`root-config --incdir` -I$(FEDRA_ROOT)/include -I$(CUDA_ROOT)/include -I$(CUDA_ROOT)/samples/common/inc -L`root-config --libdir` -L$(FEDRA_ROOT)/lib -lCore -lEve -lMathCore -lRint -lThread -lTree -lRIO -lASImage -lGpad -lHist -lGraf -lGraf3d -lcudart -lPhysics -lEdb -lEIO -lEbase -lEdr -lvt -lEmath -lAlignment -lEphys -lDataConversion

FnuMomCoord.o : $(MY_TOOL)/FnuMomCoord/src/FnuMomCoord.cpp
	g++ -c $< -w -I$(MY_TOOL)/FnuMomCoord/include `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` `root-config --glibs` `root-config --evelibs`

FnuQualityCheck.o : $(MY_TOOL)/DataQualityCheck/src/FnuQualityCheck.cpp
	g++ -c $< -w -I$(MY_TOOL)/DataQualityCheck/include `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` `root-config --glibs` `root-config --evelibs`

FnuDivideAlign.o : $(MY_TOOL)/DataQualityCheck/src/FnuDivideAlign.cu
	nvcc -c $< -w -I$(MY_TOOL)/DataQualityCheck/include -I`root-config --incdir` -I$(FEDRA_ROOT)/include -I$(CUDA_ROOT)/include -I$(CUDA_ROOT)/samples/common/inc -L`root-config --libdir` -L$(FEDRA_ROOT)/lib -lCore -lEve -lMathCore -lRint -lThread -lTree -lRIO -lASImage -lGpad -lHist -lGraf -lGraf3d -lcudart -lPhysics -lEdb -lEIO -lEbase -lEdr -lvt -lEmath -lAlignment -lEphys -lDataConversion

clean:
	$(RM) $(TARGET1)
	$(RM) $(TARGET2)
	$(RM) $(TARGET3)
	$(RM) $(TARGET4)
	$(RM) $(TARGET5)
	$(RM) $(TARGET6)
#	$(RM) $(TARGET7)
	$(RM) $(TARGET8)
