TARGET1=divide_align
TARGET2=quality_check
TARGET3=measure_momentum

FEDRALIBS := -lEIO -lEdb -lEbase -lEdr -lScan -lAlignment -lEmath -lEphys -lvt -lDataConversion
CUDA_ROOT=/usr/local/cuda
MY_TOOL=/home/kokui/LEPP/FASERnu/Tools

all: $(TARGET1) $(TARGET2) $(TARGET3)

$(TARGET1): $(TARGET1).cpp FnuDivideAlign.o
	nvcc $^ -Iinclude -I`root-config --incdir` -I$(FEDRA_ROOT)/include -I$(CUDA_ROOT)/include -I$(CUDA_ROOT)/samples/common/inc -L`root-config --libdir` -L$(FEDRA_ROOT)/lib -lCore -lEve -lMathCore -lRint -lThread -lTree -lRIO -lASImage -lGpad -lHist -lGraf -lGraf3d -lcudart -lPhysics -lEdb -lEIO -lEbase -lEdr -lvt -lEmath -lAlignment -lEphys -lDataConversion -o $@ -w

$(TARGET2): $(TARGET2).cpp FnuQualityCheck.o FnuResolution.o FnuEfficiency.o FnuPositionDistribution.o FnuAngleDistribution.o
	g++ $^ -Iinclude -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $@

$(TARGET3) : $(TARGET3).cpp FnuMomCoord.o
	g++ $^ -I$(MY_TOOL)/FnuMomCoord/include -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $@


OBJECT1=FnuMomCoord.o
OBJECT2=FnuQualityCheck.o
OBJECT3=FnuDivideAlign.o
OBJECT4=FnuResolution.o
OBJECT5=FnuEfficiency.o
OBJECT6=FnuPositionDistribution.o
OBJECT7=FnuAngleDistribution.o

$(OBJECT1) : $(MY_TOOL)/FnuMomCoord/src/FnuMomCoord.cpp
	g++ -c $< -w -I$(MY_TOOL)/FnuMomCoord/include `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` `root-config --glibs` `root-config --evelibs`

$(OBJECT2) : src/FnuQualityCheck.cpp
	g++ -c $< -w -Iinclude `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` `root-config --glibs` `root-config --evelibs`

$(OBJECT3) : src/FnuDivideAlign.cu
	nvcc -c $< -w -Iinclude -I`root-config --incdir` -I$(FEDRA_ROOT)/include -I$(CUDA_ROOT)/include -I$(CUDA_ROOT)/samples/common/inc -L`root-config --libdir` -L$(FEDRA_ROOT)/lib -lCore -lEve -lMathCore -lRint -lThread -lTree -lRIO -lASImage -lGpad -lHist -lGraf -lGraf3d -lcudart -lPhysics -lEdb -lEIO -lEbase -lEdr -lvt -lEmath -lAlignment -lEphys -lDataConversion

$(OBJECT4) : src/FnuResolution.cpp
	g++ -c $< -w -Iinclude `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` `root-config --glibs` `root-config --evelibs`

$(OBJECT5) : src/FnuEfficiency.cpp
	g++ -c $< -w -Iinclude `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` `root-config --glibs` `root-config --evelibs`

$(OBJECT6) : src/FnuPositionDistribution.cpp
	g++ -c $< -w -Iinclude `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` `root-config --glibs` `root-config --evelibs`

$(OBJECT7) : src/FnuAngleDistribution.cpp
	g++ -c $< -w -Iinclude `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` `root-config --glibs` `root-config --evelibs`

clean:
	$(RM) $(TARGET1)
	$(RM) $(TARGET2)
	$(RM) $(TARGET3)
	$(RM) $(OBJECT1)
	$(RM) $(OBJECT2)
	$(RM) $(OBJECT3)
	$(RM) $(OBJECT4)
	$(RM) $(OBJECT5)
	$(RM) $(OBJECT6)
	$(RM) $(OBJECT7)
