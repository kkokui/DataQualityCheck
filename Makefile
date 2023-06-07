TARGET1=efficiency
TARGET2=check_deltaXYdis
TARGET3=read_dxy
TARGET4=plot_sigmaPar

FEDRALIBS := -lEIO -lEdb -lEbase -lEdr -lScan -lAlignment -lEmath -lEphys -lvt -lDataConversion
CUDA_ROOT=/usr/local/cuda

all: calc_dxy efficiency check_deltaXYdis read_dxy plot_sigmaPar

calc_dxy: calc_dxy.cu
	nvcc calc_dxy.cu  -I`root-config --incdir` -I$(FEDRA_ROOT)/include -I$(CUDA_ROOT)/include -I$(CUDA_ROOT)/samples/common/inc -L`root-config --libdir` -L$(FEDRA_ROOT)/lib -lCore -lEve -lMathCore -lRint -lThread -lTree -lRIO -lASImage -lGpad -lHist -lGraf -lGraf3d -lcudart -lPhysics -lEdb -lEIO -lEbase -lEdr -lvt -lEmath -lAlignment -lEphys -lDataConversion -o calc_dxy

$(TARGET1): $(TARGET1).cpp
	g++ $(TARGET1).cpp -w -I src `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $(TARGET1)

$(TARGET2): $(TARGET2).cpp
	g++ $(TARGET2).cpp -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $(TARGET2)

$(TARGET3): $(TARGET3).cpp
	g++ $(TARGET3).cpp -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $(TARGET3)

$(TARGET4): $(TARGET4).cpp
	g++ $(TARGET4).cpp -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $(TARGET4)

clean:
	$(RM) $(TARGET)
