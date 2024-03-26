TARGET1=calc_dxy_BeforeAlign
TARGET2=fit_deltaXY
TARGET3=calc_dxy2
all: calc_dxy_BeforeAlign fit_deltaXY calc_dxy2
FEDRALIBS := -lEIO -lEdb -lEbase -lEdr -lScan -lAlignment -lEmath -lEphys -lvt -lDataConversion
$(TARGET1): $(TARGET1).cpp
	g++ $(TARGET1).cpp -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $(TARGET1)
$(TARGET2): $(TARGET2).cpp
	g++ $(TARGET2).cpp -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $(TARGET2)
$(TARGET3): $(TARGET3).cpp
	g++ $(TARGET3).cpp -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $(TARGET3)
