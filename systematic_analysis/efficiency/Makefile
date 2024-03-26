TARGET1=efficiency
TARGET2=plot_efficiency
TARGET3=efficiency_plate_pos

all: efficiency plot_efficiency $(TARGET3)
FEDRALIBS := -lEIO -lEdb -lEbase -lEdr -lScan -lAlignment -lEmath -lEphys -lvt -lDataConversion
$(TARGET1): $(TARGET1).cpp
	g++ $(TARGET1).cpp -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $(TARGET1)
$(TARGET2): $(TARGET2).cpp
	g++ $(TARGET2).cpp -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $(TARGET2)
$(TARGET3): $(TARGET3).cpp
	g++ $(TARGET3).cpp -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $(TARGET3)
