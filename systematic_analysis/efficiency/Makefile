TARGET1=efficiency_plate_pos

all: $(TARGET1)
FEDRALIBS := -lEIO -lEdb -lEbase -lEdr -lScan -lAlignment -lEmath -lEphys -lvt -lDataConversion

$(TARGET1): $(TARGET1).cpp
	g++ $(TARGET1).cpp -w `root-config --cflags` -I$(FEDRA_ROOT)/include -L$(FEDRA_ROOT)/lib  $(FEDRALIBS) `root-config --libs` -o $(TARGET1)
