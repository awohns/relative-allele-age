CC := g++

CFLAGS := -std=c++11 -Wall -O3 -march=native
CCLIBS  := -lz -lpthread

TARGET := rvage_dev3

SUBDIRS := $(shell find ./src/ -type d)
INCDIRS := $(addprefix -I, $(SUBDIRS))

SOURCES := $(shell find . -name '*.cpp')
OBJECTS := $(patsubst %.cpp, %.o, $(SOURCES))


all: $(TARGET)

%.o: %.cpp
	$(CC) $(CFLAGS) -flto $(INCDIRS) -c $< -o $@ $(CCLIBS)

$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGET) -flto $(CFLAGS) $(OBJECTS) $(CCLIBS)


.PHONY: clean

clean:
	@rm -f $(OBJECTS) $(TARGET)


