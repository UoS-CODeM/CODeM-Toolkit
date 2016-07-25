# Target name
TARGET:=CODeM

# Specify your c++ compiler here
CXX=g++

# Specify source folders
SRCDIR:= . #core misc 

# ===============================
SOURCES := $(shell find $(SRCDIR) -name '*.cpp')
OBJS:= $(SOURCES:.cpp=.o)

CXXINC=-I$(shell pwd) -I$(shell pwd)/core
CXXFLAG=-g -O2 -std=c++11 $(CXXINC)

.PHONY: clean all

all: $(TARGET)

$(TARGET): $(OBJS) 
	@echo " Linking..."
	$(CXX) $^ -o $(TARGET)

$(OBJS): %.o : %.cpp
	$(CXX) $(CXXFLAG) -o $@ -c $<

test:
	./$(TARGET) --default
clean:
	@echo " Cleaning..."; 
	rm -f $(OBJS)
	rm -f $(TARGET)
