# LBM Makefile

TARGET= LBM
CXX= g++
CXXFLAGS=-Wall -O3 -fopenmp -march=native
SRCDIR= ./src
OBJDIR= ./obj
RESTART= ./restart
SOURCE= $(wildcard $(SRCDIR)/*.cpp)
OBJECT= $(addprefix $(OBJDIR)/, $(notdir $(SOURCE:.cpp=.o)))
OUTPUT= out.$(TARGET).txt field*.vtr output0* rst stt output0*.vtm geometryflag.vtm geometryflag stop

$(TARGET): $(OBJECT)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECT) -lm

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $< 

.PHONY: clean
clean: 
	rm -rf $(TARGET) $(OBJECT) $(OUTPUT) $(OBJDIR) $(RESTART)