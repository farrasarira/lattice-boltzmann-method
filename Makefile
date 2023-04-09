# LBM Makefile

TARGET= LBM
CXX= g++
CXXFLAGS=-Wall -O3 -fopenmp -march=native
SRCDIR= ./src
OBJDIR= ./obj
RESTART= ./restart
SOURCE= $(wildcard $(SRCDIR)/*.cpp)
OBJECT= $(addprefix $(OBJDIR)/, $(notdir $(SOURCE:.cpp=.o)))
OUTPUT= out.$(TARGET).txt field*.vtr output0* rst stt output0*.vtm geometryflag.vtm geometryflag stop *.txt

CANTERA_INSTALL_ROOT=/usr
CANTERA_CORE_INCLUDES=-I$(CANTERA_INSTALL_ROOT)/include
CANTERA_EXTRA_INCLUDES=-I/usr/include/eigen3 
CANTERA_BOOST_INCLUDES=
CANTERA_SUNDIALS_INCLUDE=

CANTERA_INCLUDES=$(CANTERA_CORE_INCLUDES) $(CANTERA_SUNDIALS_INCLUDE) \
                 $(CANTERA_BOOST_INCLUDES) $(CANTERA_EXTRA_INCLUDES)

EIGEN_INCLUDES=-I./src/headers/Eigen

CANTERA_CORE_LIBS=-pthread -L/usr/lib -lcantera -lsundials_cvodes -lsundials_ida -lsundials_nvecserial -lsundials_sunlinsollapackdense -lsundials_sunlinsollapackband -lfmt -lyaml-cpp -llapack -lblas -Wl,-rpath,/usr/lib
CANTERA_EXTRA_LIBDIRS=
CANTERA_SUNDIALS_LIBS= -lsundials_cvodes -lsundials_ida -lsundials_nvecserial -lsundials_sunlinsollapackdense -lsundials_sunlinsollapackband
CANTERA_BOOST_INCLUDES=
CANTERA_BLAS_LAPACK_LIBS=-llapack -lblas

CANTERA_LIBS=$(CANTERA_CORE_LIBS) \
             $(CANTERA_EXTRA_LIBDIRS) $(CANTERA_SUNDIALS_LIBS) \
             $(CANTERA_BLAS_LAPACK_LIBS)

$(TARGET): $(OBJECT)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECT) $(CANTERA_LIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp 
	-mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $< -I./src/headers $(CANTERA_INCLUDES) $(EIGEN_INCLUDES)

.PHONY: clean
clean: 
	rm -rf $(TARGET) $(OBJECT) $(OUTPUT) $(OBJDIR) $(RESTART)