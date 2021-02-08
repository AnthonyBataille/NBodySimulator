CXX = g++

PROJDIR := $(realpath $(CURDIR)/..)
SRCDIR := $(PROJDIR)/src
BUILDDIR := $(PROJDIR)/build
HEADERSDIR :=$(PROJDIR)/headers

VERBOSE = TRUE

SDL_LIB := -lSDL2
SDL_INCLUDE := -I/usr/include/SDL2 -D_REENTRANT
INCLUDE := -I$(HEADERSDIR) -I$(SRCDIR)


CXXFLAGS := -Wall -Wextra -O3 -c -std=c++11 $(SDL_INCLUDE) $(INCLUDE) -fopenmp -DNDEBUG
LDFLAGS := $(SDL_LIB) -fopenmp

EXE = NBodySimulator
SRC = $(wildcard $(SRCDIR)/*.cpp)
OBJ = $(subst $(SRCDIR), $(BUILDDIR), $(SRC:.cpp=.o))

all: $(EXE)

$(EXE): $(OBJ)
	$(CXX) $^ $(LDFLAGS) -o $@

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $< -o $@


.PHONY:
	clean

clean:
	rm $(BUILDDIR)/*.o

mrproper:
	rm $(EXE)
