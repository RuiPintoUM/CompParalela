CPP = g++ -O3 -Wall -pg -mavx -ftree-vectorize -march=native -msse4 -Ofast 
SRCS = main.cpp fluid_solver.cpp EventManager.cpp

all:
    $(CPP) $(SRCS) -o fluid_sim

clean:
    @echo Cleaning up...
    @rm fluid_sim
	@echo Done.