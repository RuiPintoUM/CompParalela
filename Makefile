CPP = g++ -Wall -Ofast -funroll-loops -march=native -flto -mavx -ffast-math -fno-strict-aliasing 
SRCS = main.cpp fluid_solver.cpp EventManager.cpp

all:
	$(CPP) $(SRCS) -o fluid_sim

clean:
	@echo Cleaning up...
	@rm fluid_sim
	@echo Done.
