
CC = llvm-g++
CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))
LD_FLAGS := #-lOpenCL  
CC_FLAGS := -Wall -std=c++11 -pedantic -O3 #(optimierung) # -g (for debugging)

runtest.exe: obj/test.o $(OBJ_FILES)
	$(CC) $(LD_FLAGS) -o $@ $^

all: runtest.exe
	@echo "done"

obj/%.o: src/%.cpp src/%.hpp
	$(CC) $(CC_FLAGS) -c -o $@ $<
	
obj/test.o: test/test.cpp
	$(CC) $(CC_FLAGS) -c -o $@ $<

.PHONY: objects

objects: $(OBJ_FILES)
	@echo object files created	

clean:
	@rm -f obj/*.o *.exe *.o *.stackdump
	@echo project cleaned

