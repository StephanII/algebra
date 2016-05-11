
CC = llvm-g++
CPP_FILES := $(wildcard algebra/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))
LD_FLAGS := 
CC_FLAGS := -Wall -pedantic -O3 #(optimierung) # -g (for debugging)

test.exe: test.o $(OBJ_FILES)
	$(CC) $(LD_FLAGS) -o $@ $^

obj/%.o: algebra/%.cpp algebra/%.hpp
	$(CC) $(CC_FLAGS) -c -o $@ $<
	
test.o: test/test.cpp
	$(CC) $(CC_FLAGS) -c -o $@ $<

all: test.exe
	@echo "done"

clean:
	@rm -f obj/*.o *.exe *.o *.stackdump
	@echo project cleaned

