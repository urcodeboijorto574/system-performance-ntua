CFLAGS = -Wall -O2

# Detect the OS
ifeq ($(OS), Windows_NT)
    EXE=pe241_code.exe
    RM=rm
    RUN=./pe241_code.exe
else
    EXE=pe241_code
    RM=rm -f
    RUN=./pe241_code
endif

$(EXE): pe241_code.cpp
	g++ -o $(EXE) pe241_code.cpp $(CFLAGS)

run1: $(EXE)
	$(RUN)

out1: $(EXE)
	$(RUN) > output.txt

clean:
	$(RM) $(EXE) output.txt
