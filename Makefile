# This Makefile is meant to be used on Windows OS. Might modify it so it can be ran in multiple OSs.

CFLAGS = -Wall -O2

pe241_code.exe: pe241_code.cpp
	g++ -o pe241_code.exe pe241_code.cpp

run1: pe241_code.exe
	.\pe241_code.exe

out1: pe241_code.exe
	./pe241_code.exe > output.txt

clean:
	rm pe241_code.exe
