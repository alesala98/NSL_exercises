CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o class.o
	$(CC) random.o class.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
class.o : class.cpp class.h 
	$(CC) -c class.cpp -o class.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
remove:
	rm Datas/*.dat
