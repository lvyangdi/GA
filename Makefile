OBJS = graph.o traverse.o genetic.o
HEADER := $(patsubst %.o,%.h,$(OBJS))

CC = g++
CFLAGS = -O3 -std=c++11 -Wall

all: ga

ga: main.cpp $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^
clean:
	rm -f ga *.o

%.o: %.cpp $(HEADER)
	$(CC) $(CFLAGS) -c -o $@ $<
