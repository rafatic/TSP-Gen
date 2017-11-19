CC = gcc -Wall -O0 -g
LIBS = -lm
OBJETS = graphGenetic.o graphTools.o genetique.o geneticAlgo.o
EXEC= genetic

all: $(EXEC)
genetic: $(OBJETS)
	$(CC) -o $@ $(OBJETS) $(LIBS)

%.o: %.c %.h
	$(CC) -c $<

clean:
	rm -rf *.o $(EXEC)
mrproper:
	rm -rf *.o $(EXEC)
