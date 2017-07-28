CODE  = main.cpp
NAME  = trimal

CODE2 = readAl.cpp
NAME2 = readal

CODE3 = statAl.cpp
NAME3 = statal

CC    = g++
FLAGS = -Wall -O2
LIBM  = -lm
DUGB  =

all:			readal trimal statal

preview:		alignment.o statisticsGaps.o utils.o similarityMatrix.o statisticsConservation.o sequencesMatrix.o compareFiles.o

trimal:			preview
			$(CC) $(FLAGS) $(DUGB) -o $(NAME) $(CODE) $(LIBM) alignment.o statisticsGaps.o utils.o similarityMatrix.o statisticsConservation.o sequencesMatrix.o compareFiles.o

readal:			preview
			$(CC) $(FLAGS) $(DUGB) -o $(NAME2) $(CODE2) $(LIBM) alignment.o statisticsGaps.o utils.o similarityMatrix.o statisticsConservation.o sequencesMatrix.o compareFiles.o

statal:			preview
			$(CC) $(FLAGS) $(DUGB) -o $(NAME3) $(CODE3) $(LIBM) alignment.o statisticsGaps.o utils.o similarityMatrix.o statisticsConservation.o sequencesMatrix.o compareFiles.o

alignment.o:		alignment.cpp rwAlignment.cpp autAlignment.cpp alignment.h values.h defines.h
			$(CC) $(FLAGS) $(DUGB) -c alignment.cpp rwAlignment.cpp autAlignment.cpp

statisticsGaps.o:	statisticsGaps.cpp statisticsGaps.h
			$(CC) $(FLAGS) $(DUGB) -c statisticsGaps.cpp

sequencesMatrix.o:	sequencesMatrix.cpp sequencesMatrix.h
			$(CC) $(FLAGS) $(DUGB) -c sequencesMatrix.cpp

statisticsConservation.o: statisticsConservation.cpp statisticsConservation.h
			$(CC) $(FLAGS) $(DUGB) -c statisticsConservation.cpp

similarityMatrix.o:	similarityMatrix.cpp similarityMatrix.h values.h
			$(CC) $(FLAGS) $(DUGB) -c similarityMatrix.cpp

utils.o:		utils.cpp utils.h values.h  defines.h
			$(CC) $(FLAGS) $(DUGB) -c utils.cpp

compareFiles.o:		compareFiles.cpp compareFiles.h
			$(CC) $(FLAGS) $(DUGB) -c compareFiles.cpp

clean:
			-rm -f $(NAME) $(NAME2)  $(NAME3) *.o *~
