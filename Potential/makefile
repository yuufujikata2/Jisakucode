OBJS = rseq.o rsq1.o rsq2.o norm1g.o mkgfac.o fctp0.o fctp.o gintsr.o

rseq.so : $(OBJS)
	$(CC) -shared -o $@ $(OBJS) 

CC = gcc

clean :
	rm $(OBJS)

