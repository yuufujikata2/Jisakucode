OBJ_RSEQ = rseq.c rsq1.c rsq2.c norm1g.c mkgfac.c fctp0.c fctp.c gintsr.c
DIR_RSEQ = ./Rseq
OBJS = $(patsubst %.c,$(DIR_RSEQ)/%.o,$(OBJ_RSEQ) )

CC = gcc

rseq.so : $(OBJS)
	$(CC) -shared -o $@ $(OBJS)


clean :
	rm $(OBJS)

