# $Id: Makefile,v 1.8 2004/12/04 00:30:34 bsd Exp $

TARGET = spi-mcp3002

SRCS = spi-mcp3002.c
OBJS = spi-mcp3002.o


OPTIMIZE = -O3
CFLAGS   = -Wall ${OPTIMIZE}
LIBS     = -lrpip	

CC     = gcc


all :   depend
	make ${TARGET}

.c.s :
	gcc -c ${CFLAGS} -S -o $@ $<

${TARGET} : ${OBJS}
	gcc -o ${TARGET} ${OBJS} ${LIBS} 


clean :
	rm -f *.o *~ .depend Makefile.bak ${TARGET} 


depend :
	@cp Makefile Makefile.bak
	@awk '/# .Id/,/^# DO NOT DELETE/' Makefile > Makefile.new
	@${CC} ${CFLAGS} -MM ${SRCS} >> Makefile.new
	@if ! diff Makefile Makefile.new > /dev/null 2>&1 ; then \
	  mv Makefile.new Makefile; \
	else \
	  rm Makefile.new; \
	  rm Makefile.bak; \
	fi
	

# DO NOT DELETE (Dependencies follow)
spi-mcp3002.o: spi-mcp3002.c
