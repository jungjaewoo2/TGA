# You can put your build options here
all: fidlib.a

fidlib.a: fidlib.o
	$(AR) rc $@ $^

%.o: %.c fidlib.h  fidmkf.h  fidrf_cmdlist.h  fidrf_combined.h  fidrf_jit.h
	$(CC) -DT_LINUX -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o

.PHONY: all clean test

