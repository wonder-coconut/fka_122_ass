CC ?= gcc # c-compiler
CFLAGS = -Wl # compiler flags such as -O3 -W
LIBS = -lgsl -lgslcblas # libraries such as -lm -lgsl -lgslcblas

OBJECTS = C2_task2.o # object files to produce
PROGRAM = C2 # name of executable to be built

LOCAL = ${HOME}/.local
CFLAGS += \
	-Wl,-rpath,${LOCAL}/lib \
	-I${LOCAL}/include \
	-L${LOCAL}/lib \


# build .o files
%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)
# compile program
all: $(PROGRAM)
$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
