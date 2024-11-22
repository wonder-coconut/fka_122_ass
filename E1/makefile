CC = gcc

LOCAL = ${HOME}/.local

CFLAGS = \
	-Wall \
	-pedantic \
	-Werror \
	-Iinclude \
	-g \
	-Wl,-rpath,${LOCAL}/lib \
	-L${LOCAL}/lib \
	-I${LOCAL}/include

CFLAGS_OPT = \
	-O2 \
	-march=native \
	-fsanitize=address \
	-fno-omit-frame-pointer

LIB = \
    -lm \
    -lgsl \
    -lgslcblas

OBJ = \
    obj/tools.o \
	obj/velocity_verlet.o \
	obj/fput.o \
	obj/transform.o

MAIN = \
	obj/main.o

CFLAGS += $(CFLAGS_OPT)

program: obj _program

_program: $(MAIN) $(OBJ)
	$(CC) $(CFLAGS) $^ -o program $(LIB)

obj/%.o: src/%.c
	$(CC) -c $(CFLAGS) $< -o $@

obj:
	mkdir -p obj

clean:
	find -iname "*.o" -exec rm {} \;
	find -iname "*.d" -exec rm {} \;
	rm -f program
	rm -rf obj

.PHONY: clean
