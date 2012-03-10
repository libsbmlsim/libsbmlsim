PROG = ../sim
BUILD_LIB = libsbmlsim.a
BUILD_TEST = test
OBJS_COMMON = set_local_para_as_value.o get_equation.o \
calc.o calc_event.o simulate_explicit.o prepare_algebraic.o \
lu_decomposition.o lu_solve.o simulate_implicit.o \
create_myResult_content.o create_mySBML_objects.o free_mySBML_objects.o \
output_result.o print_result_list.o \
calc_temp_value.o forwarding_value.o calc_k.o print_node_type.o \
check_num.o create_calc_object_list.o str_is_number.o chomp.o \
alter_tree_structure.o check_math.o check_AST.o math_functions.o \
alg_alter_tree_structure.o ev_alter_tree_structure.o \
assignment_alter_tree_structure.o myASTNode_func.o \
prepare_reversible_fast_reaction.o calc_initial_assignment.o \
initialize_delay_val.o substitute_delay_val.o util.o
OBJS = $(OBJS_COMMON) main.o
OBJS_LIB = $(OBJS_COMMON) lib_main.o
OBJS_TEST = test.o

CFLAGS = -Wall -O2 -g

CC = gcc

.PHONY: all
all: $(PROG) $(BUILD_LIB) $(BUILD_TEST)

### make object files ###
.SUFFIXES: .o .c
.c.o:
	$(CC) $(CFLAGS) -c $<

### make executable file ###
$(PROG): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(PROG) -L. -lsbml

### make static library ###
$(BUILD_LIB): $(OBJS_LIB)
	ar rcs $(BUILD_LIB) $(OBJS_LIB)

### make test binary with library ###
$(BUILD_TEST): $(OBJS_TEST) $(BUILD_LIB)
	$(CC) $(CFLAGS) $(OBJS_TEST) -o test -L. -lsbmlsim -lsbml

### include ###
$(OBJS): header.h

### clean up! ###
.PHONY: clean
clean:
	rm -f $(PROG) $(BUILD_LIB) $(BUILD_TEST) $(OBJS) $(OBJS_LIB) $(OBJS_TEST)
