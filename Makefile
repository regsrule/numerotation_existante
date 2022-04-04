CXXFLAGS += -Wall -Wextra -Wcast-align -Wcast-qual -Wconversion -Wfloat-equal \
	    -Wformat=2 -Winit-self -Wmissing-declarations \
	    -Wmissing-include-dirs -Wpointer-arith -Wredundant-decls \
	    -Wswitch-default -Wuninitialized -Wwrite-strings \
	    -Wno-sign-conversion -Wno-unused-function \
            -Wno-missing-declarations \
            -fopenmp -std=c++14 -mcx16 -O3 -DNDEBUG 

CC = g++
INC = ./
CFLAGS += -I$(INC) -O0 -lm -lpthread -fopenmp -ltcmalloc_minimal -lnuma -DGCC -DRelease  -m64 -g 

HEADERS = \
	function.h \
	graph_cn.h \
	time_manager.h\
	community.h\
	graph_binary.h\
	graph.h\
	edge_list.h\
	rabbit_order.h\
	reorder.h\
	Util.h\
	UnitHeap.h\
	Graph_gorder.h\

OBJ = ./obj

SRC =\
	function.cpp\
	graph_cn.cpp\
	time_manager.cpp\
	main.cpp \
	community.cpp\
	graph_binary.cpp\
	graph.cpp\
	reorder.cpp\
	rabbit_order.cpp\
	edge_list.cpp\
	Util.cpp\
	UnitHeap.cpp\
	Graph_gorder.cpp\

SRC_DIRS = ./ louvain rabbit gorder

vpath %.h $(SRC_DIRS)
vpath %.cpp $(SRC_DIRS)

OBJS = $(addprefix $(OBJ)/, $(SRC:.cpp=.o))
	HEADS = $(addprefix $(INC)/, $(HEADERS))



cn-order: $(SRC) $(OBJS) $(HEADERS)
	        $(CC) $(OBJS) -o  $@ $(CFLAGS)

$(OBJ)/%.o : %.cpp
	        $(CC) $(CFLAGS) -c $< -o $@

clean:
	        rm -f $(OBJ)/*.o cn-order


run_valgrind:cn-order
	valgrind  --tool=memcheck --leak-check=yes --leak-resolution=low --show-reachable=yes ./cn-order louvain/sample_networks/archiv.txt 10 10

run_gdb:cn-order
	gdb  ./cn-order louvain/sample_networks/archiv.txt 10 10
