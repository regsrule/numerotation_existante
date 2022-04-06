CXXFLAGS += -Wall -Wextra -Wcast-align -Wcast-qual -Wconversion -Wfloat-equal \
	    -Wformat=2 -Winit-self -Wmissing-declarations \
	    -Wmissing-include-dirs -Wpointer-arith -Wredundant-decls \
	    -Wswitch-default -Wuninitialized -Wwrite-strings \
	    -Wno-sign-conversion -Wno-unused-function \
            -Wno-missing-declarations \
            -fopenmp -std=c++14 -mcx16 -O3 -DNDEBUG 

CC = g++
INC = ./
CFLAGS += -I$(INC) -g -O0 -lm -lpthread -fopenmp -ltcmalloc_minimal -lnuma
HEADERS = \
	Util.h\
	UnitHeap.h\
	Graph_gorder.h\

OBJ = ./obj

SRC =\
	main.cpp \
	Util.cpp\
	UnitHeap.cpp\
	Graph_gorder.cpp\

SRC_DIRS = ./

vpath %.h $(SRC_DIRS)
vpath %.cpp $(SRC_DIRS)

OBJS = $(addprefix $(OBJ)/, $(SRC:.cpp=.o))
	HEADS = $(addprefix $(INC)/, $(HEADERS))


gorder: $(SRC) $(OBJS) $(HEADERS)
	        $(CC) $(CFLAGS) $(OBJS) -o $@

$(OBJ)/%.o : %.cpp
	        $(CC) $(CFLAGS) -c $< -o $@

clean:
	        rm -f $(OBJ)/*.o gorder
