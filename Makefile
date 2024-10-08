CC			:= g++
CCFLAGS		:= -pedantic-errors -Wall -Wextra -Werror
CURR_DIR	:= $(shell pwd)
SRC_DIR		:= $(CURR_DIR)/src
OBJ_DIR		:= $(CURR_DIR)/obj
LIB_DIR		:= $(CURR_DIR)/lib
SRC			:= $(wildcard $(SRC_DIR)/*.cpp)
OBJ			:= $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o,$(SRC))
SFML		:= $(LIB_DIR)/SFML-2.5.1
LDFLAGS		:= -lsfml-graphics -lsfml-window -lsfml-system

TARGET 		:= main

main: $(OBJ)
	g++ -o $@ $^ $(LDFLAGS) 
	./$(TARGET)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	g++ $(CCFLAGS) -c -o $@ $<

recompile:
	make clean
	make main

clean:
	rm $(OBJ_DIR)/*
	rm $(TARGET)

