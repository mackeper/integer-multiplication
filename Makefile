# $@ : name on the left side of a : (the rule)
# $< : first name of the right side of a :
# $^ : the name of the prerequisite that caused the rule to execute.

APP_NAME = main
DEBUG_TARGET = debug
RELEASE_TARGET = release

INC_DIR=include
CXX=g++
CXXFLAGS = -I $(INC_DIR) -std=c++17 -Wfatal-errors

OBJ_DIR=obj
BIN_DIR=bin
SRC_DIR=src
LIBS=-lgmpxx -lgmp

SOURCES = $(shell find $(SRC_DIR) -name *.cpp)
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))

all: $(DEBUG_TARGET)

$(DEBUG_TARGET): CXXFLAGS += -g -fsanitize=undefined -Wall -Wconversion -pedantic
$(DEBUG_TARGET): clean $(BIN_DIR)/$(APP_NAME)

$(RELEASE_TARGET): CXXFLAGS += -O3
$(RELEASE_TARGET): clean $(BIN_DIR)/$(APP_NAME)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(BIN_DIR)/$(APP_NAME): $(OBJECTS)
	@echo $(SOURCES)
	@mkdir -p $(@D)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)