# $@ : name on the left side of a : (the rule)
# $< : first name of the right side of a :
# $^ : the name of the prerequisite that caused the rule to execute.

APP_NAME = imnln
DEBUG_TARGET = debug
RELEASE_TARGET = release
TEST_TARGET = tests
MAIN = main
MAIN_TEST = $(TEST_TARGET)
APP_NAME_TEST = $(TEST_TARGET)

INC_DIR = include
OBJ_DIR = obj
OBJ_DIR_TEST = $(OBJ_DIR)/$(TEST_TARGET)
BIN_DIR = bin
SRC_DIR = src
SRC_DIR_TEST = $(SRC_DIR)/$(TEST_TARGET)

CXX = g++
LIBS = -lgmpxx -lgmp
CXXFLAGS = -I $(INC_DIR) -std=c++17 -Wfatal-errors
DEBUG_FLAGS = -g -fsanitize=undefined -Wall -Wconversion -pedantic
RELEASE_FLAGS = -O3

SOURCES = $(shell find $(SRC_DIR) -name *.cpp ! -name $(MAIN).cpp ! -path "$(SRC_DIR_TEST)/*")
SOURCES_TEST = $(shell find $(SRC_DIR_TEST) -name *.cpp ! -name $(MAIN_TEST).cpp )
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))
OBJECTS_TEST = $(patsubst $(SRC_DIR_TEST)/%.cpp,$(OBJ_DIR_TEST)/%.o,$(SOURCES_TEST))

all: $(DEBUG_TARGET)

$(DEBUG_TARGET): CXXFLAGS += $(DEBUG_FLAGS)
$(DEBUG_TARGET): $(APP_NAME) $(TEST_TARGET)

$(RELEASE_TARGET): CXXFLAGS += $(RELEASE_FLAGS)
$(RELEASE_TARGET): clean $(APP_NAME)

$(TEST_TARGET): CXXFLAGS += $(DEBUG_FLAGS)
$(TEST_TARGET): $(BIN_DIR)/$(APP_NAME_TEST)

$(APP_NAME): $(BIN_DIR)/$(APP_NAME)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp 
	@mkdir -p $(@D)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(OBJ_DIR_TEST)/%.o: $(SRC_DIR_TEST)/%.cpp
	@mkdir -p $(@D)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(BIN_DIR)/$(APP_NAME): $(OBJECTS) $(OBJ_DIR)/$(MAIN).o
	@mkdir -p $(@D)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

$(BIN_DIR)/$(APP_NAME_TEST): $(OBJECTS_TEST) $(OBJECTS) $(OBJ_DIR_TEST)/$(MAIN_TEST).o
	@mkdir -p $(@D)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

.PHONY: clean run runtests

run: $(BIN_DIR)/$(APP_NAME)
	./$(BIN_DIR)/$(APP_NAME)

runtests: $(BIN_DIR)/$(APP_NAME_TEST)
	./$(BIN_DIR)/$(APP_NAME_TEST)

clean:
	rm -rf $(OBJ_DIR_TEST) $(OBJ_DIR) $(BIN_DIR)