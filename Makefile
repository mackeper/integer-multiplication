# $@ : name on the left side of a : (the rule)
# $< : first name of the right side of a :
# $^ : the name of the prerequisite that caused the rule to execute.

APP_NAME = imnln
DEBUG_TARGET = debug
RELEASE_TARGET = release
TEST_TARGET = tests
INT_TEST_TARGET = int_tests
MAIN = main
MAIN_TEST = $(TEST_TARGET)
MAIN_INT_TEST = $(INT_TEST_TARGET)
APP_NAME_TEST = $(TEST_TARGET)
APP_NAME_INT_TEST = $(INT_TEST_TARGET)

INC_DIR = include
OBJ_DIR = obj
OBJ_DIR_TEST = $(OBJ_DIR)/$(TEST_TARGET)
OBJ_DIR_INT_TEST = $(OBJ_DIR)/$(INT_TEST_TARGET)
BIN_DIR = bin
SRC_DIR = src
SRC_DIR_TEST = $(SRC_DIR)/$(TEST_TARGET)
SRC_DIR_INT_TEST = $(SRC_DIR)/$(INT_TEST_TARGET)

CXX = g++
LIBS = -lgmpxx -lgmp
CXXFLAGS = -I $(INC_DIR) -std=c++17 -Wfatal-errors
DEBUG_FLAGS = -g -fsanitize=undefined -Wall -Wconversion -pedantic -DDEBUG=true
RELEASE_FLAGS = -O3 -DDEBUG=false

SOURCES = $(shell find $(SRC_DIR) -name *.cpp ! -name $(MAIN).cpp ! -path "$(SRC_DIR_TEST)/*" ! -path "$(SRC_DIR_INT_TEST)/*")
SOURCES_TEST = $(shell find $(SRC_DIR_TEST) -name *.cpp ! -name $(MAIN_TEST).cpp )
SOURCES_INT_TEST = $(shell find $(SRC_DIR_INT_TEST) -name *.cpp ! -name $(MAIN_INT_TEST).cpp )
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))
OBJECTS_TEST = $(patsubst $(SRC_DIR_TEST)/%.cpp,$(OBJ_DIR_TEST)/%.o,$(SOURCES_TEST))
OBJECTS_INT_TEST = $(patsubst $(SRC_DIR_INT_TEST)/%.cpp,$(OBJ_DIR_INT_TEST)/%.o,$(SOURCES_INT_TEST))

all: CXXFLAGS += $(DEBUG_FLAGS)
all: $(DEBUG_TARGET) $(TEST_TARGET) $(INT_TEST_TARGET)

$(DEBUG_TARGET): CXXFLAGS += $(DEBUG_FLAGS)
$(DEBUG_TARGET): $(APP_NAME)

$(RELEASE_TARGET): CXXFLAGS += $(RELEASE_FLAGS)
$(RELEASE_TARGET): clean $(APP_NAME)

$(TEST_TARGET): CXXFLAGS += $(DEBUG_FLAGS)
$(TEST_TARGET): $(BIN_DIR)/$(APP_NAME_TEST)

$(INT_TEST_TARGET): CXXFLAGS += $(DEBUG_FLAGS)
$(INT_TEST_TARGET): $(BIN_DIR)/$(APP_NAME_INT_TEST)

$(APP_NAME): $(BIN_DIR)/$(APP_NAME)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp 
	@mkdir -p $(@D)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(OBJ_DIR_TEST)/%.o: $(SRC_DIR_TEST)/%.cpp
	@mkdir -p $(@D)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(OBJ_DIR_INT_TEST)/%.o: $(SRC_DIR_INT_TEST)/%.cpp
	@mkdir -p $(@D)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(BIN_DIR)/$(APP_NAME): $(OBJECTS) $(OBJ_DIR)/$(MAIN).o
	@mkdir -p $(@D)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

$(BIN_DIR)/$(APP_NAME_TEST): $(OBJECTS_TEST) $(OBJECTS) $(OBJ_DIR_TEST)/$(MAIN_TEST).o
	@mkdir -p $(@D)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

$(BIN_DIR)/$(APP_NAME_INT_TEST): $(OBJECTS_INT_TEST) $(OBJECTS) $(OBJ_DIR_INT_TEST)/$(MAIN_INT_TEST).o
	@mkdir -p $(@D)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

.PHONY: clean run runtests runinttests runalltests

run: $(BIN_DIR)/$(APP_NAME)
	./$(BIN_DIR)/$(APP_NAME)

runtests: $(BIN_DIR)/$(APP_NAME_TEST)
	./$(BIN_DIR)/$(APP_NAME_TEST)

runinttests: $(BIN_DIR)/$(APP_NAME_INT_TEST)
	./$(BIN_DIR)/$(APP_NAME_INT_TEST)

runalltests: $(BIN_DIR)/$(APP_NAME_INT_TEST) $(BIN_DIR)/$(APP_NAME_TEST)
	./$(BIN_DIR)/$(APP_NAME_TEST)
	./$(BIN_DIR)/$(APP_NAME_INT_TEST)

clean:
	rm -rf $(OBJ_DIR_TEST) $(OBJ_DIR) $(BIN_DIR)