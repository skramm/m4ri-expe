COLOR_1=-e "\e[1;33m"
COLOR_2=-e "\e[1;34m"
COLOR_OFF="\e[0m"


SRC_DIR=.
BIN_DIR=build
OBJ_DIR=build

SRC_FILES=$(wildcard *.cpp)
HEADERS=$(wildcard *.h*)

EXEC_FILES=$(patsubst %.cpp,$(BIN_DIR)/%,$(SRC_FILES))
#EXEC_FILES=$(patsubst $(SRC_DIR)/%.cpp,$(BIN_DIR)/%,$(SRC_FILES))

LDFLAGS +=-lm4ri

all: $(EXEC_FILES)
	@echo "done target @<"


# generic compile rule
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	@echo $(COLOR_2) " - Compiling app file $<." $(COLOR_OFF)
	$(CXX) -o $@ -c $< $(CFLAGS)

# linking
# -s option: strip symbol (don't add if debugging)
$(BIN_DIR)/%: $(OBJ_DIR)/%.o
	@echo $(COLOR_1) " - Link demo $@." $(COLOR_OFF)
	$(CXX) -o $@ $<  $(LDFLAGS)

show: $(SRC_FILES)
	@echo SRC_FILES=$(SRC_FILES)
#	@echo OBJ_FILES=$(OBJ_FILES)
	@echo EXEC_FILES=$(EXEC_FILES)
