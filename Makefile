# Generated by Kimi
# 定义编译器
CXX = clang++

# 定义编译选项，可以根据需要修改
CXXFLAGS = -Wall -O0 -std=c++17 -m64 -mavx2 -fopenmp -D BREAKDOWN_TIME -D STRONG_CHECKSUM

# 定义源文件目录
SRC_DIR = source

# 定义目标文件目录
BIN_DIR = bin

# 自动生成源文件列表
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BIN_DIR)/%,$(SRCS))

# 默认目标
all: $(OBJS)
	@echo "All programs compiled successfully."

# 编译单个源文件
$(BIN_DIR)/%: $(SRC_DIR)/%.cpp
	@mkdir -p $(BIN_DIR) # 确保目标目录存在
	@echo "Compiling $<..."
	@$(CXX) $(CXXFLAGS) -o $@ $<

# 清理所有编译生成的文件
clean:
	@rm -rf $(BIN_DIR)/*

# 运行所有程序
run: all
	@for bin in $(BIN_DIR)/*; do \
		echo "------------ Running $$bin ------------"; \
		"$$bin" 233; \
	done

# 防止与文件名冲突
.PHONY: all clean run