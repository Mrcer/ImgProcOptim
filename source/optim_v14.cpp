/* Unimplemented Warning: 仅初始化 im2col 就已经难敌常规优化，因此不再继续往下写 
 * + 实现 im2col 优化访存连续性
 */
#include <chrono>
#include <cstddef>
#include <iostream>
#include <random>

#define ELEMENT(buffer, i, j, dim1) (buffer[(i) * (dim1) + (j)])

void print(unsigned char *buffer, size_t dim0, size_t dim1)
{
    for (size_t i = 0; i < dim0; ++i)
    {
        for (size_t j = 0; j < dim1; ++j)
        {
            std::cout << static_cast<int>(ELEMENT(buffer, i, j, dim1)) << " ";
        }
    }
}

class FigureProcessor
{
private:
    unsigned char *figure;
    unsigned char *result;
    unsigned char *im2col_buffer;
    const size_t size;
    const size_t buffer_size;

public:
    FigureProcessor(size_t size, size_t seed = 0) : size(size), buffer_size(size * size)
    {

        // !!! Please do not modify the following code !!!
        // 如果你需要修改内存的数据结构，请不要修改初始化的顺序和逻辑
        // 助教可能会通过指定某个初始化seed 的方式来验证你的代码
        // 如果你修改了初始化的顺序，可能会导致你的代码无法通过测试
        std::random_device rd;
        std::mt19937_64 gen(seed == 0 ? rd() : seed);
        std::uniform_int_distribution<unsigned char> distribution(0, 255);
        // !!! ----------------------------------------- !!!

        // 两个数组的初始化在这里，可以改动，但请注意 gen 的顺序是从上到下从左到右即可。
        figure = static_cast<unsigned char *>(malloc(buffer_size));
        result = static_cast<unsigned char *>(malloc(buffer_size));
        im2col_buffer = static_cast<unsigned char *>(malloc(9 * buffer_size));
        for (size_t i = 0; i < buffer_size; ++i)
        {
            figure[i] = static_cast<unsigned char>(distribution(gen));
            result[i] = 0;
        }
    }

    ~FigureProcessor()
    {
        free(figure);
        free(result);
    }

    void im2colInit()
    {
        unsigned char *src1;
        unsigned char *src2 = figure;
        unsigned char *src3 = figure + size;
        unsigned char *dst = im2col_buffer;
        unsigned char *dst_end;
#define APPEND(val) *dst = (val), dst += 1;

        APPEND(0)
        APPEND(0)
        APPEND(0)
        APPEND(0)
        APPEND(*src2)
        APPEND(*src3)
        APPEND(0)
        APPEND(*(src2 + 1))
        APPEND(*(src3 + 1))
        src2 += 1, src3 += 1;
        dst_end = dst + 9 * (size - 2);
        while(dst != dst_end)
        {
            APPEND(*(src2 - 1))
            APPEND(*(src2 - 1))
            APPEND(*(src3 - 1))
            APPEND(*src2)
            APPEND(*src2)
            APPEND(*src3)
            APPEND(*(src2 + 1))
            APPEND(*(src2 + 1))
            APPEND(*(src3 + 1))
            src2 += 1, src3 += 1;
        }
        APPEND(0)
        APPEND(*(src2 - 1))
        APPEND(*(src3 - 1))
        APPEND(0)
        APPEND(*src2)
        APPEND(*src3)
        APPEND(0)
        APPEND(0)
        APPEND(0)
        src2 += 1, src3 += 1;

        src1 = figure;
        for(int i = 1 ; i < size - 1 ; ++i)
        {
            APPEND(*src1)
            APPEND(*src2)
            APPEND(*src3)
            APPEND(*src1)
            APPEND(*src2)
            APPEND(*src3)
            APPEND(*(src1 + 1))
            APPEND(*(src2 + 1))
            APPEND(*(src3 + 1))
            src1 += 1, src2 += 1, src3 += 1;

            dst_end = dst + 9 * (size - 2);
            while (dst != dst_end)
            {
                APPEND(*(src1 - 1))
                APPEND(*(src2 - 1))
                APPEND(*(src3 - 1))
                APPEND(*src1)
                APPEND(*src2)
                APPEND(*src3)
                APPEND(*(src1 + 1))
                APPEND(*(src2 + 1))
                APPEND(*(src3 + 1))
                src1 += 1, src2 += 1, src3 += 1;
            }

            APPEND(*(src1 - 1))
            APPEND(*(src2 - 1))
            APPEND(*(src3 - 1))
            APPEND(*src1)
            APPEND(*src2)
            APPEND(*src3)
            APPEND(*src1)
            APPEND(*src2)
            APPEND(*src3)
            src1 += 1, src2 += 1, src3 += 1;
        }

        APPEND(0)
        APPEND(0)
        APPEND(0)
        APPEND(*src1)
        APPEND(*src2)
        APPEND(0)
        APPEND(*(src1 + 1))
        APPEND(*(src2 + 1))
        APPEND(0)
        src1 += 1, src2 += 1;
        dst_end = dst + 9 * (size - 2);
        while (dst != dst_end)
        {
            APPEND(*(src1 - 1))
            APPEND(*(src2 - 1))
            APPEND(*(src2 - 1))
            APPEND(*src1)
            APPEND(*src2)
            APPEND(*src2)
            APPEND(*(src1 + 1))
            APPEND(*(src2 + 1))
            APPEND(*(src2 + 1))
            src1 += 1, src2 += 1;
        }
        APPEND(*(src1 - 1))
        APPEND(*(src2 - 1))
        APPEND(0)
        APPEND(*src1)
        APPEND(*src2)
        APPEND(0)
        APPEND(0)
        APPEND(0)
        APPEND(0)
#undef APPEND
    }

    // Gaussian filter
    // [[1, 2, 1], [2, 4, 2], [1, 2, 1]] / 16
    // FIXME: Feel free to optimize this function
    // Hint: You can use SIMD instructions to optimize this function
    void gaussianFilter()
    {
        im2colInit();
        // 仅初始化 im2col 就已经难敌常规优化，因此不再继续往下写
    }

    // Power law transformation
    // FIXME: Feel free to optimize this function
    // Hint: LUT to optimize this function?
    void powerLawTransformation()
    {
        constexpr float gamma = 0.5f;

        for (size_t i = 0; i < buffer_size; ++i)
        {
            if (figure[i] == 0)
            {
                result[i] = 0;
                continue;
            }
            float normalized = figure[i] / 255.0f;
            result[i] = static_cast<unsigned char>(
                255.0f * std::pow(normalized, gamma) + 0.5f);
        }
    }

    // Run benchmark
    unsigned int calcChecksum()
    {
        unsigned int sum = 0;
        constexpr size_t mod = 1000000007;
        for (size_t i = 0; i < buffer_size; ++i)
        {
            sum += result[i];
            sum %= mod;
        }
        return sum;
    }
    void runBenchmark()
    {
#ifdef DEBUG
        print(figure, size, size);
        std::cout << "\n";
#endif
        auto start = std::chrono::high_resolution_clock::now();
        gaussianFilter();
        auto middle = std::chrono::high_resolution_clock::now();
#ifdef DEBUG
        // print(result, size, size);
        print(im2col_buffer, buffer_size, 9);
        std::cout << "\n";
#endif
        unsigned int sum = calcChecksum();
        auto middle2 = std::chrono::high_resolution_clock::now();
        powerLawTransformation();
        auto end = std::chrono::high_resolution_clock::now();
#ifdef DEBUG
        print(result, size, size);
        std::cout << "\n";
#endif
        sum += calcChecksum();
        sum %= 1000000007;
        std::cout << "Checksum: " << sum << "\n";
        auto milliseconds_1 = std::chrono::duration_cast<std::chrono::milliseconds>(middle - start);
        auto milliseconds_2 = std::chrono::duration_cast<std::chrono::milliseconds>(end - middle2);
        auto milliseconds = milliseconds_1 + milliseconds_2;
#ifdef BREAKDOWN_TIME
        std::cout << "Benchmark time 1: " << milliseconds_1.count() << " ms\n";
        std::cout << "Benchmark time 2: " << milliseconds_2.count() << " ms\n";
#endif
        std::cout << "Benchmark time: " << milliseconds.count() << " ms\n";
    }
};

// Main function
// !!! Please do not modify the main function !!!
int main(int argc, const char **argv)
{
#ifdef DEBUG
    constexpr size_t size = 3;
#else
    constexpr size_t size = 16384;
#endif
    FigureProcessor processor(size, argc > 1 ? std::stoul(argv[1]) : 0);
    processor.runBenchmark();
    return 0;
}
