/*
 * + 使用指针分配内存，避免 vector 重复扩容开销和二级指针访存
 */
#include <chrono>
#include <cstddef>
#include <iostream>
#include <random>
#include <cmath>

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

    // Gaussian filter
    // [[1, 2, 1], [2, 4, 2], [1, 2, 1]] / 16
    // FIXME: Feel free to optimize this function
    // Hint: You can use SIMD instructions to optimize this function
    void gaussianFilter()
    {
        // 处理内部区域
        for (size_t i = 1; i < size - 1; ++i)
        {
            for (size_t j = 1; j < size - 1; ++j)
            {
                ELEMENT(result, i, j, size) =
                    (ELEMENT(figure, i - 1, j - 1, size) + 2 * ELEMENT(figure, i - 1, j, size) +
                     ELEMENT(figure, i - 1, j + 1, size) + 2 * ELEMENT(figure, i, j - 1, size) +
                     4 * ELEMENT(figure, i, j, size) + 2 * ELEMENT(figure, i, j + 1, size) + ELEMENT(figure, i + 1, j - 1, size) +
                     2 * ELEMENT(figure, i + 1, j, size) + ELEMENT(figure, i + 1, j + 1, size)) /
                    16;
            }
        }

        for (size_t i = 1; i < size - 1; ++i)
        {
            ELEMENT(result, i, 0, size) =
                (ELEMENT(figure, i - 1, 0, size) + 2 * ELEMENT(figure, i - 1, 0, size) +
                 ELEMENT(figure, i - 1, 1, size) + 2 * ELEMENT(figure, i, 0, size) +
                 4 * ELEMENT(figure, i, 0, size) + 2 * ELEMENT(figure, i, 1, size) +
                 ELEMENT(figure, i + 1, 0, size) + 2 * ELEMENT(figure, i + 1, 0, size) +
                 ELEMENT(figure, i + 1, 1, size)) /
                16;

            ELEMENT(result, i, size - 1, size) =
                (ELEMENT(figure, i - 1, size - 2, size) + 2 * ELEMENT(figure, i - 1, size - 1, size) +
                 ELEMENT(figure, i - 1, size - 1, size) + 2 * ELEMENT(figure, i, size - 2, size) +
                 4 * ELEMENT(figure, i, size - 1, size) + 2 * ELEMENT(figure, i, size - 1, size) +
                 ELEMENT(figure, i + 1, size - 2, size) + 2 * ELEMENT(figure, i + 1, size - 1, size) +
                 ELEMENT(figure, i + 1, size - 1, size)) /
                16;
        }

        for (size_t j = 1; j < size - 1; ++j)
        {
            ELEMENT(result, 0, j, size) =
                (ELEMENT(figure, 0, j - 1, size) + 2 * ELEMENT(figure, 0, j, size) +
                 ELEMENT(figure, 0, j + 1, size) + 2 * ELEMENT(figure, 0, j - 1, size) +
                 4 * ELEMENT(figure, 0, j, size) + 2 * ELEMENT(figure, 0, j + 1, size) +
                 ELEMENT(figure, 1, j - 1, size) + 2 * ELEMENT(figure, 1, j, size) +
                 ELEMENT(figure, 1, j + 1, size)) /
                16;

            ELEMENT(result, size - 1, j, size) =
                (ELEMENT(figure, size - 2, j - 1, size) + 2 * ELEMENT(figure, size - 2, j, size) +
                 ELEMENT(figure, size - 2, j + 1, size) + 2 * ELEMENT(figure, size - 1, j - 1, size) +
                 4 * ELEMENT(figure, size - 1, j, size) + 2 * ELEMENT(figure, size - 1, j + 1, size) +
                 ELEMENT(figure, size - 1, j - 1, size) + 2 * ELEMENT(figure, size - 1, j, size) +
                 ELEMENT(figure, size - 1, j + 1, size)) /
                16;
        }

        // 处理四个角点
        // 左上角
        ELEMENT(result, 0, 0, size) = (4 * ELEMENT(figure, 0, 0, size) + 2 * ELEMENT(figure, 0, 1, size) +
                                       2 * ELEMENT(figure, 1, 0, size) + ELEMENT(figure, 1, 1, size)) /
                                      9;
        // 右上角
        ELEMENT(result, 0, size - 1, size) =
            (4 * ELEMENT(figure, 0, size - 1, size) + 2 * ELEMENT(figure, 0, size - 2, size) +
             2 * ELEMENT(figure, 1, size - 1, size) + ELEMENT(figure, 1, size - 2, size)) /
            9;

        // 左下角
        ELEMENT(result, size - 1, 0, size) =
            (4 * ELEMENT(figure, size - 1, 0, size) + 2 * ELEMENT(figure, size - 1, 1, size) +
             2 * ELEMENT(figure, size - 2, 0, size) + ELEMENT(figure, size - 2, 1, size)) /
            9;

        // 右下角
        ELEMENT(result, size - 1, size - 1, size) =
            (4 * ELEMENT(figure, size - 1, size - 1, size) + 2 * ELEMENT(figure, size - 1, size - 2, size) +
             2 * ELEMENT(figure, size - 2, size - 1, size) + ELEMENT(figure, size - 2, size - 2, size)) /
            9;
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
    // Run benchmark
    unsigned int calcStrongChecksum()
    {
        unsigned int sum = 0;
        constexpr size_t mod = 1000000007;
        for (size_t i = 0; i < buffer_size; ++i)
        {
            sum += result[i] * int(10 * sin(i/97.0));
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
        print(result, size, size);
        std::cout << "\n";
#endif
        unsigned int sum = calcChecksum();
        unsigned int strongSum = calcStrongChecksum();
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
        strongSum += calcStrongChecksum();
        strongSum %= 1000000007;
#ifdef STRONG_CHECKSUM
        std::cout << "Strong checksum: " << strongSum << "\n";
#endif
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
