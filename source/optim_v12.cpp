/*
 * + 基于 v11，修改 gaussian filter 算法，减少访存和运算量
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
        unsigned char *src1;
        unsigned char *src2 = figure;
        unsigned char *src3 = figure + size;
        unsigned char *dst = result;
        unsigned char *dst_end;
        int a1,a2,a3,b1,b2,b3,c1,c2,c3;     // 使用寄存器避免重复访存（如果有足够的寄存器），abc表示列，123表示行

        // 第一行 ------------------
        // 左上角元素
        b2 = *src2; ++src2; c2 = *src2; ++src2;
        b3 = *src3; ++src3; c3 = *src3; ++src3;
        *dst = (4 * b2 + 2 * c2 + 2 * b3 + c3) / 9;
        ++dst;
        // 第一行中间元素
        dst_end = dst + size - 2;
        while (dst != dst_end)
        {
            a2 = b2; b2 = c2; c2 = *src2; ++src2;
            a3 = b3; b3 = c3; c3 = *src3; ++src3;
            // 根据原代码，卷积核第一行用第二行元素重复，因此合并，下面不再赘述
            *dst = (3 * a2 + 6 * b2 + 3 * c2 + a3 + 2 * b3 + c3) / 16;
            ++dst;
        }
        // 右上角元素
        a2 = b2; b2 = c2;
        a3 = b3; b3 = c3;
        *dst = (2 * a2 + 4 * b2 + a3 + 2 * b3) / 9;
        ++dst;

        // 中间行 -----------------
        src1 = figure;  // 该循环每次迭代前所有 src 均指向 figure 的第一列
        for(int i = 1 ; i < size - 1 ; ++i)
        {
            // 最左列元素
            b1 = *src1; ++src1; c1 = *src1; ++src1;
            b2 = *src2; ++src2; c2 = *src2; ++src2;
            b3 = *src3; ++src3; c3 = *src3; ++src3;
            *dst = (3 * b1 + c1 + 6 * b2 + 2 * c2 + 3 * b3 + c3) / 16;
            ++dst;

            // 中间行中间列元素
            dst_end = dst + size - 2;
            while (dst != dst_end)
            {
                a1 = b1; b1 = c1; c1 = *src1; ++src1;
                a2 = b2; b2 = c2; c2 = *src2; ++src2;
                a3 = b3; b3 = c3; c3 = *src3; ++src3;
                *dst = (a1 + 2 * b1 + c1 + 2 * a2 + 4 * b2 + 2 * c2 + a3 + 2 * b3 + c3) / 16;
                ++dst;
            }

            // 最右列元素
            a1 = b1; b1 = c1;
            a2 = b2; b2 = c2;
            a3 = b3; b3 = c3;
            *dst = (a1 + 3 * b1 + 2 * a2 + 6 * b2 + a3 + 3 * b3) / 16;
            ++dst;
        }

        // 最后一行 ---------------
        // 左下角元素
        b1 = *src1; ++src1; c1 = *src1; ++src1;
        b2 = *src2; ++src2; c2 = *src2; ++src2;
        *dst = (2 * b1 + c1 + 4 * b2 + 2 * c2) / 9;
        ++dst;
        // 最后一行中间元素
        dst_end = dst + size - 2;
        while (dst != dst_end)
        {
            a1 = b1; b1 = c1; c1 = *src1; ++src1;
            a2 = b2; b2 = c2; c2 = *src2; ++src2;
            *dst = (a1 + 2 * b1 + c1 + 3 * a2 + 6 * b2 + 3 * c2) / 16;
            ++dst;
        }
        // 右下角元素
        a1 = b1; b1 = c1;
        a2 = b2; b2 = c2;
        *dst = (a1 + 2 * b1 + 2 * a2 + 4 * b2) / 9;
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
        print(result, size, size);
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
