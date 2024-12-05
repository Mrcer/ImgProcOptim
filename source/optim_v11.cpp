/* Warnning: 该代码存在一个 bug，仅用于演示性能
 * + 修改 gaussian filter 算法，尝试减少访存和运算量
 */
#include <chrono>
#include <cstddef>
#include <iostream>
#include <random>

#include <cstdint>
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
    uint16_t *result_buffer;

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
        result_buffer = static_cast<uint16_t *>(malloc(buffer_size * 2));

        for (size_t i = 0; i < buffer_size; ++i)
        {
            figure[i] = static_cast<unsigned char>(distribution(gen));
            result_buffer[i] = 0;
        }
    }

    ~FigureProcessor()
    {
        free(figure);
        free(result);
    }

    // 累加一行卷积核计算结果（不包括除法）
    void gaussianFilterStepAdd(unsigned char *src, uint16_t *dst, const size_t iter_i, const unsigned char factor)
    {
        for (int i = 0; i < iter_i; ++i)
        {
            // 最左列元素
            unsigned char a;
            unsigned char b = *src;
            unsigned char c = *(src + 1);
            *dst += (3 * b + c) * factor; // a + 2 * b + c = 3 * b + c, when a = b

            uint16_t *dst_end = dst + size - 1;
            for (src += 2, dst += 1; dst < dst_end; src += 1, dst += 1)
            {
                a = b;
                b = c;
                c = *src;
                *dst += (a + 2 * b + c) * factor;
            }

            // 最右列元素
            a = b;
            b = c;
            *dst += (a + 3 * b) * factor; // a + 2 * b + c = a + 3 * b, when c = b

            // 保证 src 和 dst 指针指向下一行第一元素
            dst += 1;
        }
    }

    // 角点除以 9，其他元素除以 16
    void gaussianFilterStepDiv()
    {
        uint16_t *src = result_buffer;
        unsigned char *dst = result;
        unsigned char *dst_end;

#define DIV_9        \
    *dst = *src / 9; \
    dst += 1;        \
    src += 1;

#define DIV_16                                \
    for (; dst < dst_end; dst += 1, src += 1) \
    {                                         \
        *dst = *src / 16;                     \
    }

        // 左上角
        DIV_9

        // 第一列
        dst_end = dst + size - 2;
        DIV_16

        // 右上角
        DIV_9

        // 中间
        dst_end = dst + size * (size - 2);
        DIV_16

        // 左下角
        DIV_9

        // 最后一列
        dst_end = dst + size - 2;
        DIV_16

        // 右下角
        DIV_9
    }

    // Gaussian filter
    // [[1, 2, 1], [2, 4, 2], [1, 2, 1]] / 16
    // FIXME: Feel free to optimize this function
    // Hint: You can use SIMD instructions to optimize this function
    void gaussianFilter()
    {
        // Bug：第一行和最后一行少加了 figure 外元素
        // 卷积核第一行
        gaussianFilterStepAdd(figure, result_buffer + size, size - 1, 1);
        // 卷积核第二行
        gaussianFilterStepAdd(figure, result_buffer, size, 2);
        // 卷积核第三行
        gaussianFilterStepAdd(figure + size, result_buffer, size - 1, 1);

        gaussianFilterStepDiv();
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
