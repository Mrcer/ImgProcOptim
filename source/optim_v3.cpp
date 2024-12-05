/*
 * + 使用 AVX2 指令集优化高斯滤波器和查表
 */
#include <chrono>
#include <cstddef>
#include <iostream>
#include <random>
#include <cstdlib>
#include <immintrin.h>
#include <cstdint>

#define ELEMENT(buffer, i, j, dim1) (buffer[(i) * (dim1) + (j)])

void print(uint16_t *buffer, size_t dim0, size_t dim1)
{
    for (size_t i = 0; i < dim0; ++i)
    {
        for (size_t j = 0; j < dim1; ++j)
        {
            std::cout << static_cast<int>(ELEMENT(buffer, i, j, dim1)) << " ";
        }
        std::cout << std::endl;
    }
}

class FigureProcessor
{
private:
    uint16_t *figure;
#define ELEMENT_FIGURE(i, j) ELEMENT(figure, i + 1, j + 1, size + 2)
    uint16_t *result;
    const size_t size;
    const size_t buffer_size;

    struct _POWER_LAW_LUT {
        constexpr _POWER_LAW_LUT(): table() {
            for (int i = 0; i < 256; ++i) {
                table[i] = static_cast<unsigned>(255 * sqrt(i / 255.0) + 0.5);
            }
        }
        constexpr float sqrt(float x) {
            float val = x;
            float last = 0;
            do {
                last = val;
                val = (val + x / val) / 2;
            } while (abs(val - last) > 0.0001);
            return val;
        }
        constexpr float abs(float x) const {
            return (x < 0)? -x : x;
        }
        unsigned table[256];
    } POWER_LAW_LUT;

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
        figure = static_cast<uint16_t *>(malloc(sizeof(uint16_t) * (buffer_size + 4 * (size + 1))));    // 增加 padding
        result = static_cast<uint16_t *>(malloc(sizeof(uint16_t) * (buffer_size)));

#ifdef DEBUG
        int _debug_cnt = 0;
#define GENRAND() (_debug_cnt++ % 7)
#else
#define GENRAND() static_cast<uint16_t>(distribution(gen))
#endif
#define APPEND(x) *figure_ptr = (x), figure_ptr += 1;
#define SET_LASTROW(x) *(figure_ptr - size - 2) = (x);
#define SET_NEXTROW(x) *(figure_ptr + size + 2) = (x);
        auto figure_ptr = figure;
        uint16_t rand;
        APPEND(0)
        figure_ptr += size;
        APPEND(0)
        rand = GENRAND();
        APPEND(rand)
        SET_LASTROW(rand)
        APPEND(rand)
        for(size_t i = 1; i < size - 1; ++i)
        {
            auto rand = GENRAND();
            SET_LASTROW(rand)
            APPEND(rand)
        }
        rand = GENRAND();
        SET_LASTROW(rand)
        APPEND(rand)
        APPEND(rand)

        for(size_t i = 1; i < size - 1; ++i)
        {
            rand = GENRAND();
            APPEND(rand)
            APPEND(rand)
            for(size_t j = 1 ; j < size - 1 ; ++j)
            {
                APPEND(GENRAND())
            }
            rand = GENRAND();
            APPEND(rand)
            APPEND(rand)
        }

        SET_NEXTROW(0)
        rand = GENRAND();
        APPEND(rand)
        SET_NEXTROW(rand)
        APPEND(rand)
        for(size_t i = 1; i < size - 1; ++i)
        {
            rand = GENRAND();
            SET_NEXTROW(rand)
            APPEND(rand)
        }
        rand = GENRAND();
        SET_NEXTROW(rand)
        APPEND(rand)
        SET_NEXTROW(0)
        APPEND(rand)
#undef GENRAND
#undef APPEND
#undef SET_LASTROW
#undef SET_NEXTROW
    }

    ~FigureProcessor()
    {
        free(figure);
        free(result);
    }

    // 使用 AVX2 指令集对以 src 为左上角的 3x18 元素计算高斯，结果为 16 个元素存入 dst
    void gaussianFilterBlock(uint16_t *src, uint16_t *dst)
    {
        __v16hi a2,a3,b1,b2,b3,c1,c2,c3,result;
        // a1 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src));
        result = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src));
        b1 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + 1));
        c1 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + 2));
        a2 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + size + 2));
        b2 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + size + 3));
        c2 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + size + 4));
        a3 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + 2 * (size + 2)));
        b3 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + 2 * (size + 2) + 1));
        c3 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + 2 * (size + 2) + 2));

        // result = _mm256_xor_ps(result, result);
        // result = _mm256_add_epi16(result, a1);  // result += K[0][0]
        b1 = _mm256_slli_epi16(b1, 1);
        result = _mm256_add_epi16(result, b1);  // result += 2 * K[0][1]
        result = _mm256_add_epi16(result, c1);  // result += K[0][0]
        a2 = _mm256_slli_epi16(a2, 1);
        result = _mm256_add_epi16(result, a2);  // result += 2 * K[1][0]
        b2 = _mm256_slli_epi16(b2, 2);
        result = _mm256_add_epi16(result, b2);  // result += 4 * K[1][1]
        c2 = _mm256_slli_epi16(c2, 1);
        result = _mm256_add_epi16(result, c2);  // result += 2 * K[1][2]
        result = _mm256_add_epi16(result, a3);  // result += K[2][0]
        b3 = _mm256_slli_epi16(b3, 1);
        result = _mm256_add_epi16(result, b3);  // result += 2 * K[2][1]
        result = _mm256_add_epi16(result, c3);  // result += K[2][2]
        result = _mm256_srli_epi16(result, 4);           // result /= 16
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(dst), result);
    }

    // Gaussian filter
    // [[1, 2, 1], [2, 4, 2], [1, 2, 1]] / 16
    // FIXME: Feel free to optimize this function
    // Hint: You can use SIMD instructions to optimize this function
    void gaussianFilter()
    {
        if(size % 16 != 0 || size < 16) {
            std::cerr << "Error: size must be a multiple of 16 and greater than or equal to 16\n";
            return;
        }
        
        for (size_t i = 0; i < size; ++i)
        {
            auto src = &ELEMENT(figure, i, 0, size + 2);
            auto dst = &ELEMENT(result, i, 0, size);
            for (size_t j = 0; j < size; j += 16)
            {
                gaussianFilterBlock(src, dst);
                src += 16;
                dst += 16;
            }
        }

        // 处理四个角点
        // 左上角
        ELEMENT(result, 0, 0, size) = (4 * ELEMENT_FIGURE(0, 0) + 2 * ELEMENT_FIGURE(0, 1) +
                                       2 * ELEMENT_FIGURE(1, 0) + ELEMENT_FIGURE(1, 1)) /
                                      9;
        // 右上角
        ELEMENT(result, 0, size - 1, size) =
            (4 * ELEMENT_FIGURE(0, size - 1) + 2 * ELEMENT_FIGURE(0, size - 2) +
             2 * ELEMENT_FIGURE(1, size - 1) + ELEMENT_FIGURE(1, size - 2)) /
            9;

        // 左下角
        ELEMENT(result, size - 1, 0, size) =
            (4 * ELEMENT_FIGURE(size - 1, 0) + 2 * ELEMENT_FIGURE(size - 1, 1) +
             2 * ELEMENT_FIGURE(size - 2, 0) + ELEMENT_FIGURE(size - 2, 1)) /
            9;

        // 右下角
        ELEMENT(result, size - 1, size - 1, size) =
            (4 * ELEMENT_FIGURE(size - 1, size - 1) + 2 * ELEMENT_FIGURE(size - 1, size - 2) +
             2 * ELEMENT_FIGURE(size - 2, size - 1) + ELEMENT_FIGURE(size - 2, size - 2)) /
            9;
    }

    // Power law transformation
    // FIXME: Feel free to optimize this function
    // Hint: LUT to optimize this function?
    void powerLawTransformation()
    {
        if(size % 16 != 0 || size < 16) {
            std::cerr << "Error: size must be a multiple of 16 and greater than or equal to 16\n";
            return;
        }

        for (size_t i = 0; i < size; ++i)
        {
            uint16_t *result_ptr = result + i * size;
            uint16_t *figure_ptr = &ELEMENT_FIGURE(i, 0);
            auto end = result_ptr + size;
            const __v8su zero = _mm256_setzero_si256();
            while(result_ptr < end)
            {
                __v16hi load_reg, store_reg;
                __v8su idx_hi, idx_lo, val_hi, val_lo;
                unsigned *lut = POWER_LAW_LUT.table;
                load_reg = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(figure_ptr));
                idx_lo = _mm256_unpacklo_epi16(load_reg, zero);
                idx_hi = _mm256_unpackhi_epi16(load_reg, zero);
                val_lo = _mm256_i32gather_epi32(lut, idx_lo, 4);
                val_hi = _mm256_i32gather_epi32(lut, idx_hi, 4);
                store_reg = _mm256_packus_epi32(val_lo, val_hi);
                _mm256_storeu_si256(reinterpret_cast<__m256i *>(result_ptr), store_reg);
                result_ptr += 16;
                figure_ptr += 16;
            }
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
        print(figure, size + 2, size + 2);
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
    constexpr size_t size = 16;
#else
    constexpr size_t size = 16384;
#endif
    FigureProcessor processor(size, argc > 1 ? std::stoul(argv[1]) : 0);
    processor.runBenchmark();
    return 0;
}
