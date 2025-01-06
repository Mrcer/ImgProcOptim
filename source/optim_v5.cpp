#include <chrono>
#include <cstddef>
#include <iostream>
#include <random>
#include <cstdlib>
#include <immintrin.h>
#include <cstdint>

#define ELEMENT(buffer, i, j, dim1) (buffer[(i) * (dim1) + (j)])

class FigureProcessor
{
private:
    uint8_t *figure;
    uint8_t *result;
    const size_t size;
    const size_t buffer_size;
#define ELEMENT_FIGURE(i, j) ELEMENT(figure, i + 1, j + 1, size + 2)
#define ELEMENT_RESULT(i, j) ELEMENT(result, i + 1, j + 1, size + 2)

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
    FigureProcessor(size_t size, size_t seed = 0) : size(size), buffer_size((size * size + 4 * (size + 1)))
    {
        if(size % 16 != 0 || size < 128)
        {
            std::cerr << "Error: size must be a multiple of 16 and >= 128" << std::endl;
            exit(1);
        }

        // !!! Please do not modify the following code !!!
        // 如果你需要修改内存的数据结构，请不要修改初始化的顺序和逻辑
        // 助教可能会通过指定某个初始化seed 的方式来验证你的代码
        // 如果你修改了初始化的顺序，可能会导致你的代码无法通过测试
        std::random_device rd;
        std::mt19937_64 gen(seed == 0 ? rd() : seed);
        std::uniform_int_distribution<unsigned char> distribution(0, 255);
        // !!! ----------------------------------------- !!!

        // 两个数组的初始化在这里，可以改动，但请注意 gen 的顺序是从上到下从左到右即可。
        figure = static_cast<uint8_t *>(malloc(sizeof(uint8_t) * buffer_size));    // 增加 padding
        result = static_cast<uint8_t *>(malloc(sizeof(uint8_t) * buffer_size));

#define GENRAND() static_cast<uint8_t>(distribution(gen))
#define APPEND(x) *figure_ptr = (x), figure_ptr += 1;
#define SET_LASTROW(x) *(figure_ptr - size - 2) = (x);
#define SET_NEXTROW(x) *(figure_ptr + size + 2) = (x);
        auto figure_ptr = figure;
        uint8_t rand;
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
    inline void gaussianFilterBlock(uint8_t *src, uint8_t *dst)
    {
        __v16hu a2,a3,b1,b2,b3,c1,c2,c3,result;
        __m128i load_buffer;

        // a1 is assigned to result
        load_buffer = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src));  // SSE2
        result = _mm256_cvtepu8_epi16(load_buffer);
        load_buffer = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src + 1));
        b1 = _mm256_cvtepu8_epi16(load_buffer);
        load_buffer = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src + 2));
        c1 = _mm256_cvtepu8_epi16(load_buffer);
        load_buffer = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src + size + 2));
        a2 = _mm256_cvtepu8_epi16(load_buffer);
        load_buffer = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src + size + 3));
        b2 = _mm256_cvtepu8_epi16(load_buffer);
        load_buffer = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src + size + 4));
        c2 = _mm256_cvtepu8_epi16(load_buffer);
        load_buffer = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src + ((size + 2) << 1)));
        a3 = _mm256_cvtepu8_epi16(load_buffer);
        load_buffer = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src + ((size + 2) << 1) + 1));
        b3 = _mm256_cvtepu8_epi16(load_buffer);
        load_buffer = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src + ((size + 2) << 1) + 2));
        c3 = _mm256_cvtepu8_epi16(load_buffer);

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
        result = _mm256_srli_epi16(result, 4);  // result /= 16
        
        static const __m256i zero = _mm256_setzero_si256();     // AVX
        __v32qu buffer256 = _mm256_packus_epi16(result, zero);
        buffer256 = _mm256_permute4x64_epi64(buffer256, 0xD8);     // 3120 => 0xD8
        __v16qu buffer128 = _mm256_extracti128_si256(buffer256, 0);
        _mm_storeu_si128(reinterpret_cast<__m128i *>(dst), buffer128); // SSE2
    }

    // Gaussian filter
    // [[1, 2, 1], [2, 4, 2], [1, 2, 1]] / 16
    // FIXME: Feel free to optimize this function
    // Hint: You can use SIMD instructions to optimize this function
    void gaussianFilter()
    {
        #pragma omp parallel for schedule(static) proc_bind(spread) num_threads(4)
        for (size_t i = 0; i < size; ++i)
        {
            auto src = &ELEMENT(figure, i, 0, size + 2);
            auto dst = &ELEMENT_RESULT(i, 0);
            for (size_t j = 0; j < size; j += 16)
            {
                gaussianFilterBlock(src, dst);
                src += 16;
                dst += 16;
            }
        }

        // 处理四个角点
        // 左上角
        ELEMENT_RESULT(0, 0) = ((ELEMENT_FIGURE(0, 0) << 2) + (ELEMENT_FIGURE(0, 1) << 1) +
                                       (ELEMENT_FIGURE(1, 0) << 1) + ELEMENT_FIGURE(1, 1)) /
                                      9;
        // 右上角
        ELEMENT_RESULT(0, size - 1) =
            ((ELEMENT_FIGURE(0, size - 1) << 2) + (ELEMENT_FIGURE(0, size - 2) << 1) +
             (ELEMENT_FIGURE(1, size - 1) << 1) + ELEMENT_FIGURE(1, size - 2)) /
            9;

        // 左下角
        ELEMENT_RESULT(size - 1, 0) =
            ((ELEMENT_FIGURE(size - 1, 0) << 2) + (ELEMENT_FIGURE(size - 1, 1) << 1) +
             (ELEMENT_FIGURE(size - 2, 0) << 1) + ELEMENT_FIGURE(size - 2, 1)) /
            9;

        // 右下角
        ELEMENT_RESULT(size - 1, size - 1) =
            ((ELEMENT_FIGURE(size - 1, size - 1) << 2) + (ELEMENT_FIGURE(size - 1, size - 2) << 1) +
             (ELEMENT_FIGURE(size - 2, size - 1) << 1) + ELEMENT_FIGURE(size - 2, size - 2)) /
            9;
    }

    // Power law transformation
    // FIXME: Feel free to optimize this function
    // Hint: LUT to optimize this function?
    void powerLawTransformation()
    {
        uint8_t *figure_valid_aligned_begin = reinterpret_cast<uint8_t *>(reinterpret_cast<size_t>(&ELEMENT_FIGURE(0, 0)) & ~0xF);
        uint8_t *result_valid_aligned_begin = reinterpret_cast<uint8_t *>(reinterpret_cast<size_t>(&ELEMENT_RESULT(0, 0)) & ~0xF);
        uint8_t *result_valid_aligned_end = reinterpret_cast<uint8_t *>((reinterpret_cast<size_t>(&ELEMENT_RESULT(size - 1, size - 1)) & ~0xF) + 16);
        size_t loop_times = (result_valid_aligned_end - result_valid_aligned_begin) >> 4;
        #pragma omp parallel for schedule(static) proc_bind(spread) num_threads(4)
        for (size_t i = 0; i < loop_times; ++i)
        {
            static const __m256i zero = _mm256_setzero_si256();
            size_t offset = i << 4;
            uint8_t *result_ptr = result_valid_aligned_begin + offset;
            uint8_t *figure_ptr = figure_valid_aligned_begin + offset;
            __v16qu load_reg;
            __v8su idx_hi, idx_lo, val_hi, val_lo;
            __v16hu val;
            unsigned *lut = POWER_LAW_LUT.table;

            load_reg = _mm_stream_load_si128(reinterpret_cast<const __m128i *>(figure_ptr));   // SSE4.1
            idx_lo = _mm256_cvtepu8_epi32(load_reg);
            load_reg = _mm_srli_si128(load_reg, 8);     // SSE2
            idx_hi = _mm256_cvtepu8_epi32(load_reg);
            val_lo = _mm256_i32gather_epi32(lut, idx_lo, 4);
            val_hi = _mm256_i32gather_epi32(lut, idx_hi, 4);
            val = _mm256_packus_epi32(val_lo, val_hi);
            val = _mm256_permute4x64_epi64(val, 0xD8);     // Aware that _mm256_packus_epi32 runs in two parts

            __v32qu buffer256 = _mm256_packus_epi16(val, zero);
            buffer256 = _mm256_permute4x64_epi64(buffer256, 0xD8);     // 3120 => 0xD8
            __v16qu buffer128 = _mm256_extracti128_si256(buffer256, 0);
            _mm_stream_si128(reinterpret_cast<__m128i *>(result_ptr), buffer128); // SSE2
        }
    }

    // Run benchmark
    unsigned int calcChecksum()
    {
        unsigned int sum = 0;
        constexpr size_t mod = 1000000007;
        for (size_t i = 0; i < size; ++i)
        {
            for (size_t j = 0; j < size; ++j)
            {
                sum += ELEMENT_RESULT(i, j);
                sum %= mod;
            }
        }
        return sum;
    }
    
    void runBenchmark()
    {
        auto start = std::chrono::high_resolution_clock::now();
        gaussianFilter();
        auto middle = std::chrono::high_resolution_clock::now();

        unsigned int sum = calcChecksum();

        auto middle2 = std::chrono::high_resolution_clock::now();
        powerLawTransformation();
        auto end = std::chrono::high_resolution_clock::now();

        sum += calcChecksum();
        sum %= 1000000007;
        std::cout << "Checksum: " << sum << "\n";

        auto milliseconds =
            std::chrono::duration_cast<std::chrono::milliseconds>(middle - start) +
            std::chrono::duration_cast<std::chrono::milliseconds>(end - middle2);
        std::cout << "Benchmark time: " << milliseconds.count() << " ms\n";
    }
};

// Main function
// !!! Please do not modify the main function !!!
int main(int argc, const char **argv)
{
    constexpr size_t size = 16384;
    FigureProcessor processor(size, argc > 1 ? std::stoul(argv[1]) : 0);
    processor.runBenchmark();
    return 0;
}
