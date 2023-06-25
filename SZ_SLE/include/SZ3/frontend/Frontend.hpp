#ifndef SZ3_FRONTEND_INTERFACE
#define SZ3_FRONTEND_INTERFACE
/**
 * Frontend is the combination of Predictor and Quantizer
 * For compression, it takes the original data as input, and outputs integer values
 * which will be used for lossless compression by the Encoder and Lossless modules.
 */

#include "SZ3/def.hpp"
#include <vector>

namespace SZ {


    namespace concepts {

        template<class T, uint N>
        class FrontendInterface {
        public:

            virtual ~FrontendInterface() = default;

            virtual std::vector<int> compress(T *data, size_t x, size_t y, size_t z, size_t szBlk) = 0;
            // virtual std::vector<int> compress(T *data) = 0;

            // virtual T *decompress(std::vector<int> &quant_inds, T *dec_data, size_t x, size_t y, size_t z) = 0;
            virtual T *decompress(std::vector<int> &quant_inds, T *dec_data) = 0;

            virtual int *save(uchar *&c, size_t &regCnt, size_t szBlk) = 0;

            virtual void load(const uchar *&c, size_t &remaining_length, size_t x, size_t y, size_t z, size_t &regCnt, std::vector<int> &quant_inds, bool real) = 0;
            // virtual void load(const uchar *&c, size_t &remaining_length) = 0;

            virtual size_t size_est() = 0;

            virtual int get_radius() const = 0;

            virtual size_t get_blkSize() = 0;

            virtual size_t get_blkNum() const = 0;

            virtual size_t *get_meta() = 0;

            virtual size_t get_ifTree() = 0;

            virtual size_t get_num_elements() const = 0;

            virtual void print() = 0;

            virtual void clear() = 0;
        };

    }

}

#endif
