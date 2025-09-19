#include <bits/stdc++.h>
#include <omp.h>
using namespace std;

// SM4 S-box 
const uint8_t sbox[256] = {
    0xd6, 0x90, 0xe9, 0xfe, 0xcc, 0xe1, 0x3d, 0xb7, 0x16, 0xb6, 0x14, 0xc2, 0x28, 0xfb, 0x2c, 0x05,
    0x2b, 0x67, 0x9a, 0x76, 0x2a, 0xbe, 0x04, 0xc3, 0xaa, 0x44, 0x13, 0x26, 0x49, 0x86, 0x06, 0x99,
    0x9c, 0x42, 0x50, 0xf4, 0x91, 0xef, 0x98, 0x7a, 0x33, 0x54, 0x0b, 0x43, 0xed, 0xcf, 0xac, 0x62,
    0xe4, 0xb3, 0x1c, 0xa9, 0xc9, 0x08, 0xe8, 0x95, 0x80, 0xdf, 0x94, 0xfa, 0x75, 0x8f, 0x3f, 0xa6,
    0x47, 0x07, 0xa7, 0xfc, 0xf3, 0x73, 0x17, 0xba, 0x83, 0x59, 0x3c, 0x19, 0xe6, 0x85, 0x4f, 0xa8,
    0x68, 0x6b, 0x81, 0xb2, 0x71, 0x64, 0xda, 0x8b, 0xf8, 0xeb, 0x0f, 0x4b, 0x70, 0x56, 0x9d, 0x35,
    0x1e, 0x24, 0x0e, 0x5e, 0x63, 0x58, 0xd1, 0xa2, 0x25, 0x22, 0x7c, 0x3b, 0x01, 0x21, 0x78, 0x87,
    0xd4, 0x00, 0x46, 0x57, 0x9f, 0xd3, 0x27, 0x52, 0x4c, 0x36, 0x02, 0xe7, 0xa0, 0xc4, 0xc8, 0x9e,
    0xea, 0xbf, 0x8a, 0xd2, 0x40, 0xc7, 0x38, 0xb5, 0xa3, 0xf7, 0xf2, 0xce, 0xf9, 0x61, 0x15, 0xa1,
    0xe0, 0xae, 0x5d, 0xa4, 0x9b, 0x34, 0x1a, 0x55, 0xad, 0x93, 0x32, 0x30, 0xf5, 0x8c, 0xb1, 0xe3,
    0x1d, 0xf6, 0xe2, 0x2e, 0x82, 0x66, 0xca, 0x60, 0xc0, 0x29, 0x23, 0xab, 0x0d, 0x53, 0x4e, 0x6f,
    0xd5, 0xdb, 0x37, 0x45, 0xde, 0xfd, 0x8e, 0x2f, 0x03, 0xff, 0x6a, 0x72, 0x6d, 0x6c, 0x5b, 0x51,
    0x8d, 0x1b, 0xaf, 0x92, 0xbb, 0xdd, 0xbc, 0x7f, 0x11, 0xd9, 0x5c, 0x41, 0x1f, 0x10, 0x5a, 0xd8,
    0x0a, 0xc1, 0x31, 0x88, 0xa5, 0xcd, 0x7b, 0xbd, 0x2d, 0x74, 0xd0, 0x12, 0xb8, 0xe5, 0xb4, 0xb0,
    0x89, 0x69, 0x97, 0x4a, 0x0c, 0x96, 0x77, 0x7e, 0x65, 0xb9, 0xf1, 0x09, 0xc5, 0x6e, 0xc6, 0x84,
    0x18, 0xf0, 0x7d, 0xec, 0x3a, 0xdc, 0x4d, 0x20, 0x79, 0xee, 0x5f, 0x3e, 0xd7, 0xcb, 0x39, 0x48
};

// Fixed key constants so that we can use in key expansion for initial XOR
const uint32_t FK[4] = {0xa3b1bac6, 0x56aa3350, 0x677d9197, 0xb27022dc};

// Round constants are precomputed for each round in key/enc
uint32_t CK[32];

// LAT for S-box
int lat[256][256];

// Left rotate
uint32_t rol(uint32_t x, int k) {
    return ((x << k) & 0xFFFFFFFFu) | (x >> (32 - k));
}

// Right rotate
uint32_t ror(uint32_t x, int k) {
    return ((x >> k) & 0xFFFFFFFFu) | (x << (32 - k));
}

// L transform enc/dec
uint32_t SM4_L(uint32_t x) {
    return x ^ rol(x, 2) ^ rol(x, 10) ^ rol(x, 18) ^ rol(x, 24);
}

// L' transform key exp
uint32_t SM4_LP(uint32_t x) {
    // XOR with rotations 13,23 to mix key bits differently
    return x ^ rol(x, 13) ^ rol(x, 23);
}

// LT transform analysis
uint32_t SM4_LT(uint32_t x) {
    return x ^ ror(x, 2) ^ ror(x, 10) ^ ror(x, 18) ^ ror(x, 24);
}

// Tau S-box apply
uint32_t tau(uint32_t x) {
    uint32_t res = 0;
    // Apply S-box to each byte
    for (int i = 0; i < 4; ++i) {
        uint8_t b = (x >> ((3 - i) * 8)) & 0xFF;
        res |= ((uint32_t)sbox[b] << ((3 - i) * 8));
    }
    return res;
}

// Init CK -? compute round constants from formula
void init_ck() {
    // Compute CK[i] bytes as (4*i + j)*7 % 256 for variety
    for (int i = 0; i < 32; ++i) {
        uint32_t ck = 0;
        for (int j = 0; j < 4; ++j) {
            int byte_val = ((4 * i + j) * 7) % 256;
            ck |= ((uint32_t)byte_val << ((3 - j) * 8));
        }
        CK[i] = ck;
    }
}

// Key expandgenerate round keys from master key
void key_expand(const uint32_t mk[4], uint32_t rk[32]) {
    uint32_t k[36];
    // Init k[0-3] = mk ^ FK to start schedule
    for (int i = 0; i < 4; ++i) {
        k[i] = mk[i] ^ FK[i];
    }
    // Generate rk using XOR, tau, LP to derive subkeys
    for (int i = 0; i < 32; ++i) {
        uint32_t tmp = k[i + 1] ^ k[i + 2] ^ k[i + 3] ^ CK[i];
        uint32_t t = tau(tmp);
        uint32_t lp = SM4_LP(t);
        k[i + 4] = k[i] ^ lp;
        rk[i] = k[i + 4];
    }
}

// SM4 encrypt: process 128-bit block to ciphertext
void sm4_encrypt(const uint32_t plain[4], const uint32_t rk[32], uint32_t cipher[4]) {
    uint32_t x[36];
    // Load plain into state
    for (int i = 0; i < 4; ++i) {
        x[i] = plain[i];
    }
    // 32 rounds: XOR rk, tau for confusion, L for diffusion
    for (int i = 0; i < 32; ++i) {
        uint32_t tmp = x[i + 1] ^ x[i + 2] ^ x[i + 3] ^ rk[i];
        uint32_t t = tau(tmp);
        uint32_t l = SM4_L(t);
        x[i + 4] = x[i] ^ l;
    }
    // Reverse final state for output order
    for (int i = 0; i < 4; ++i) {
        cipher[i] = x[35 - i];
    }
}

// SM4 decrypt: process 128-bit block to plaintext
void sm4_decrypt(const uint32_t cipher[4], const uint32_t rk[32], uint32_t plain[4]) {
    uint32_t x[36];
    // Load cipher into state
    for (int i = 0; i < 4; ++i) {
        x[i] = cipher[i];
    }
    // 32 rounds with reversed rk to undo encryption
    for (int i = 0; i < 32; ++i) {
        uint32_t tmp = x[i + 1] ^ x[i + 2] ^ x[i + 3] ^ rk[31 - i];
        uint32_t t = tau(tmp);
        uint32_t l = SM4_L(t);
        x[i + 4] = x[i] ^ l;
    }
    // Reverse final state for output order
    for (int i = 0; i < 4; ++i) {
        plain[i] = x[35 - i];
    }
}

// Precompute LAT: build table for S-box biases
void precompute_lat() {
    // Zero LAT array
    memset(lat, 0, sizeof(lat));
    // Count parity agreements for each mask pair
    for (int x = 0; x < 256; ++x) {
        uint8_t sx = sbox[x];
        for (int a = 0; a < 256; ++a) {
            int parity_xa = __builtin_parity((uint32_t)(x & a));
            for (int b = 0; b < 256; ++b) {
                int parity_sxb = __builtin_parity((uint32_t)(sx & b));
                if (parity_xa == parity_sxb) {
                    ++lat[a][b];
                }
            }
        }
    }
}

// Multi-round bias-> approx using piling-up lemma
double compute_multi_round_bias(double single_bias) {
    const int num_rounds = 6;
    // Piling-up: 2^{n-1} * bias^n for combined prob
    return pow(2.0, num_rounds - 1) * pow(single_bias, num_rounds);
}

int main() {
    // Init CK constants
    init_ck();
    // Precompute LAT for analysis
    precompute_lat();

    // Test encrypt with sample data
    uint32_t key[4] = {0x01234567, 0x89abcdef, 0xfedcba98, 0x76543210};
    uint32_t rk[32];
    key_expand(key, rk);
    uint32_t plain[4] = {0x01234567, 0x89abcdef, 0xfedcba98, 0x76543210};
    uint32_t cipher[4];
    sm4_encrypt(plain, rk, cipher);
    printf("Encryption test: %08x %08x %08x %08x\n", cipher[0], cipher[1], cipher[2], cipher[3]);

    uint32_t dec_plain[4];
    sm4_decrypt(cipher, rk, dec_plain);
    printf("Decryption test: %08x %08x %08x %08x\n", dec_plain[0], dec_plain[1], dec_plain[2], dec_plain[3]);

    // Search max bias for 3 active S-boxes in single round
    double max_bias = 0.0;
    uint32_t best_gamma = 0;
    uint64_t active_three = 0;
    const uint64_t total = 1ULL << 32;
    uint64_t progress_counter = 0;

    // Parallel search over all gamma masks
    #pragma omp parallel reduction(+:active_three, progress_counter)
    {
        double local_max_bias = 0.0;
        uint32_t local_best_gamma = 0;
        #pragma omp for schedule(dynamic)
        for (uint64_t gg = 0; gg < total; ++gg) {
            uint32_t gamma = static_cast<uint32_t>(gg);
            // Compute delta as LT(gamma) for output mask
            uint32_t delta = SM4_LT(gamma);
            uint8_t alpha[4], beta[4];
            int num_active = 0;
            double prod_eps = 1.0;
            bool valid = true;
            // Split masks to byte level for S-box
            for (int i = 0; i < 4; ++i) {
                alpha[i] = (gamma >> ((3 - i) * 8)) & 0xFF;
                beta[i] = (delta >> ((3 - i) * 8)) & 0xFF;
                if (alpha[i] == 0 && beta[i] == 0) {
                    continue;
                }
                if (alpha[i] == 0 || beta[i] == 0) {
                    valid = false;
                    break;
                }
                // Get bias eps from LAT entry
                int entry = lat[alpha[i]][beta[i]];
                double eps_i = fabs(static_cast<double>(entry) / 256.0 - 0.5);
                if (eps_i == 0.0) {
                    valid = false;
                    break;
                }
                prod_eps *= eps_i;
                ++num_active;
            }
            // If exactly 3 active and valid, compute round bias
            if (valid && num_active == 3) {
                ++active_three;
                double bias = pow(2.0, static_cast<double>(num_active) - 1.0) * prod_eps;
                if (bias > local_max_bias) {
                    local_max_bias = bias;
                    local_best_gamma = gamma;
                }
            }

            // Update progress counter
            if ((gg % (1ULL << 20)) == 0) {
                progress_counter += (1ULL << 20);
            }
        }
        // Merge local max to global
        #pragma omp critical
        {
            if (local_max_bias > max_bias) {
                max_bias = local_max_bias;
                best_gamma = local_best_gamma;
            }
        }
    }
    uint64_t total_checked = total;
    double log_bias = log2(max_bias);
    double perc = static_cast<double>(progress_counter) / static_cast<double>(total) * 100.0;
    printf("Progress: %llu / %llu (%.2f%%), Active 3: %llu\n", progress_counter, total, perc, active_three);
    printf("Best gamma: 0x%08x\n", best_gamma);
    printf("Max single-round bias: %e (2^{%.4f})\n", max_bias, log_bias);
    printf("Total checked: %llu, with 3 active: %llu\n", total_checked, active_three);

    // Print biases for selected rounds
    const int approximated_rounds[] = {5, 6, 10, 11, 15, 16};
    set<int> approx_set(approximated_rounds, approximated_rounds + 6);
    for (int r = 1; r <= 18; ++r) {
        double round_bias = approx_set.count(r) ? max_bias : 0.0;
        double round_log_bias = approx_set.count(r) ? log_bias : -INFINITY;
        printf("Bias for round %d: %e (2^{%.4f})\n", r, round_bias, round_log_bias);
    }

    // Compute multi-round bias over 6 rounds
    double multi_bias = compute_multi_round_bias(max_bias);
    double multi_log_bias = log2(multi_bias);
    printf("Final multi round bias (piling up over 6 rounds): %e (2^{%.4f})\n", multi_bias, multi_log_bias);

    // Verify L(LT(gamma)) == gamma for consistency
    uint32_t delta = SM4_LT(best_gamma);
    uint32_t check = SM4_L(delta);
    if (check == best_gamma) {
        printf("Verification: L(L^T(gamma)) == gamma? OK if masks consistent.\n");
    }

    return 0;
}