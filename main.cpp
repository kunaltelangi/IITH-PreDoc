#include <bits/stdc++.h>
using namespace std;
using u32 = uint32_t;
using u64 = uint64_t;

int desired_active_count = 3;            // -1 = don't care
int requirements[4] = {0, 1, 1, 1};      // per-byte: -1 don't care, 0 inactive, 1 active
bool require_both_nonzero = true;        // active iff (alpha!=0 && beta!=0) when true

const int approximated_rounds[] = {5, 6, 10, 11, 15, 16};
const int num_approximated_rounds = sizeof(approximated_rounds) / sizeof(approximated_rounds[0]);
const int START = 2;
const int END   = 19;

// If you want to filter gammas during search by per-approximated-round parity bits,
// set require_mask_on_approximated = true and fill required_mask_bits[] with 0/1 or -1 (don't care).
// By default we DON'T filter (search full gamma space) and then apply the best_gamma afterwards.
const bool require_mask_on_approximated = false;
int required_mask_bits[num_approximated_rounds] = { -1, -1, -1, -1, -1, -1 };

// sm4 sbox
static const uint8_t sbox[256] = {
    0xd6,0x90,0xe9,0xfe,0xcc,0xe1,0x3d,0xb7,0x16,0xb6,0x14,0xc2,0x28,0xfb,0x2c,0x05,
    0x2b,0x67,0x9a,0x76,0x2a,0xbe,0x04,0xc3,0xaa,0x44,0x13,0x26,0x49,0x86,0x06,0x99,
    0x9c,0x42,0x50,0xf4,0x91,0xef,0x98,0x7a,0x33,0x54,0x0b,0x43,0xed,0xcf,0xac,0x62,
    0xe4,0xb3,0x1c,0xa9,0xc9,0x08,0xe8,0x95,0x80,0xdf,0x94,0xfa,0x75,0x8f,0x3f,0xa6,
    0x47,0x07,0xa7,0xfc,0xf3,0x73,0x17,0xba,0x83,0x59,0x3c,0x19,0xe6,0x85,0x4f,0xa8,
    0x68,0x6b,0x81,0xb2,0x71,0x64,0xda,0x8b,0xf8,0xeb,0x0f,0x4b,0x70,0x56,0x9d,0x35,
    0x1e,0x24,0x0e,0x5e,0x63,0x58,0xd1,0xa2,0x25,0x22,0x7c,0x3b,0x01,0x21,0x78,0x87,
    0xd4,0x00,0x46,0x57,0x9f,0xd3,0x27,0x52,0x4c,0x36,0x02,0xe7,0xa0,0xc4,0xc8,0x9e,
    0xea,0xbf,0x8a,0xd2,0x40,0xc7,0x38,0xb5,0xa3,0xf7,0xf2,0xce,0xf9,0x61,0x15,0xa1,
    0xe0,0xae,0x5d,0xa4,0x9b,0x34,0x1a,0x55,0xad,0x93,0x32,0x30,0xf5,0x8c,0xb1,0xe3,
    0x1d,0xf6,0xe2,0x2e,0x82,0x66,0xca,0x60,0xc0,0x29,0x23,0xab,0x0d,0x53,0x4e,0x6f,
    0xd5,0xdb,0x37,0x45,0xde,0xfd,0x8e,0x2f,0x03,0xff,0x6a,0x72,0x6d,0x6c,0x5b,0x51,
    0x8d,0x1b,0xaf,0x92,0xbb,0xdd,0xbc,0x7f,0x11,0xd9,0x5c,0x41,0x1f,0x10,0x5a,0xd8,
    0x0a,0xc1,0x31,0x88,0xa5,0xcd,0x7b,0xbd,0x2d,0x74,0xd0,0x12,0xb8,0xe5,0xb4,0xb0,
    0x89,0x69,0x97,0x4a,0x0c,0x96,0x77,0x7e,0x65,0xb9,0xf1,0x09,0xc5,0x6e,0xc6,0x84,
    0x18,0xf0,0x7d,0xec,0x3a,0xdc,0x4d,0x20,0x79,0xee,0x5f,0x3e,0xd7,0xcb,0x39,0x48
};

struct RoundInfo {
    int round;
    u32 F_input;
    u32 P_out;
    bool is_approximated;
    u32 gamma;
    u32 delta;
    int dot_on_tmp;
    int dot_on_pout;
};

//LAT Table
static int lat_tbl[256][256];

// Rotate right
static inline u32 ror32(u32 x, int k) {
    return (x >> k) | (x << (32 - k));
}

// Rotate left
static inline u32 rotl(u32 x, int k) { 
    return (u32)((x << k) | (x >> (32 - k))); 
}

static inline u32 L_key(u32 b) { 
    return b ^ rotl(b, 13) ^ rotl(b, 23); 
}

static inline u32 L_cipher(u32 b) { 
    return b ^ rotl(b, 2) ^ rotl(b, 10) ^ rotl(b, 18) ^ rotl(b, 24); 
}

// Non-linear transformation (tau)
static inline u32 tau(u32 x) {
    u32 r = 0;
    for (int i = 0; i < 4; ++i) {
        uint8_t byte = (x >> (8 * i)) & 0xFF;
        r |= ((u32)sbox[byte]) << (8 * i);
    }
    return r;
}

// Fixed constants used in SM4 key schedule
const u32 FK[4] = {0xa3b1bac6u, 0x56aa3350u, 0x677d9197u, 0xb27022dcu};

// Compute round constant CK[i] for SM4 key schedule
static inline u32 compute_ck(int i) {
    u32 ck = 0;
    for (int j = 0; j < 4; ++j) {
        uint8_t b = ((4 * i + j) * 7) & 0xff;
        ck |= ((u32)b) << (8 * j);
    }
    return ck;
}

// SM4 key schedule
vector<u32> key_schedule(const u32 MK[4]) {
    u32 K[36];
    for (int i = 0; i < 4; ++i) K[i] = MK[i] ^ FK[i];
    vector<u32> rk(32);
    for (int i = 0; i < 32; ++i) {
        u32 tmp = K[i+1] ^ K[i+2] ^ K[i+3] ^ compute_ck(i);
        K[i+4] = K[i] ^ L_key(tau(tmp));
        rk[i] = K[i+4];
    }
    return rk;
}

// SM4 encryption
void sm4_encrypt(const u32 pt[4], const u32 key[4], u32 ct[4]) {
    auto rk = key_schedule(key);
    u32 X[36];
    for (int i = 0; i < 4; ++i) X[i] = pt[i];
    for (int i = 0; i < 32; ++i) {
        u32 tmp = X[i+1] ^ X[i+2] ^ X[i+3] ^ rk[i];
        X[i+4] = X[i] ^ L_cipher(tau(tmp));
    }
    ct[0]=X[35]; ct[1]=X[34]; ct[2]=X[33]; ct[3]=X[32];
}

// SM4 decryption
void sm4_decrypt(const u32 ct[4], const u32 key[4], u32 pt[4]) {
    auto rk = key_schedule(key);
    reverse(rk.begin(), rk.end()); // reverse round keys for decryption
    u32 X[36];
    for (int i = 0; i < 4; ++i) X[i] = ct[i];
    for (int i = 0; i < 32; ++i) {
        u32 tmp = X[i+1] ^ X[i+2] ^ X[i+3] ^ rk[i];
        X[i+4] = X[i] ^ L_cipher(tau(tmp));
    }
    pt[0]=X[35]; pt[1]=X[34]; pt[2]=X[33]; pt[3]=X[32];
}

// Linear transpose 
static inline u32 SM4_LT(u32 x) { 
    return x ^ ror32(x,2) ^ ror32(x,10) ^ ror32(x,18) ^ ror32(x,24); 
}

// Precompute LAT 
void precompute_lat() {
    memset(lat_tbl, 0, sizeof(lat_tbl));
    for (int x = 0; x < 256; ++x) {
        uint8_t sx = sbox[x];
        for (int a = 0; a < 256; ++a) {
            int px = __builtin_parity((unsigned)(x & a));
            for (int b = 0; b < 256; ++b) {
                int ps = __builtin_parity((unsigned)(sx & b));
                if (px == ps) ++lat_tbl[a][b];
            }
        }
    }
}

inline bool is_active_byte(uint8_t alpha, uint8_t beta) {
    if (require_both_nonzero) return (alpha != 0 && beta != 0);
    else return (alpha != 0 || beta != 0);
}

bool check_gamma_and_compute_bias(u32 gamma, int &out_num_active, double &out_prod_eps) {
    u32 delta = SM4_LT(gamma);
    out_num_active = 0;
    out_prod_eps = 1.0;
    for (int i = 0; i < 4; ++i) {
        uint8_t alpha = (gamma >> ((3 - i) * 8)) & 0xFF;
        uint8_t beta  = (delta >> ((3 - i) * 8)) & 0xFF;
        if (requirements[i] == 0) {
            if (!(alpha == 0 && beta == 0)) return false;
            continue;
        }
        if (requirements[i] == 1) {
            if (!is_active_byte(alpha, beta)) return false;
            int entry = lat_tbl[alpha][beta];
            double eps = fabs(entry / 256.0 - 0.5);
            if (eps == 0.0) return false;
            out_prod_eps *= eps;
            ++out_num_active;
            continue;
        }
        /*
        if (is_active_byte(alpha, beta)) {
            int entry = lat_tbl[alpha][beta];
            double eps = fabs(entry / 256.0 - 0.5);
            if (eps == 0.0) return false;
            out_prod_eps *= eps;
            ++out_num_active;
        }
            */
    }
    if (desired_active_count >= 0 && out_num_active != desired_active_count) return false;
    return true;
}

/*
static inline int parity32(u32 x) {
    return __builtin_parity((unsigned)x);
}
*/
string hex32(u32 v) {
    ostringstream ss;
    ss << "0x" << hex << nouppercase << setfill('0') << setw(8) << (v & 0xFFFFFFFFu);
    return ss.str();
}

void linear_trail(u32 best_gamma, const u32 pt[4], const u32 key[4]) {
    auto rk = key_schedule(key);
    vector<u32> X(36);
    for (int i = 0; i < 4; ++i) X[i] = pt[i];

    vector<RoundInfo> info(32);
    u32 delta_best = SM4_LT(best_gamma);

    for (int i = 0; i < 32; ++i) {
        u32 tmp = X[i+1] ^ X[i+2] ^ X[i+3] ^ rk[i];
        u32 pout = X[i] ^ L_cipher(tau(tmp));
        X[i+4] = pout;

        RoundInfo ri;
        ri.round = i+1;
        ri.F_input = tmp;
        ri.P_out = X[i+4];
        ri.is_approximated = false;
        ri.gamma = 0;
        ri.delta = 0;
        ri.dot_on_tmp = -1;
        ri.dot_on_pout = -1;

        for (int j = 0; j < (int)num_approximated_rounds; ++j) {
            if (ri.round == approximated_rounds[j]) {
                ri.is_approximated = true;
                ri.gamma = best_gamma;
                ri.delta = delta_best;
                ri.dot_on_tmp = ri.gamma & ri.F_input;
                ri.dot_on_pout = ri.gamma & ri.P_out;
                break;
            }
        }
        info[i] = ri;
    }

    cout << "\n=== 18-Round Linear Trail Applied (Rounds 2..19) ===\n\n";
    cout << "Selected gamma = " << hex32(best_gamma) << "    delta = " << hex32(delta_best) << "\n\n";
/*
    cout << "Approximated-round mask dot-product (Γ · tmp) decomposition:\n\n";
    for (int j = 0; j < (int)num_approximated_rounds; ++j) {
        int round_no = approximated_rounds[j];
        if (round_no < 1 || round_no > 32) continue;
        const RoundInfo &ri = info[round_no - 1];
        u32 tmp = ri.F_input;
        int dot_tmp = ri.dot_on_tmp;
        int dot_pout = ri.dot_on_pout;
        int i = round_no - 1;

        int d1 = parity32(best_gamma & X[i+1]);
        int d2 = parity32(best_gamma & X[i+2]);
        int d3 = parity32(best_gamma & X[i+3]);
        int drk = parity32(best_gamma & rk[i]);
        int check = d1 ^ d2 ^ d3 ^ drk;

        cout << "Round " << setw(2) << round_no << ":\n";
        cout << "  F-input (tmp) = " << hex32(tmp) << "\n";
        cout << "  contributions: parity(Gamma & X[" << (i+1) << "])=" << d1
             << "  parity(Gamma & X[" << (i+2) << "])=" << d2
             << "  parity(Gamma & X[" << (i+3) << "])=" << d3
             << "  parity(Gamma & RK[" << i << "])=" << drk << "\n";
        cout << "  XOR(check) = " << check << "   (should equal Gamma·tmp = " << dot_tmp << ")\n";
        cout << "  Gamma·tmp   = " << dot_tmp << "    Gamma·P_out = " << dot_pout << "\n\n";
    }
        */

    const int C1 = 6, C2 = 13, C3 = 16, C4 = 12, C5 = 12, C6 = 14;
    auto draw_line = [&]() {
        cout << '+'
             << string(C1+1, '-') << '+'
             << string(C2+1, '-') << '+'
             << string(C3+1, '-') << '+'
             << string(C4+1, '-') << '+'
             << string(C5+1, '-') << '+'
             << string(C6+1, '-') << '+'
             << '\n';
    };

    draw_line();
    cout << "| " << setw(C1) << left << "Rnd"
         << "| " << setw(C2) << left << "Approx?"
         << "| " << setw(C3) << left << "F_input"
         << "| " << setw(C4) << left << "Gamma"
         << "| " << setw(C5) << left << "Delta"
         << "| " << setw(C6) << left << "P_out"
         << "|\n";
    draw_line();
    for (int r = START; r <= END; ++r) {
        const RoundInfo &ri = info[r-1];
        cout << "| " << setw(C1) << right << ri.round;
        cout << "| " << setw(C2) << left << (ri.is_approximated ? "YES" : "NO");
        if (ri.is_approximated) {
            cout << "| " << setw(C3) << left << hex32(ri.F_input);
            cout << "| " << setw(C4) << left << hex32(ri.gamma);
            cout << "| " << setw(C5) << left << hex32(ri.delta);
        } else {
            cout << "| " << setw(C3) << left << "-"
                 << "| " << setw(C4) << left << "-"
                 << "| " << setw(C5) << left << "-";
        }
        cout << "| " << setw(C6) << left << hex32(ri.P_out) << "|\n";
    }
    draw_line();
    cout << "\nPropagation display complete.\n\n";
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    u32 key[4] = {0x01234567u, 0x89abcdefu, 0xfedcba98u, 0x76543210u};
    u32 pt[4]  = {0x01234567u, 0x89abcdefu, 0xfedcba98u, 0x76543210u};
    u32 ct[4], dec[4];

    cout << "SM4 LAT search + linear trail propagation\n";
    cout << "Configuration: desired_active_count=" << desired_active_count << ", requirements={";
    for (int i = 0; i < 4; ++i) {
        if (requirements[i] == -1) cout << 'X'; else cout << requirements[i];
        if (i < 3) cout << ' ';
    }
    cout << "}\n";
    //cout << "Mask filtering on approximated rounds: " << (require_mask_on_approximated ? "ENABLED" : "DISABLED") << "\n";
    /*
    if (require_mask_on_approximated) {
        cout << " Per-approximated-round required bits: ";
        for (int j = 0; j < (int)num_approximated_rounds; ++j) {
            cout << required_mask_bits[j] << (j+1 < (int)num_approximated_rounds ? ' ' : '\n');
        }
    }
    */
    cout << "\nPrecomputing LAT...\n";
    precompute_lat();
    cout << "LAT computed.\n\n";

 
    auto rk_pre = key_schedule(key);
    vector<u32> Xpre(36);
    for (int i = 0; i < 4; ++i) Xpre[i] = pt[i];
    vector<u32> tmp_for_round(32, 0);
    for (int i = 0; i < 32; ++i) {
        u32 tmp = Xpre[i+1] ^ Xpre[i+2] ^ Xpre[i+3] ^ rk_pre[i];
        tmp_for_round[i] = tmp;
        Xpre[i+4] = Xpre[i] ^ L_cipher(tau(tmp));
    }

    double max_bias = 0.0;
    u32 best_gamma = 0;
    u64 matched_count = 0;
    const u64 TOTAL = (1ULL << 32);
    const u64 STEP  = (1ULL << 22);

    for (u64 gg = 0; gg < TOTAL; ++gg) {
        u32 gamma = (u32)gg;
/*
        if (require_mask_on_approximated) {
            bool ok = true;
            for (int j = 0; j < (int)num_approximated_rounds; ++j) {
                if (required_mask_bits[j] == -1) continue;
                int round_no = approximated_rounds[j];
                int idx = round_no - 1;
                int bit = parity32(gamma & tmp_for_round[idx]);
                if (bit != required_mask_bits[j]) { ok = false; break; }
            }
            if (!ok) continue;
        }
*/
        int num_active = 0;
        double prod_eps = 1.0;
        if (check_gamma_and_compute_bias(gamma, num_active, prod_eps)) {
            double bias = pow(2.0, num_active - 1.0) * prod_eps;
            ++matched_count;
            if (bias > max_bias) {
                 max_bias = bias; 
                 best_gamma = gamma; 
            }
        }

        if ((gg & (STEP - 1)) == 0) {
            double pct = (double)gg / (double)TOTAL * 100.0;
            cout << "Progress: " << dec << (uint64_t)gg << "/" << TOTAL << " (" << fixed << setprecision(2) << pct << "%) matches: " << matched_count << "\r";
            cout.flush();
        }
    }

    cout << "\nSearch finished. Matches: " << matched_count << "\n";
    if (matched_count == 0) return 0;

    cout << "Best gamma: " << hex32(best_gamma) << "\n";
    u32 delta_best = SM4_LT(best_gamma);
    cout << "alpha bytes: ";
    for (int i = 0; i < 4; ++i)
        cout << hex << setw(2) << setfill('0') << ((best_gamma >> ((3 - i) * 8)) & 0xFF) << ' ';
    cout << dec << setfill(' ') << "\n";
    cout << "beta  bytes: ";
    for (int i = 0; i < 4; ++i)
        cout << hex << setw(2) << setfill('0') << ((delta_best >> ((3 - i) * 8)) & 0xFF) << ' ';
    cout << dec << setfill(' ') << "\n\n";

    cout << "Max bias for single round = " << scientific << max_bias << " (2^" << fixed << setprecision(4) << log2(max_bias) << ")\n";
    double all_round_bias = pow(2.0, 5.0) * pow(max_bias, 6.0);
    cout << "18-round bias using piling-up lemma: " << scientific << all_round_bias << " (2^" << fixed << setprecision(4) << log2(all_round_bias) << ")\n\n";

    linear_trail(best_gamma, pt, key);

    return 0;
}
