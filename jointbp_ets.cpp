#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <deque>
#include <limits>
#include <system_error>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#if defined(_WIN32)
#include <io.h>
#else
#include <unistd.h>
#endif

namespace {

struct APM {
    long long a = -1;
    long long b = -1;
};

struct Params {
    long long P = 0;
    int J = 0;
    int L = 0;
    int L2 = 0;
    std::vector<long long> a;
    std::vector<long long> b;
    std::vector<long long> c;
    std::vector<long long> d;
};

using Msg = std::array<double, 4>;

static long long mod_norm(long long x, long long m) {
    long long r = x % m;
    if (r < 0) r += m;
    return r;
}

static long long egcd(long long a, long long b, long long &x, long long &y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    long long x1 = 0, y1 = 0;
    long long g = egcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

static std::optional<long long> modinv(long long a, long long m) {
    a = mod_norm(a, m);
    if (m == 1) return 0;
    long long x = 0, y = 0;
    long long g = egcd(a, m, x, y);
    if (g != 1) return std::nullopt;
    return mod_norm(x, m);
}

using BlockMatrix = std::vector<std::vector<std::optional<APM>>>;

static int mod_index(int x, int m) {
    int r = x % m;
    if (r < 0) r += m;
    return r;
}

static BlockMatrix build_block_matrix_hx(
    long long P,
    int J,
    int L2,
    const std::vector<long long> &a,
    const std::vector<long long> &b,
    const std::vector<long long> &c,
    const std::vector<long long> &d
) {
    int L = 2 * L2;
    BlockMatrix blocks(J, std::vector<std::optional<APM>>(L));
    for (int r = 0; r < J; ++r) {
        for (int cc = 0; cc < L; ++cc) {
            if (cc < L2) {
                int idx = mod_index(cc - r, L2);
                blocks[r][cc] = APM{mod_norm(a[idx], P), mod_norm(b[idx], P)};
            } else {
                int t = cc - L2;
                int idx = mod_index(t - r, L2);
                blocks[r][cc] = APM{mod_norm(c[idx], P), mod_norm(d[idx], P)};
            }
        }
    }
    return blocks;
}

static BlockMatrix build_block_matrix_hz(
    long long P,
    int J,
    int L2,
    const std::vector<long long> &a,
    const std::vector<long long> &b,
    const std::vector<long long> &c,
    const std::vector<long long> &d
) {
    int L = 2 * L2;
    BlockMatrix blocks(J, std::vector<std::optional<APM>>(L));
    for (int r = 0; r < J; ++r) {
        for (int cc = 0; cc < L; ++cc) {
            if (cc < L2) {
                int k = mod_index(cc - r, L2);
                int idx = mod_index(-k, L2);
                auto invc = modinv(mod_norm(c[idx], P), P);
                if (!invc) continue;
                APM blk;
                blk.a = *invc;
                blk.b = mod_norm(-(*invc) * mod_norm(d[idx], P), P);
                blocks[r][cc] = blk;
            } else {
                int t = cc - L2;
                int k = mod_index(t - r, L2);
                int idx = mod_index(-k, L2);
                auto inva = modinv(mod_norm(a[idx], P), P);
                if (!inva) continue;
                APM blk;
                blk.a = *inva;
                blk.b = mod_norm(-(*inva) * mod_norm(b[idx], P), P);
                blocks[r][cc] = blk;
            }
        }
    }
    return blocks;
}

static bool read_params(const std::string &path, Params &out) {
    std::ifstream in(path);
    if (!in) return false;
    unsigned long long seed = 0;
    if (!(in >> out.P >> out.J >> out.L)) return false;
    if (!(in >> seed)) return false;
    if (out.L % 2 != 0) return false;
    out.L2 = out.L / 2;

    out.a.resize(out.L2);
    out.b.resize(out.L2);
    out.c.resize(out.L2);
    out.d.resize(out.L2);
    for (int i = 0; i < out.L2; ++i) if (!(in >> out.a[i])) return false;
    for (int i = 0; i < out.L2; ++i) if (!(in >> out.b[i])) return false;
    for (int i = 0; i < out.L2; ++i) if (!(in >> out.c[i])) return false;
    for (int i = 0; i < out.L2; ++i) if (!(in >> out.d[i])) return false;
    return true;
}

static std::string format_affine_pairs(
    const std::vector<long long> &a,
    const std::vector<long long> &b
) {
    std::ostringstream oss;
    oss << "[";
    size_t n = std::min(a.size(), b.size());
    for (size_t i = 0; i < n; ++i) {
        if (i > 0) oss << " ";
        oss << "(" << a[i] << "," << b[i] << ")";
    }
    oss << "]";
    return oss.str();
}

static std::string format_double_fixed(double value, int precision) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}

static std::string format_scaled_unit(double value, double scale, const std::string &suffix) {
    double scaled = value / scale;
    double abs_scaled = std::fabs(scaled);
    int precision = 0;
    if (abs_scaled < 10.0) {
        precision = 2;
    } else if (abs_scaled < 100.0) {
        precision = 1;
    }
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << scaled;
    return oss.str() + suffix;
}

static std::string format_qbps(double qbps) {
    double abs_qbps = std::fabs(qbps);
    if (abs_qbps >= 1e9) {
        return format_scaled_unit(qbps, 1e9, "G");
    }
    if (abs_qbps >= 1e6) {
        return format_scaled_unit(qbps, 1e6, "M");
    }
    if (abs_qbps >= 1e3) {
        return format_scaled_unit(qbps, 1e3, "k");
    }
    return format_scaled_unit(qbps, 1.0, "");
}

static std::string format_latency(double seconds) {
    double abs_sec = std::fabs(seconds);
    if (abs_sec < 1e-3) {
        return format_scaled_unit(seconds, 1e-6, "us");
    }
    return format_scaled_unit(seconds, 1e-3, "ms");
}

static std::string format_girth(int girth) {
    if (girth <= 0) return "inf";
    return std::to_string(girth);
}

static std::string format_deg_triplet(int min_v, double avg_v, int max_v) {
    return std::to_string(min_v) + "/" + format_double_fixed(avg_v, 2) + "/" +
           std::to_string(max_v);
}

static const char *kColorReset = "\033[0m";
static const char *kColorBold = "\033[1m";
static const char *kColorDim = "\033[2m";
static const char *kColorCyan = "\033[38;5;33m";
static const char *kColorMagenta = "\033[38;5;135m";
static const char *kColorGreen = "\033[38;5;34m";
static const char *kColorYellow = "\033[38;5;214m";
static const char *kColorSky = "\033[38;5;33m";
static const char *kColorLavender = "\033[38;5;69m";
static const char *kColorPurple = "\033[38;5;99m";
static const char *kColorPink = "\033[38;5;168m";
static const char *kColorPeach = "\033[38;5;209m";
static constexpr int kBoxWidth = 72;

static bool env_no_color() {
    return std::getenv("NO_COLOR") != nullptr;
}

static bool is_tty_stdout() {
#if defined(_WIN32)
    return _isatty(_fileno(stdout));
#else
    return isatty(fileno(stdout));
#endif
}

static bool is_tty_stderr() {
#if defined(_WIN32)
    return _isatty(_fileno(stderr));
#else
    return isatty(fileno(stderr));
#endif
}

static bool use_color_stdout() {
    if (env_no_color()) return false;
    return is_tty_stdout();
}

static bool use_color_stderr() {
    if (env_no_color()) return false;
    return is_tty_stderr();
}

static bool use_color_for_stream(std::ostream &os) {
    if (&os == &std::cout) return use_color_stdout();
    if (&os == &std::cerr) return use_color_stderr();
    return false;
}

static void set_color(std::ostream &os, bool use_color, const char *code) {
    if (use_color) os << code;
}

static bool g_inline_progress_enabled = false;
static bool g_progress_inline_active = false;
static size_t g_progress_block_lines = 0;

static void finish_progress_inline() {
    if (!g_progress_inline_active) return;
    std::cout << "\n";
    std::cout.flush();
    g_progress_inline_active = false;
    g_progress_block_lines = 0;
}

static void print_gradient_line(
    std::ostream &os,
    const std::string &line,
    const std::vector<const char *> &colors,
    bool use_color
) {
    if (!use_color || colors.empty()) {
        os << line << "\n";
        return;
    }
    int n = static_cast<int>(line.size());
    int segments = static_cast<int>(colors.size());
    os << kColorBold;
    for (int s = 0; s < segments; ++s) {
        int start = (n * s) / segments;
        int end = (n * (s + 1)) / segments;
        if (start >= end) continue;
        os << colors[s] << line.substr(static_cast<size_t>(start),
                                       static_cast<size_t>(end - start));
    }
    os << kColorReset << "\n";
}

static void print_ascii_banner(std::ostream &os) {
    bool use_color = use_color_for_stream(os);
    const std::vector<const char *> gradient = {
        kColorSky, kColorLavender, kColorPurple, kColorPink, kColorPeach
    };
    print_gradient_line(os, "JJJJJ  OOO  III  N   N TTTTT BBBB  PPPP      EEEE TTTTT SSSS",
                        gradient, use_color);
    print_gradient_line(os, "  J   O   O  I   NN  N   T   B   B P   P     E      T  S",
                        gradient, use_color);
    print_gradient_line(os, "  J   O   O  I   N N N   T   BBBB  PPPP      EEEE   T   SSS",
                        gradient, use_color);
    print_gradient_line(os, "J J   O   O  I   N  NN   T   B   B P         E      T     S",
                        gradient, use_color);
    print_gradient_line(os, "JJJ   OOO  III  N   N   T   BBBB  P         EEEE   T  SSSS",
                        gradient, use_color);
    os << "\n";
}

static void print_box_header(std::ostream &os, const std::string &title) {
    bool use_color = use_color_for_stream(os);
    std::string line(static_cast<size_t>(kBoxWidth), '-');
    std::string label = " " + title;
    if (static_cast<int>(label.size()) > kBoxWidth) {
        label.resize(static_cast<size_t>(kBoxWidth));
    }
    set_color(os, use_color, kColorSky);
    os << "+" << line << "+\n";
    set_color(os, use_color, kColorSky);
    os << "|";
    set_color(os, use_color, kColorBold);
    set_color(os, use_color, kColorMagenta);
    os << label;
    set_color(os, use_color, kColorSky);
    os << std::string(static_cast<size_t>(kBoxWidth - label.size()), ' ') << "|\n";
    set_color(os, use_color, kColorSky);
    os << "+" << line << "+\n";
    set_color(os, use_color, kColorReset);
}

static void print_tips(std::ostream &os, const std::vector<std::string> &tips) {
    if (tips.empty()) return;
    bool use_color = use_color_for_stream(os);
    set_color(os, use_color, kColorBold);
    set_color(os, use_color, kColorMagenta);
    os << "Tips for getting started:";
    set_color(os, use_color, kColorReset);
    os << "\n";
    for (size_t i = 0; i < tips.size(); ++i) {
        set_color(os, use_color, kColorSky);
        os << (i + 1) << ".";
        set_color(os, use_color, kColorReset);
        os << " " << tips[i] << "\n";
    }
    os << "\n";
}

static void print_kv_compact(
    std::ostream &os,
    const std::vector<std::pair<std::string, std::string>> &items,
    size_t max_width
) {
    std::string line;
    for (const auto &kv : items) {
        std::string token = kv.first + "=" + kv.second;
        if (line.empty()) {
            line = token;
            continue;
        }
        if (line.size() + 2 + token.size() <= max_width) {
            line += "  " + token;
        } else {
            os << line << "\n";
            line = token;
        }
    }
    if (!line.empty()) {
        os << line << "\n";
    }
}

static void print_kv_compact_prefixed(
    std::ostream &os,
    const std::string &prefix,
    const std::vector<std::pair<std::string, std::string>> &items,
    size_t max_width
) {
    std::string line;
    std::string indent(prefix.size(), ' ');
    for (const auto &kv : items) {
        std::string token = kv.first + "=" + kv.second;
        if (line.empty()) {
            line = prefix + token;
            continue;
        }
        if (line.size() + 2 + token.size() <= max_width) {
            line += "  " + token;
        } else {
            os << line << "\n";
            line = indent + token;
        }
    }
    if (!line.empty()) {
        os << line << "\n";
    }
}

static std::vector<std::string> build_kv_lines(
    const std::vector<std::pair<std::string, std::string>> &items,
    size_t max_width
) {
    std::vector<std::string> lines;
    std::string line;
    for (const auto &kv : items) {
        std::string token = kv.first + "=" + kv.second;
        if (line.empty()) {
            line = token;
            continue;
        }
        if (line.size() + 2 + token.size() <= max_width) {
            line += "  " + token;
        } else {
            lines.push_back(line);
            line = token;
        }
    }
    if (!line.empty()) {
        lines.push_back(line);
    }
    return lines;
}

static std::string format_kv_cell(
    const std::string &key,
    const std::string &value,
    int key_width,
    int value_width
) {
    std::ostringstream oss;
    std::string key_out = key;
    if (static_cast<int>(key_out.size()) > key_width) {
        key_out.resize(static_cast<size_t>(key_width));
    }
    oss << std::left << std::setw(key_width) << key_out
        << "="
        << std::right << std::setw(value_width) << value;
    return oss.str();
}

static std::vector<std::string> build_kv_columns(
    const std::vector<std::pair<std::string, std::string>> &items,
    int columns,
    int key_width,
    int value_width,
    int gap
) {
    std::vector<std::string> lines;
    std::string line;
    int col = 0;
    std::string gap_str(static_cast<size_t>(gap), ' ');
    for (const auto &kv : items) {
        std::string cell = format_kv_cell(kv.first, kv.second, key_width, value_width);
        if (col > 0) {
            line += gap_str;
        }
        line += cell;
        col++;
        if (col == columns) {
            lines.push_back(line);
            line.clear();
            col = 0;
        }
    }
    if (!line.empty()) {
        lines.push_back(line);
    }
    return lines;
}

static std::vector<std::string> wrap_text_line(const std::string &line, size_t width) {
    std::vector<std::string> out;
    std::string rest = line;
    while (rest.size() > width) {
        size_t cut = rest.rfind(' ', width);
        if (cut == std::string::npos || cut == 0) {
            cut = width;
        }
        out.push_back(rest.substr(0, cut));
        size_t next = cut;
        while (next < rest.size() && rest[next] == ' ') ++next;
        rest = rest.substr(next);
    }
    if (!rest.empty()) {
        out.push_back(rest);
    }
    return out;
}

static size_t print_boxed_block(
    std::ostream &os,
    const std::string &title,
    const std::vector<std::string> &lines
) {
    bool use_color = use_color_for_stream(os);
    std::string border(static_cast<size_t>(kBoxWidth), '-');
    std::string label = " " + title;
    if (label.size() > static_cast<size_t>(kBoxWidth)) {
        label.resize(static_cast<size_t>(kBoxWidth));
    }
    size_t printed = 0;
    set_color(os, use_color, kColorSky);
    os << "+" << border << "+\n";
    ++printed;
    set_color(os, use_color, kColorSky);
    os << "|";
    set_color(os, use_color, kColorBold);
    set_color(os, use_color, kColorMagenta);
    os << label;
    set_color(os, use_color, kColorSky);
    os << std::string(static_cast<size_t>(kBoxWidth) - label.size(), ' ') << "|\n";
    ++printed;
    set_color(os, use_color, kColorSky);
    os << "+" << border << "+\n";
    ++printed;
    set_color(os, use_color, kColorReset);

    for (const auto &raw : lines) {
        auto wrapped = wrap_text_line(raw, kBoxWidth);
        for (const auto &line : wrapped) {
            os << "|" << line;
            if (line.size() < static_cast<size_t>(kBoxWidth)) {
                os << std::string(static_cast<size_t>(kBoxWidth) - line.size(), ' ');
            }
            os << "|\n";
            ++printed;
        }
    }

    set_color(os, use_color, kColorSky);
    os << "+" << border << "+\n";
    set_color(os, use_color, kColorReset);
    ++printed;
    return printed;
}

static void print_progress_block(const std::vector<std::string> &lines) {
    if (!g_inline_progress_enabled) {
        print_boxed_block(std::cout, "PROGRESS", lines);
        std::cout.flush();
        return;
    }
    if (g_progress_inline_active && g_progress_block_lines > 0) {
        std::cout << "\033[" << g_progress_block_lines << "A\r\033[J";
    }
    size_t printed = print_boxed_block(std::cout, "PROGRESS", lines);
    std::cout.flush();
    g_progress_inline_active = true;
    g_progress_block_lines = printed;
}

static void print_stdout_header(const std::string &title) {
    bool use_color = use_color_for_stream(std::cout);
    finish_progress_inline();
    print_ascii_banner(std::cout);
    set_color(std::cout, use_color, kColorBold);
    set_color(std::cout, use_color, kColorMagenta);
    std::cout << "AI Report:";
    set_color(std::cout, use_color, kColorReset);
    std::cout << " ";
    set_color(std::cout, use_color, kColorBold);
    set_color(std::cout, use_color, kColorCyan);
    std::cout << title;
    set_color(std::cout, use_color, kColorReset);
    std::cout << "\n";
}

static void print_stdout_section(const std::string &title) {
    std::cout << "\n";
    finish_progress_inline();
    print_box_header(std::cout, title);
}

static void print_help_section(const std::string &title) {
    std::cerr << "\n";
    print_box_header(std::cerr, title);
}

static std::vector<std::vector<int>> build_check_neighbors(
    const BlockMatrix &blocks,
    long long P
) {
    int J = static_cast<int>(blocks.size());
    int L = static_cast<int>(blocks[0].size());
    int m = J * static_cast<int>(P);
    std::vector<std::vector<int>> checks(m, std::vector<int>(L, -1));
    for (int r = 0; r < J; ++r) {
        for (int x = 0; x < P; ++x) {
            int check = r * static_cast<int>(P) + x;
            for (int c = 0; c < L; ++c) {
                const auto &blk = blocks[r][c];
                if (!blk) continue;
                long long y = mod_norm(blk->a * x + blk->b, P);
                checks[check][c] = c * static_cast<int>(P) + static_cast<int>(y);
            }
        }
    }
    return checks;
}

static void normalize_msg(Msg &m) {
    double sum = m[0] + m[1] + m[2] + m[3];
    if (sum <= 0.0) {
        m = {0.25, 0.25, 0.25, 0.25};
        return;
    }
    for (double &v : m) v /= sum;
}

static int xbit(int state) {
    return (state == 1 || state == 3) ? 1 : 0;
}

static int zbit(int state) {
    return (state == 2 || state == 3) ? 1 : 0;
}

static std::vector<int> diff_indices_bits(
    const std::vector<int> &truth,
    const std::vector<int> &est,
    bool use_xbit
);
static int unsatisfied_count(
    const std::vector<int> &hat,
    const std::vector<int> &target
);
static std::vector<int> unsatisfied_indices(
    const std::vector<int> &hat,
    const std::vector<int> &target
);
static void print_index_list(const std::string &label, const std::vector<int> &idxs);
static void write_index_list(std::ostream &os, const std::string &label, const std::vector<int> &idxs);
static std::vector<int> sorted_unique(std::vector<int> v);
static std::vector<std::string> sorted_unique_strings(std::vector<std::string> v);

struct ETSIndex;
struct EdgeRef;

struct FlipPPBest {
    bool has = false;
    int weight = std::numeric_limits<int>::max();
    std::vector<int> vars;
    std::vector<int> delta;
    std::string tag;
};

struct FlipSetRef {
    std::string tag;
    const std::vector<int> *vars = nullptr;
};

struct ETSFixResult {
    bool applied = false;
    int size = 0;
    std::string label;
};

static std::string format_ets_label(const ETSFixResult &res) {
    if (!res.applied) return "none";
    if (!res.label.empty()) return res.label;
    if (res.size > 0) return std::to_string(res.size);
    return "none";
}

struct CycleEntry {
    std::array<int, 5> vars{};
    std::array<int, 5> checks{};
};

struct CycleList {
    std::vector<CycleEntry> cycles8;
    std::vector<CycleEntry> cycles10;
    bool have8 = false;
    bool have10 = false;
};

struct CycleIndex {
    CycleList x;
    CycleList z;
    bool any_loaded = false;
};

struct OSDPPInfo {
    bool computed = false;
    bool solvable = false;
    bool unique = false;
    int rank = 0;
    int rows = 0;
    int cols = 0;
    std::vector<int> delta;
};

struct OSDMinWeightCandidate {
    bool has = false;
    int k = 0;
    int dof = 0;
    int rows = 0;
    int cols = 0;
    int weight = 0;
    std::vector<uint64_t> delta_bits;
};

struct ETSIndexRef {
    int size = 0;
    std::string label;
    const ETSIndex *x = nullptr;
    const ETSIndex *z = nullptr;
};

struct PPEarlyContext {
    bool enable = false;
    bool verbose = false;
    std::vector<ETSIndexRef> ets;
    const CycleIndex *cycles = nullptr;
};

static std::vector<FlipSetRef> build_flip_sets(
    const std::vector<int> &union_vars,
    const std::vector<std::vector<int>> &window,
    const std::vector<int> &last_vars
);
static bool attempt_ets_fix_side(
    const std::string &label,
    const ETSIndex &ets_index,
    const std::vector<std::vector<int>> &checks,
    const std::vector<std::vector<EdgeRef>> &var_to_checks,
    const std::vector<int> &target_syndrome,
    bool use_xbit,
    std::vector<int> &est,
    bool verbose
);
static ETSFixResult attempt_ets_fix_side_multi(
    const std::string &label,
    const std::vector<ETSIndexRef> &ets_list,
    const std::vector<std::vector<int>> &checks,
    const std::vector<std::vector<EdgeRef>> &var_to_checks,
    const std::vector<int> &target_syndrome,
    bool use_xbit,
    std::vector<int> &est,
    bool verbose
);
static bool apply_flip_pp_best(
    const std::string &label,
    const FlipPPBest &best,
    const std::vector<std::vector<int>> &checks,
    const std::vector<int> &target_syndrome,
    bool use_xbit,
    std::vector<int> &est,
    bool verbose
);
static bool attempt_flip_pp_side(
    const std::string &label,
    const std::vector<int> &flip_vars_in,
    const std::vector<std::vector<int>> &checks,
    const std::vector<int> &target_syndrome,
    bool use_xbit,
    std::vector<int> &est,
    bool verbose,
    FlipPPBest *best,
    const std::string &set_tag
);
static bool attempt_osd_pp_side(
    const std::string &label,
    const std::vector<double> &abs_llr,
    const std::vector<std::vector<int>> &checks,
    const std::vector<int> &target_syndrome,
    bool use_xbit,
    std::vector<int> &est,
    bool verbose,
    const std::vector<int> *truth_err
);

struct EdgeRef {
    int check = -1;
    int pos = -1;
};

struct RowBasis {
    int nvars = 0;
    int words = 0;
    std::vector<std::vector<uint64_t>> rows;
    std::vector<int> pivot_col;
    std::vector<int> order;
};

static int highest_bit_index(const std::vector<uint64_t> &row) {
    for (int w = static_cast<int>(row.size()) - 1; w >= 0; --w) {
        uint64_t v = row[w];
        if (v == 0) continue;
        int bit = 63 - __builtin_clzll(v);
        return w * 64 + bit;
    }
    return -1;
}

static RowBasis build_row_basis(
    const std::vector<std::vector<int>> &checks,
    int nvars
) {
    RowBasis basis;
    basis.nvars = nvars;
    basis.words = (nvars + 63) / 64;
    std::vector<int> pivot_to_row(nvars, -1);
    for (const auto &chk : checks) {
        std::vector<uint64_t> row(basis.words, 0);
        for (int v : chk) {
            if (v < 0 || v >= nvars) continue;
            int w = v >> 6;
            int b = v & 63;
            row[w] |= (1ULL << b);
        }
        while (true) {
            int pivot = highest_bit_index(row);
            if (pivot < 0) break;
            int existing = pivot_to_row[pivot];
            if (existing < 0) {
                int idx = static_cast<int>(basis.rows.size());
                basis.rows.push_back(std::move(row));
                basis.pivot_col.push_back(pivot);
                pivot_to_row[pivot] = idx;
                break;
            }
            const auto &base = basis.rows[existing];
            for (int w = 0; w < basis.words; ++w) {
                row[w] ^= base[w];
            }
        }
    }
    basis.order.resize(basis.rows.size());
    for (size_t i = 0; i < basis.rows.size(); ++i) {
        basis.order[i] = static_cast<int>(i);
    }
    std::sort(basis.order.begin(), basis.order.end(),
              [&](int a, int b) { return basis.pivot_col[a] > basis.pivot_col[b]; });
    return basis;
}

static void fill_diff_bits(
    std::vector<uint64_t> &out,
    const std::vector<int> &truth,
    const std::vector<int> &est,
    bool use_xbit
) {
    std::fill(out.begin(), out.end(), 0);
    size_t n = std::min(truth.size(), est.size());
    for (size_t i = 0; i < n; ++i) {
        int t = use_xbit ? xbit(truth[i]) : zbit(truth[i]);
        int e = use_xbit ? xbit(est[i]) : zbit(est[i]);
        if ((t ^ e) != 0) {
            out[i >> 6] |= (1ULL << (i & 63));
        }
    }
}

static bool in_row_space(const RowBasis &basis, std::vector<uint64_t> &vec) {
    if (basis.words == 0) return vec.empty();
    for (int idx : basis.order) {
        int pivot = basis.pivot_col[idx];
        int w = pivot >> 6;
        uint64_t mask = 1ULL << (pivot & 63);
        if ((vec[w] & mask) == 0) continue;
        const auto &row = basis.rows[idx];
        for (int i = 0; i < basis.words; ++i) {
            vec[i] ^= row[i];
        }
    }
    for (uint64_t v : vec) {
        if (v != 0) return false;
    }
    return true;
}

static std::vector<std::vector<EdgeRef>> build_var_adjacency(
    const std::vector<std::vector<int>> &checks,
    int nvars
) {
    std::vector<std::vector<EdgeRef>> var_to_checks(nvars);
    for (int c = 0; c < static_cast<int>(checks.size()); ++c) {
        for (int pos = 0; pos < static_cast<int>(checks[c].size()); ++pos) {
            int v = checks[c][pos];
            if (v < 0 || v >= nvars) continue;
            var_to_checks[v].push_back(EdgeRef{c, pos});
        }
    }
    return var_to_checks;
}

struct DegreeStats {
    int min = 0;
    int max = 0;
    double avg = 0.0;
    long long edges = 0;
};

static DegreeStats compute_var_degree_stats(const std::vector<std::vector<EdgeRef>> &var_to_checks) {
    DegreeStats stats;
    if (var_to_checks.empty()) return stats;
    stats.min = std::numeric_limits<int>::max();
    long long sum = 0;
    for (const auto &edges : var_to_checks) {
        int deg = static_cast<int>(edges.size());
        stats.min = std::min(stats.min, deg);
        stats.max = std::max(stats.max, deg);
        sum += deg;
    }
    stats.edges = sum;
    stats.avg = static_cast<double>(sum) / static_cast<double>(var_to_checks.size());
    if (stats.min == std::numeric_limits<int>::max()) stats.min = 0;
    return stats;
}

static DegreeStats compute_check_degree_stats(const std::vector<std::vector<int>> &checks) {
    DegreeStats stats;
    if (checks.empty()) return stats;
    stats.min = std::numeric_limits<int>::max();
    long long sum = 0;
    for (const auto &chk : checks) {
        int deg = 0;
        for (int v : chk) {
            if (v >= 0) ++deg;
        }
        stats.min = std::min(stats.min, deg);
        stats.max = std::max(stats.max, deg);
        sum += deg;
    }
    stats.edges = sum;
    stats.avg = static_cast<double>(sum) / static_cast<double>(checks.size());
    if (stats.min == std::numeric_limits<int>::max()) stats.min = 0;
    return stats;
}

struct OrthCheckResult {
    bool orthogonal = true;
    long long odd_pairs = 0;
    int first_x = -1;
    int first_z = -1;
    int first_overlap = 0;
};

static OrthCheckResult check_hx_hz_orthogonality(
    const std::vector<std::vector<int>> &x_checks,
    const std::vector<std::vector<int>> &z_checks,
    int nvars
) {
    OrthCheckResult res;
    if (x_checks.empty() || z_checks.empty() || nvars <= 0) return res;
    std::vector<int> marks(nvars, -1);
    int stamp = 0;
    for (int xi = 0; xi < static_cast<int>(x_checks.size()); ++xi) {
        ++stamp;
        for (int v : x_checks[xi]) {
            if (v >= 0 && v < nvars) marks[v] = stamp;
        }
        for (int zi = 0; zi < static_cast<int>(z_checks.size()); ++zi) {
            int parity = 0;
            int overlap = 0;
            for (int v : z_checks[zi]) {
                if (v >= 0 && v < nvars && marks[v] == stamp) {
                    parity ^= 1;
                    ++overlap;
                }
            }
            if (parity != 0) {
                res.orthogonal = false;
                res.odd_pairs += 1;
                if (res.first_x < 0) {
                    res.first_x = xi;
                    res.first_z = zi;
                    res.first_overlap = overlap;
                }
            }
        }
    }
    return res;
}

static int compute_girth(
    const std::vector<std::vector<int>> &checks,
    const std::vector<std::vector<EdgeRef>> &var_to_checks
) {
    int nvars = static_cast<int>(var_to_checks.size());
    int nchecks = static_cast<int>(checks.size());
    if (nvars <= 0 || nchecks <= 0) return 0;
    int total = nvars + nchecks;
    int best = std::numeric_limits<int>::max();
    std::vector<int> dist(total, 0);
    std::vector<int> parent(total, -1);
    std::vector<int> seen(total, 0);
    std::vector<int> queue(total, 0);
    int stamp = 0;

    for (int start = 0; start < nvars; ++start) {
        if (best == 4) break;
        ++stamp;
        int qh = 0;
        int qt = 0;
        queue[qt++] = start;
        seen[start] = stamp;
        dist[start] = 0;
        parent[start] = -1;
        while (qh < qt) {
            int u = queue[qh++];
            int du = dist[u];
            if (2 * du >= best) continue;
            if (u < nvars) {
                for (const auto &edge : var_to_checks[u]) {
                    int c = edge.check;
                    if (c < 0 || c >= nchecks) continue;
                    int v = nvars + c;
                    if (seen[v] != stamp) {
                        seen[v] = stamp;
                        dist[v] = du + 1;
                        parent[v] = u;
                        queue[qt++] = v;
                    } else if (parent[u] != v) {
                        int cycle = du + dist[v] + 1;
                        if (cycle < best) best = cycle;
                    }
                }
            } else {
                int c = u - nvars;
                if (c < 0 || c >= nchecks) continue;
                for (int v_id : checks[c]) {
                    if (v_id < 0 || v_id >= nvars) continue;
                    int v = v_id;
                    if (seen[v] != stamp) {
                        seen[v] = stamp;
                        dist[v] = du + 1;
                        parent[v] = u;
                        queue[qt++] = v;
                    } else if (parent[u] != v) {
                        int cycle = du + dist[v] + 1;
                        if (cycle < best) best = cycle;
                    }
                }
            }
        }
    }
    if (best == std::numeric_limits<int>::max()) return 0;
    return best;
}

template <typename LogFn>
static void log_trapset(
    char side_label,
    const std::vector<int> &vars_in,
    const std::vector<int> &odd_checks_in,
    const std::vector<std::vector<int>> &checks,
    const std::vector<std::vector<EdgeRef>> &var_to_checks,
    const CycleList *cycles,
    LogFn log_line
) {
    if (vars_in.empty()) return;
    std::vector<int> vars = vars_in;
    std::sort(vars.begin(), vars.end());
    vars.erase(std::unique(vars.begin(), vars.end()), vars.end());
    std::vector<int> odd_checks = odd_checks_in;
    std::sort(odd_checks.begin(), odd_checks.end());
    odd_checks.erase(std::unique(odd_checks.begin(), odd_checks.end()), odd_checks.end());

    int nvars_total = static_cast<int>(var_to_checks.size());
    std::vector<char> var_in(nvars_total, 0);
    for (int v : vars) {
        if (v >= 0 && v < nvars_total) var_in[v] = 1;
    }

    std::vector<int> check_list;
    std::vector<char> check_in(checks.size(), 0);
    for (int v : vars) {
        if (v < 0 || v >= nvars_total) continue;
        for (const auto &edge : var_to_checks[v]) {
            int c = edge.check;
            if (c < 0 || c >= static_cast<int>(checks.size())) continue;
            if (!check_in[c]) {
                check_in[c] = 1;
                check_list.push_back(c);
            }
        }
    }
    std::sort(check_list.begin(), check_list.end());

    std::vector<char> odd_in(checks.size(), 0);
    for (int c : odd_checks) {
        if (c >= 0 && c < static_cast<int>(checks.size())) odd_in[c] = 1;
    }

    std::vector<int> var_pos(nvars_total, -1);
    for (size_t i = 0; i < vars.size(); ++i) {
        int v = vars[i];
        if (v >= 0 && v < nvars_total) var_pos[v] = static_cast<int>(i);
    }
    std::vector<int> check_pos(checks.size(), -1);
    for (size_t i = 0; i < check_list.size(); ++i) {
        int c = check_list[i];
        if (c >= 0 && c < static_cast<int>(checks.size())) check_pos[c] = static_cast<int>(i);
    }

    auto v_sym = [&](int v) {
        int idx = (v >= 0 && v < nvars_total) ? var_pos[v] : -1;
        if (idx < 0) return std::string("v?");
        return std::string("v") + std::to_string(idx);
    };
    auto c_sym = [&](int c) {
        int idx = (c >= 0 && c < static_cast<int>(check_pos.size())) ? check_pos[c] : -1;
        if (idx < 0) return std::string("c?");
        return std::string("c") + std::to_string(idx);
    };

    char side_lower = (side_label == 'Z') ? 'z' : 'x';
    {
        std::ostringstream oss;
        oss << "trap_" << side_lower << "_map vars_count=" << vars.size() << " vars=[";
        for (size_t i = 0; i < vars.size(); ++i) {
            if (i) oss << ",";
            oss << "v" << i << "=" << vars[i];
        }
        oss << "] checks_count=" << check_list.size() << " checks=[";
        for (size_t i = 0; i < check_list.size(); ++i) {
            if (i) oss << ",";
            oss << "c" << i << "=" << check_list[i];
        }
        oss << "]";
        log_line(oss.str());
    }

    {
        std::ostringstream oss;
        oss << "trap_" << side_lower
            << " a=" << vars.size()
            << " b=" << odd_checks.size()
            << " vars_count=" << vars.size()
            << " vars=[";
        for (size_t i = 0; i < vars.size(); ++i) {
            if (i) oss << ",";
            oss << "v" << i;
        }
        oss << "] checks_count=" << check_list.size()
            << " checks=[";
        for (size_t i = 0; i < check_list.size(); ++i) {
            if (i) oss << ",";
            oss << "c" << i;
        }
        oss << "] odd_checks_count=" << odd_checks.size()
            << " odd_checks=[";
        for (size_t i = 0; i < odd_checks.size(); ++i) {
            if (i) oss << ",";
            oss << c_sym(odd_checks[i]);
        }
        oss << "]";
        log_line(oss.str());
    }

    for (int c : check_list) {
        std::vector<int> cols;
        std::string bits(vars.size(), '0');
        for (int v : checks[c]) {
            if (v < 0 || v >= nvars_total) continue;
            int pos = var_pos[v];
            if (pos >= 0) {
                cols.push_back(v);
                bits[static_cast<size_t>(pos)] = '1';
            }
        }
        std::sort(cols.begin(), cols.end());
        std::ostringstream oss;
        oss << "trap_" << side_lower << "_row c=" << c_sym(c)
            << " odd=" << (odd_in[c] ? 1 : 0)
            << " cols_count=" << cols.size()
            << " cols=[";
        for (size_t i = 0; i < cols.size(); ++i) {
            if (i) oss << ",";
            oss << v_sym(cols[i]);
        }
        oss << "] bits=" << bits;
        log_line(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "trap_" << side_lower << "_matrix cols=[";
        for (size_t i = 0; i < vars.size(); ++i) {
            if (i) oss << " ";
            oss << "v" << i;
        }
        oss << "]";
        log_line(oss.str());
    }
    for (int c : check_list) {
        std::string bits(vars.size(), '0');
        for (int v : checks[c]) {
            if (v < 0 || v >= nvars_total) continue;
            int pos = var_pos[v];
            if (pos >= 0) bits[static_cast<size_t>(pos)] = '1';
        }
        std::ostringstream oss;
        oss << "trap_" << side_lower << "_matrix_row c=" << c_sym(c)
            << " odd=" << (odd_in[c] ? 1 : 0) << " bits=[";
        for (size_t i = 0; i < bits.size(); ++i) {
            if (i) oss << " ";
            oss << bits[i];
        }
        oss << "]";
        log_line(oss.str());
    }

    auto in_var = [&](int v) -> bool {
        return v >= 0 && v < nvars_total && var_in[v];
    };
    auto in_check = [&](int c) -> bool {
        return c >= 0 && c < static_cast<int>(check_in.size()) && check_in[c];
    };

    std::vector<CycleEntry> cycles8;
    std::vector<CycleEntry> cycles10;
    if (cycles) {
        for (const auto &entry : cycles->cycles8) {
            bool ok = true;
            for (int i = 0; i < 4; ++i) {
                if (!in_var(entry.vars[i]) || !in_check(entry.checks[i])) {
                    ok = false;
                    break;
                }
            }
            if (ok) cycles8.push_back(entry);
        }
        for (const auto &entry : cycles->cycles10) {
            bool ok = true;
            for (int i = 0; i < 5; ++i) {
                if (!in_var(entry.vars[i]) || !in_check(entry.checks[i])) {
                    ok = false;
                    break;
                }
            }
            if (ok) cycles10.push_back(entry);
        }
    }

    {
        std::ostringstream oss;
        oss << "trap_" << side_lower << "_cycle8_count=" << cycles8.size();
        log_line(oss.str());
    }
    for (size_t i = 0; i < cycles8.size(); ++i) {
        const auto &entry = cycles8[i];
        std::ostringstream oss;
        oss << "trap_" << side_lower << "_cycle8 id=" << i
            << " vars_order_count=4 vars_order=[";
        for (int j = 0; j < 4; ++j) {
            if (j) oss << ",";
            oss << v_sym(entry.vars[j]);
        }
        oss << "] checks_order_count=4 checks_order=[";
        for (int j = 0; j < 4; ++j) {
            if (j) oss << ",";
            oss << c_sym(entry.checks[j]);
        }
        oss << "] path=[";
        for (int j = 0; j < 4; ++j) {
            if (j) oss << ",";
            oss << v_sym(entry.vars[j]) << "," << c_sym(entry.checks[j]);
        }
        oss << "]";
        log_line(oss.str());
    }

    {
        std::ostringstream oss;
        oss << "trap_" << side_lower << "_cycle10_count=" << cycles10.size();
        log_line(oss.str());
    }
    for (size_t i = 0; i < cycles10.size(); ++i) {
        const auto &entry = cycles10[i];
        std::ostringstream oss;
        oss << "trap_" << side_lower << "_cycle10 id=" << i
            << " vars_order_count=5 vars_order=[";
        for (int j = 0; j < 5; ++j) {
            if (j) oss << ",";
            oss << v_sym(entry.vars[j]);
        }
        oss << "] checks_order_count=5 checks_order=[";
        for (int j = 0; j < 5; ++j) {
            if (j) oss << ",";
            oss << c_sym(entry.checks[j]);
        }
        oss << "] path=[";
        for (int j = 0; j < 5; ++j) {
            if (j) oss << ",";
            oss << v_sym(entry.vars[j]) << "," << c_sym(entry.checks[j]);
        }
        oss << "]";
        log_line(oss.str());
    }

    auto log_links = [&](const std::vector<CycleEntry> &list, int count, int length) {
        for (size_t i = 0; i < list.size(); ++i) {
            for (size_t j = i + 1; j < list.size(); ++j) {
                std::vector<int> shared_vars;
                std::vector<int> shared_checks;
                for (int a = 0; a < count; ++a) {
                    int v = list[i].vars[a];
                    for (int b = 0; b < count; ++b) {
                        if (v == list[j].vars[b]) {
                            shared_vars.push_back(v);
                            break;
                        }
                    }
                }
                for (int a = 0; a < count; ++a) {
                    int c = list[i].checks[a];
                    for (int b = 0; b < count; ++b) {
                        if (c == list[j].checks[b]) {
                            shared_checks.push_back(c);
                            break;
                        }
                    }
                }
                std::vector<int> shared_vars_local;
                std::vector<int> shared_checks_local;
                for (int v : shared_vars) {
                    int idx = (v >= 0 && v < nvars_total) ? var_pos[v] : -1;
                    if (idx >= 0) shared_vars_local.push_back(idx);
                }
                for (int c : shared_checks) {
                    int idx = (c >= 0 && c < static_cast<int>(check_pos.size())) ? check_pos[c] : -1;
                    if (idx >= 0) shared_checks_local.push_back(idx);
                }
                std::sort(shared_vars_local.begin(), shared_vars_local.end());
                shared_vars_local.erase(std::unique(shared_vars_local.begin(),
                                                    shared_vars_local.end()),
                                        shared_vars_local.end());
                std::sort(shared_checks_local.begin(), shared_checks_local.end());
                shared_checks_local.erase(std::unique(shared_checks_local.begin(),
                                                      shared_checks_local.end()),
                                          shared_checks_local.end());

                std::ostringstream oss;
                oss << "trap_" << side_lower
                    << "_link a=" << length << ":" << i
                    << " b=" << length << ":" << j
                    << " shared_vars_count=" << shared_vars_local.size()
                    << " shared_vars=[";
                for (size_t k = 0; k < shared_vars_local.size(); ++k) {
                    if (k) oss << ",";
                    oss << "v" << shared_vars_local[k];
                }
                oss << "] shared_checks_count=" << shared_checks_local.size()
                    << " shared_checks=[";
                for (size_t k = 0; k < shared_checks_local.size(); ++k) {
                    if (k) oss << ",";
                    oss << "c" << shared_checks_local[k];
                }
                oss << "]";
                log_line(oss.str());
            }
        }
    };
    if (cycles8.size() > 1) log_links(cycles8, 4, 8);
    if (cycles10.size() > 1) log_links(cycles10, 5, 10);
}

struct CycleCoverItem {
    uint64_t mask = 0;
    int length = 0;
    int id = -1;
    int weight = 0;
};

static int popcount_u64(uint64_t v) {
    return static_cast<int>(__builtin_popcountll(v));
}

static bool cycle_cover_dfs(
    const std::vector<CycleCoverItem> &items,
    const std::vector<std::vector<int>> &items_by_bit,
    uint64_t target_mask,
    int depth_left,
    uint64_t covered,
    int max_cycle_vars,
    std::vector<char> &used,
    std::vector<int> &out
) {
    if (covered == target_mask) return true;
    if (depth_left == 0) return false;
    uint64_t remaining = target_mask & ~covered;
    if (remaining == 0) return true;
    if (popcount_u64(remaining) > depth_left * max_cycle_vars) return false;
    int bit = static_cast<int>(__builtin_ctzll(remaining));
    const auto &cands = items_by_bit[bit];
    for (int idx : cands) {
        if (used[static_cast<size_t>(idx)]) continue;
        uint64_t new_covered = covered | items[idx].mask;
        if (new_covered == covered) continue;
        uint64_t next_remaining = target_mask & ~new_covered;
        if (popcount_u64(next_remaining) > (depth_left - 1) * max_cycle_vars) continue;
        used[static_cast<size_t>(idx)] = 1;
        out.push_back(idx);
        if (cycle_cover_dfs(items, items_by_bit, target_mask, depth_left - 1,
                            new_covered, max_cycle_vars, used, out)) {
            return true;
        }
        out.pop_back();
        used[static_cast<size_t>(idx)] = 0;
    }
    return false;
}

template <typename LogFn>
static void maybe_log_diff_trapset_cycle_union(
    char side_label,
    const std::vector<int> &diff_in,
    const std::vector<int> &us,
    const std::vector<std::vector<int>> &checks,
    const std::vector<std::vector<EdgeRef>> &var_to_checks,
    const CycleIndex *cycles,
    LogFn log_line
) {
    if (!cycles || diff_in.empty()) return;
    std::vector<int> vars = diff_in;
    std::sort(vars.begin(), vars.end());
    vars.erase(std::unique(vars.begin(), vars.end()), vars.end());
    if (vars.empty()) return;
    if (vars.size() > 30) return;
    int nvars_total = static_cast<int>(var_to_checks.size());
    std::vector<char> var_in(nvars_total, 0);
    for (int v : vars) {
        if (v >= 0 && v < nvars_total) var_in[v] = 1;
    }
    std::vector<char> check_in(checks.size(), 0);
    for (int v : vars) {
        if (v < 0 || v >= nvars_total) continue;
        for (const auto &edge : var_to_checks[v]) {
            int c = edge.check;
            if (c >= 0 && c < static_cast<int>(checks.size())) {
                check_in[c] = 1;
            }
        }
    }
    std::vector<int> var_pos(nvars_total, -1);
    for (size_t i = 0; i < vars.size(); ++i) {
        int v = vars[i];
        if (v >= 0 && v < nvars_total) var_pos[v] = static_cast<int>(i);
    }
    const CycleList *cycle_list = (side_label == 'Z') ? &cycles->z : &cycles->x;
    if (!cycle_list) return;
    std::vector<CycleEntry> cycles8;
    std::vector<CycleEntry> cycles10;
    for (const auto &entry : cycle_list->cycles8) {
        bool ok = true;
        for (int i = 0; i < 4; ++i) {
            int v = entry.vars[i];
            int c = entry.checks[i];
            if (v < 0 || v >= nvars_total || !var_in[v]) {
                ok = false;
                break;
            }
            if (c < 0 || c >= static_cast<int>(check_in.size()) || !check_in[c]) {
                ok = false;
                break;
            }
        }
        if (ok) cycles8.push_back(entry);
    }
    for (const auto &entry : cycle_list->cycles10) {
        bool ok = true;
        for (int i = 0; i < 5; ++i) {
            int v = entry.vars[i];
            int c = entry.checks[i];
            if (v < 0 || v >= nvars_total || !var_in[v]) {
                ok = false;
                break;
            }
            if (c < 0 || c >= static_cast<int>(check_in.size()) || !check_in[c]) {
                ok = false;
                break;
            }
        }
        if (ok) cycles10.push_back(entry);
    }
    if (cycles8.empty() && cycles10.empty()) return;

    std::vector<CycleCoverItem> items;
    items.reserve(cycles8.size() + cycles10.size());
    auto add_cycle = [&](const CycleEntry &entry, int length, int id, int count) {
        uint64_t mask = 0;
        for (int i = 0; i < count; ++i) {
            int v = entry.vars[i];
            if (v < 0 || v >= nvars_total) return;
            int pos = var_pos[v];
            if (pos < 0) return;
            mask |= (1ULL << pos);
        }
        if (mask == 0) return;
        items.push_back(CycleCoverItem{mask, length, id, popcount_u64(mask)});
    };
    for (size_t i = 0; i < cycles8.size(); ++i) {
        add_cycle(cycles8[i], 8, static_cast<int>(i), 4);
    }
    for (size_t i = 0; i < cycles10.size(); ++i) {
        add_cycle(cycles10[i], 10, static_cast<int>(i), 5);
    }
    if (items.empty()) return;

    int var_count = static_cast<int>(vars.size());
    uint64_t target_mask = (1ULL << var_count) - 1ULL;
    uint64_t union_mask = 0;
    for (const auto &item : items) {
        union_mask |= item.mask;
    }
    if (union_mask != target_mask) return;
    int total_cycles = static_cast<int>(items.size());
    int max_cycle_vars = 0;
    for (const auto &item : items) {
        if (item.weight > max_cycle_vars) max_cycle_vars = item.weight;
    }
    if (max_cycle_vars <= 0) return;

    std::vector<int> cover_indices;
    bool found = false;
    if (var_count <= max_cycle_vars * 6) {
        std::vector<std::vector<int>> items_by_bit(var_count);
        for (size_t idx = 0; idx < items.size(); ++idx) {
            uint64_t mask = items[idx].mask;
            while (mask) {
                int bit = static_cast<int>(__builtin_ctzll(mask));
                items_by_bit[bit].push_back(static_cast<int>(idx));
                mask &= (mask - 1);
            }
        }
        for (int bit = 0; bit < var_count; ++bit) {
            if (items_by_bit[bit].empty()) {
                items_by_bit.clear();
                break;
            }
            auto &list = items_by_bit[bit];
            std::sort(list.begin(), list.end(),
                      [&](int a, int b) { return items[a].weight > items[b].weight; });
        }
        if (!items_by_bit.empty()) {
            for (int depth = 1; depth <= 6 && !found; ++depth) {
                cover_indices.clear();
                std::vector<char> used(items.size(), 0);
                if (cycle_cover_dfs(items, items_by_bit, target_mask, depth, 0, max_cycle_vars,
                                    used, cover_indices)) {
                    found = true;
                }
            }
        }
    }
    if (found) {
        log_trapset(side_label, vars, us, checks, var_to_checks, cycle_list, log_line);
        char side_lower = (side_label == 'Z') ? 'z' : 'x';
        std::ostringstream oss;
        oss << "trap_" << side_lower << "_cycle_union_count=" << cover_indices.size()
            << " cycles=[";
        for (size_t i = 0; i < cover_indices.size(); ++i) {
            if (i) oss << ",";
            const auto &item = items[cover_indices[i]];
            oss << item.length << ":" << item.id;
        }
        oss << "]";
        log_line(oss.str());
        return;
    }
    if (total_cycles < 7) return;

    struct CycleUnionInfo {
        int length = 0;
        int count = 0;
        std::array<int, 5> vars{};
        std::array<int, 5> checks{};
    };
    std::vector<CycleUnionInfo> union_cycles;
    union_cycles.reserve(cycles8.size() + cycles10.size());
    for (const auto &entry : cycles8) {
        CycleUnionInfo info;
        info.length = 8;
        info.count = 4;
        info.vars = entry.vars;
        info.checks = entry.checks;
        union_cycles.push_back(info);
    }
    for (const auto &entry : cycles10) {
        CycleUnionInfo info;
        info.length = 10;
        info.count = 5;
        info.vars = entry.vars;
        info.checks = entry.checks;
        union_cycles.push_back(info);
    }

    int ncycles = static_cast<int>(union_cycles.size());
    std::vector<std::vector<int>> adj(ncycles);
    int link_count = 0;
    auto shares = [&](const CycleUnionInfo &a, const CycleUnionInfo &b) {
        for (int i = 0; i < a.count; ++i) {
            int v = a.vars[i];
            for (int j = 0; j < b.count; ++j) {
                if (v == b.vars[j]) return true;
            }
        }
        for (int i = 0; i < a.count; ++i) {
            int c = a.checks[i];
            for (int j = 0; j < b.count; ++j) {
                if (c == b.checks[j]) return true;
            }
        }
        return false;
    };
    for (int i = 0; i < ncycles; ++i) {
        for (int j = i + 1; j < ncycles; ++j) {
            if (shares(union_cycles[i], union_cycles[j])) {
                adj[i].push_back(j);
                adj[j].push_back(i);
                link_count++;
            }
        }
    }

    std::vector<int> comp_id(ncycles, -1);
    struct CompInfo {
        int size = 0;
        int count8 = 0;
        int count10 = 0;
    };
    std::vector<CompInfo> comps;
    for (int i = 0; i < ncycles; ++i) {
        if (comp_id[i] != -1) continue;
        CompInfo info;
        std::vector<int> stack;
        stack.push_back(i);
        int comp_index = static_cast<int>(comps.size());
        comp_id[i] = comp_index;
        while (!stack.empty()) {
            int v = stack.back();
            stack.pop_back();
            info.size++;
            if (union_cycles[v].length == 8) info.count8++;
            if (union_cycles[v].length == 10) info.count10++;
            for (int nb : adj[v]) {
                if (comp_id[nb] == -1) {
                    comp_id[nb] = comp_index;
                    stack.push_back(nb);
                }
            }
        }
        comps.push_back(info);
    }
    std::sort(comps.begin(), comps.end(),
              [](const CompInfo &a, const CompInfo &b) { return a.size > b.size; });

    char side_lower = (side_label == 'Z') ? 'z' : 'x';
    {
        std::ostringstream oss;
        oss << "trap_" << side_lower << "_cycle_union_summary total_cycles=" << total_cycles
            << " cycle8=" << cycles8.size() << " cycle10=" << cycles10.size();
        log_line(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "trap_" << side_lower << "_cycle_union_links count=" << link_count;
        log_line(oss.str());
    }
    {
        std::ostringstream oss;
        oss << "trap_" << side_lower << "_cycle_union_components count=" << comps.size()
            << " sizes=[";
        for (size_t i = 0; i < comps.size(); ++i) {
            if (i) oss << ",";
            oss << comps[i].size;
        }
        oss << "]";
        log_line(oss.str());
    }
    for (size_t i = 0; i < comps.size(); ++i) {
        std::ostringstream oss;
        oss << "trap_" << side_lower << "_cycle_union_component id=" << i
            << " size=" << comps[i].size
            << " cycle8=" << comps[i].count8
            << " cycle10=" << comps[i].count10;
        log_line(oss.str());
    }
}

template <typename LogFn>
static void maybe_log_trapset(
    char side_label,
    const std::vector<int> &us,
    bool ets_applied,
    const std::vector<int> &diff,
    const std::vector<std::vector<int>> &checks,
    const std::vector<std::vector<EdgeRef>> &var_to_checks,
    const CycleIndex *cycles,
    LogFn log_line
) {
    if (ets_applied) return;
    if (us.size() != 2) return;
    if (diff.empty() || diff.size() > 12) return;
    const CycleList *cycle_list = nullptr;
    if (cycles) {
        cycle_list = (side_label == 'Z') ? &cycles->z : &cycles->x;
    }
    log_trapset(side_label, diff, us, checks, var_to_checks, cycle_list, log_line);
}

static std::vector<int> compute_syndrome(
    const std::vector<std::vector<int>> &checks,
    const std::vector<int> &err,
    bool use_xbit
) {
    std::vector<int> synd(checks.size(), 0);
    for (int c = 0; c < static_cast<int>(checks.size()); ++c) {
        int parity = 0;
        for (int v : checks[c]) {
            int bit = use_xbit ? xbit(err[v]) : zbit(err[v]);
            parity ^= bit;
        }
        synd[c] = parity;
    }
    return synd;
}

static Msg multiply_msg(const Msg &a, const Msg &b) {
    return Msg{a[0] * b[0], a[1] * b[1], a[2] * b[2], a[3] * b[3]};
}

static int max_rank_by_llr(
    const std::vector<int> &diff,
    const std::vector<double> &abs_llr
) {
    if (diff.empty()) return -1;
    std::vector<double> sorted = abs_llr;
    std::sort(sorted.begin(), sorted.end());
    int max_rank = 0;
    for (int v : diff) {
        if (v < 0 || v >= static_cast<int>(abs_llr.size())) continue;
        double val = abs_llr[v];
        int rank = static_cast<int>(std::upper_bound(sorted.begin(), sorted.end(), val) - sorted.begin());
        if (rank > max_rank) max_rank = rank;
    }
    return max_rank;
}

static const char *tf(bool v) {
    return v ? "true" : "false";
}

static Msg det_msg_xbit(int bit) {
    return bit ? Msg{0.0, 0.5, 0.0, 0.5} : Msg{0.5, 0.0, 0.5, 0.0};
}

static Msg det_msg_zbit(int bit) {
    return bit ? Msg{0.0, 0.0, 0.5, 0.5} : Msg{0.5, 0.5, 0.0, 0.0};
}

static void freeze_side_msgs(
    const std::vector<std::vector<int>> &checks,
    std::vector<std::vector<Msg>> &v2c,
    std::vector<std::vector<Msg>> &c2v,
    const std::vector<int> &est,
    bool use_xbit
) {
    for (int c = 0; c < static_cast<int>(checks.size()); ++c) {
        int deg = static_cast<int>(checks[c].size());
        for (int i = 0; i < deg; ++i) {
            int v = checks[c][i];
            int bit = use_xbit ? xbit(est[v]) : zbit(est[v]);
            Msg m = use_xbit ? det_msg_xbit(bit) : det_msg_zbit(bit);
            v2c[c][i] = m;
            c2v[c][i] = m;
        }
    }
}

static Msg scaled_mix(const Msg &new_msg, const Msg &old_msg, double damping) {
    if (damping <= 0.0) return new_msg;
    Msg out;
    for (int i = 0; i < 4; ++i) {
        out[i] = (1.0 - damping) * new_msg[i] + damping * old_msg[i];
    }
    normalize_msg(out);
    return out;
}

struct JointBPResult {
    std::vector<int> est;
    int iterations = 0;
    bool syndrome_match = false;
    bool bp_syndrome_match = false;
    bool pp_success = false;
    bool pp_success_ets = false;
    bool pp_success_flip = false;
    bool pp_used_ets_any = false;
    bool pp_used_ets6 = false;
    bool pp_used_ets12 = false;
    std::vector<std::string> pp_used_ets_labels;
    std::vector<double> abs_llr_x;
    std::vector<double> abs_llr_z;
    std::vector<int> flip_x;
    std::vector<int> flip_z;
    std::vector<int> flip_x_history;
    std::vector<int> flip_z_history;
    std::vector<std::vector<int>> flip_x_window;
    std::vector<std::vector<int>> flip_z_window;
};

static JointBPResult joint_bp_decode(
    const std::vector<std::vector<int>> &x_checks,
    const std::vector<std::vector<int>> &z_checks,
    const std::vector<std::vector<EdgeRef>> &var_to_x,
    const std::vector<std::vector<EdgeRef>> &var_to_z,
    const std::vector<int> &sx,
    const std::vector<int> &sz,
    const Msg &prior,
    int max_iter,
    int flip_hist_window,
    bool freeze_syn,
    double damping,
    bool verbose,
    bool verbose_all,
    bool collect_llr_final,
    const PPEarlyContext *pp_ctx,
    const std::vector<int> *truth_err,
    const RowBasis *hx_basis,
    const RowBasis *hz_basis,
    std::vector<std::string> *iter_logs,
    std::vector<std::string> *pp_logs
) {
    int mX = static_cast<int>(x_checks.size());
    int mZ = static_cast<int>(z_checks.size());
    int nvars = static_cast<int>(var_to_x.size());
    bool pp_verbose = (pp_ctx && pp_ctx->verbose);

    std::vector<std::vector<Msg>> x_v2c(mX), x_c2v(mX);
    for (int c = 0; c < mX; ++c) {
        int deg = static_cast<int>(x_checks[c].size());
        x_v2c[c].assign(deg, prior);
        x_c2v[c].assign(deg, Msg{0.25, 0.25, 0.25, 0.25});
    }

    std::vector<std::vector<Msg>> z_v2c(mZ), z_c2v(mZ);
    for (int c = 0; c < mZ; ++c) {
        int deg = static_cast<int>(z_checks[c].size());
        z_v2c[c].assign(deg, prior);
        z_c2v[c].assign(deg, Msg{0.25, 0.25, 0.25, 0.25});
    }

    std::vector<int> est(nvars, 0);
    std::vector<int> prev_est;
    std::vector<int> prev_prev_est;
    std::vector<int> last_flip_x;
    std::vector<int> last_flip_z;
    std::vector<double> abs_llr_x;
    std::vector<double> abs_llr_z;
    std::deque<std::vector<int>> flip_hist_x_window;
    std::deque<std::vector<int>> flip_hist_z_window;
    bool flip_hist_active = false;
    bool freeze_x_active = false;
    bool freeze_z_active = false;
    int stall_same = 0;
    bool stall_active = false;
    bool oscillation_active = false;
    bool us_oscillation_active = false;
    std::deque<int> us_signal_hist;
    const int stall_window = 1;
    const int hist_window = (flip_hist_window > 0) ? flip_hist_window : 0;
    const int pp_warmup_iters = 100;
    const int us_osc_window = 100;
    const int us_osc_min_turns = 2;
    const int us_osc_min_range = 10;
    const double us_osc_min_range_frac = 0.2;
    const int pp_us_limit = 20;
    const int osd_us_limit = 20;
    bool pp_success = false;
    bool pp_success_ets = false;
    bool pp_success_flip = false;
    bool pp_used_ets_any = false;
    bool pp_used_ets6 = false;
    bool pp_used_ets12 = false;
    std::vector<std::string> pp_used_ets_labels;
    bool keep_llr = collect_llr_final;

    auto record_ets_label = [&](const ETSFixResult &res) {
        if (!res.applied) return;
        if (!res.label.empty()) {
            pp_used_ets_labels.push_back(res.label);
        } else if (res.size > 0) {
            pp_used_ets_labels.push_back(std::to_string(res.size));
        }
    };

    auto detect_us_oscillation = [&](const std::deque<int> &hist,
                                     int &range_out,
                                     int &turns_out,
                                     int &trend_out,
                                     double &mean_out) -> bool {
        if (static_cast<int>(hist.size()) < us_osc_window) return false;
        int min_v = hist.front();
        int max_v = hist.front();
        long long sum = 0;
        for (int v : hist) {
            if (v < min_v) min_v = v;
            if (v > max_v) max_v = v;
            sum += v;
        }
        mean_out = sum / static_cast<double>(hist.size());
        range_out = max_v - min_v;
        if (mean_out <= 0.0) return false;
        int min_range = std::max(us_osc_min_range,
                                 static_cast<int>(mean_out * us_osc_min_range_frac));
        if (range_out < min_range) return false;
        trend_out = hist.back() - hist.front();
        int prev_sign = 0;
        turns_out = 0;
        for (size_t i = 1; i < hist.size(); ++i) {
            int diff = hist[i] - hist[i - 1];
            int sign = (diff > 0) ? 1 : (diff < 0) ? -1 : 0;
            if (sign != 0) {
                if (prev_sign != 0 && sign != prev_sign) {
                    turns_out++;
                }
                prev_sign = sign;
            }
        }
        return turns_out >= us_osc_min_turns;
    };

    auto make_result = [&](int iters, bool bp_ok, bool final_ok) {
        JointBPResult out;
        out.est = est;
        out.iterations = iters;
        out.syndrome_match = final_ok;
        out.bp_syndrome_match = bp_ok;
        out.pp_success = pp_success;
        out.pp_success_ets = pp_success_ets;
        out.pp_success_flip = pp_success_flip;
        out.pp_used_ets_any = pp_used_ets_any;
        out.pp_used_ets6 = pp_used_ets6;
        out.pp_used_ets12 = pp_used_ets12;
        if (pp_success) {
            out.pp_used_ets_labels = sorted_unique_strings(pp_used_ets_labels);
        }
        if (keep_llr && abs_llr_x.size() == static_cast<size_t>(nvars)) {
            out.abs_llr_x = abs_llr_x;
            out.abs_llr_z = abs_llr_z;
        }
        out.flip_x = last_flip_x;
        out.flip_z = last_flip_z;
        if (flip_hist_active && hist_window > 0) {
            std::vector<char> flip_hist_x(nvars, 0);
            std::vector<char> flip_hist_z(nvars, 0);
            for (const auto &vec : flip_hist_x_window) {
                for (int v : vec) {
                    if (v >= 0 && v < nvars) flip_hist_x[v] = 1;
                }
            }
            for (const auto &vec : flip_hist_z_window) {
                for (int v : vec) {
                    if (v >= 0 && v < nvars) flip_hist_z[v] = 1;
                }
            }
            out.flip_x_history.reserve(nvars / 8);
            out.flip_z_history.reserve(nvars / 8);
            for (int i = 0; i < nvars; ++i) {
                if (flip_hist_x[i]) out.flip_x_history.push_back(i);
            }
            for (int i = 0; i < nvars; ++i) {
                if (flip_hist_z[i]) out.flip_z_history.push_back(i);
            }
            out.flip_x_window.assign(flip_hist_x_window.begin(), flip_hist_x_window.end());
            out.flip_z_window.assign(flip_hist_z_window.begin(), flip_hist_z_window.end());
        }
        return out;
    };

    struct StatusMetrics {
        bool syn_all = false;
        int us_x = 0;
        int us_z = 0;
        bool have_diff = false;
        int diff_x_count = -1;
        int diff_z_count = -1;
        int min_diff_x_count = -1;
        int min_diff_z_count = -1;
        bool have_row = false;
        bool diff_row_x = false;
        bool diff_row_z = false;
        bool have_llr_rank = false;
        int diff_x_rank_max = -1;
        int diff_z_rank_max = -1;
        int flip_x_count = -1;
        int flip_z_count = -1;
        int min_us_x_nonzero = -1;
        int min_us_z_nonzero = -1;
        int min_diff_x_rank = -1;
        int min_diff_z_rank = -1;
    };

    auto compute_status_metrics = [&](const std::vector<int> &cur_est,
                                      const std::vector<int> &sx_hat_local,
                                      const std::vector<int> &sz_hat_local,
                                      const std::vector<int> *prev_for_flip,
                                      bool collect_llr_flag,
                                      StatusMetrics &out) {
        out.syn_all = (sx_hat_local == sx) && (sz_hat_local == sz);
        out.us_x = unsatisfied_count(sx_hat_local, sx);
        out.us_z = unsatisfied_count(sz_hat_local, sz);
        out.flip_x_count = -1;
        out.flip_z_count = -1;
        if (prev_for_flip && prev_for_flip->size() == cur_est.size()) {
            int fx = 0;
            int fz = 0;
            for (size_t i = 0; i < cur_est.size(); ++i) {
                if (xbit((*prev_for_flip)[i]) != xbit(cur_est[i])) fx++;
                if (zbit((*prev_for_flip)[i]) != zbit(cur_est[i])) fz++;
            }
            out.flip_x_count = fx;
            out.flip_z_count = fz;
        }
        out.have_diff = (truth_err != nullptr && truth_err->size() == cur_est.size());
        out.diff_x_count = -1;
        out.diff_z_count = -1;
        out.have_row = false;
        out.diff_row_x = false;
        out.diff_row_z = false;
        out.have_llr_rank = false;
        out.diff_x_rank_max = -1;
        out.diff_z_rank_max = -1;
        if (!out.have_diff) return;

        out.diff_x_count = 0;
        out.diff_z_count = 0;
        out.have_row = (hx_basis != nullptr && hz_basis != nullptr &&
                        hx_basis->nvars == nvars && hz_basis->nvars == nvars);
        out.have_llr_rank = collect_llr_flag &&
                            abs_llr_x.size() == static_cast<size_t>(nvars) &&
                            abs_llr_z.size() == static_cast<size_t>(nvars);
        std::vector<double> sorted_x;
        std::vector<double> sorted_z;
        if (out.have_llr_rank) {
            sorted_x = abs_llr_x;
            sorted_z = abs_llr_z;
            std::sort(sorted_x.begin(), sorted_x.end());
            std::sort(sorted_z.begin(), sorted_z.end());
        }
        std::vector<uint64_t> diff_x_bits;
        std::vector<uint64_t> diff_z_bits;
        if (out.have_row) {
            diff_x_bits.assign(hz_basis->words, 0);
            diff_z_bits.assign(hx_basis->words, 0);
        }
        size_t n = std::min(truth_err->size(), cur_est.size());
        for (size_t i = 0; i < n; ++i) {
            int dx = xbit((*truth_err)[i]) ^ xbit(cur_est[i]);
            int dz = zbit((*truth_err)[i]) ^ zbit(cur_est[i]);
            if (dx) {
                out.diff_x_count++;
                if (out.have_row) {
                    diff_x_bits[i >> 6] |= (1ULL << (i & 63));
                }
                if (out.have_llr_rank) {
                    double val = abs_llr_x[i];
                    int rank = static_cast<int>(std::upper_bound(sorted_x.begin(),
                                                                 sorted_x.end(),
                                                                 val) - sorted_x.begin());
                    if (rank > out.diff_x_rank_max) out.diff_x_rank_max = rank;
                }
            }
            if (dz) {
                out.diff_z_count++;
                if (out.have_row) {
                    diff_z_bits[i >> 6] |= (1ULL << (i & 63));
                }
                if (out.have_llr_rank) {
                    double val = abs_llr_z[i];
                    int rank = static_cast<int>(std::upper_bound(sorted_z.begin(),
                                                                 sorted_z.end(),
                                                                 val) - sorted_z.begin());
                    if (rank > out.diff_z_rank_max) out.diff_z_rank_max = rank;
                }
            }
        }
        if (out.have_row) {
            out.diff_row_x = in_row_space(*hz_basis, diff_x_bits);
            out.diff_row_z = in_row_space(*hx_basis, diff_z_bits);
        }
    };

    auto format_status_en = [&](int iter_idx, const StatusMetrics &s) {
        std::ostringstream oss;
        oss << "iter=" << iter_idx
            << " | syndrome_match=" << tf(s.syn_all)
            << " | us=(" << s.us_x << "," << s.us_z << ")";
        oss << " | us_min=(";
        if (s.min_us_x_nonzero >= 0) {
            oss << s.min_us_x_nonzero;
        } else {
            oss << "NA";
        }
        oss << ",";
        if (s.min_us_z_nonzero >= 0) {
            oss << s.min_us_z_nonzero;
        } else {
            oss << "NA";
        }
        oss << ") | diff=(";
        if (s.have_diff) {
            oss << s.diff_x_count << "," << s.diff_z_count << ")";
        } else {
            oss << "NA,NA)";
        }
        oss << " | diff_min=(";
        if (s.min_diff_x_count >= 0) {
            oss << s.min_diff_x_count;
        } else {
            oss << "NA";
        }
        oss << ",";
        if (s.min_diff_z_count >= 0) {
            oss << s.min_diff_z_count;
        } else {
            oss << "NA";
        }
        oss << ")";
        if (s.have_diff) {
            if (s.have_row) {
                oss << " | diff_row=(" << tf(s.diff_row_x) << "," << tf(s.diff_row_z) << ")";
            } else {
                oss << " | diff_row=(NA,NA)";
            }
        } else {
            oss << " | diff_row=(NA,NA)";
        }
        oss << " | diff_llr=(";
        if (s.have_diff && s.have_llr_rank) {
            if (s.diff_x_rank_max >= 0) {
                oss << s.diff_x_rank_max;
            } else {
                oss << "NA";
            }
            oss << ",";
            if (s.diff_z_rank_max >= 0) {
                oss << s.diff_z_rank_max;
            } else {
                oss << "NA";
            }
        } else {
            oss << "NA,NA";
        }
        oss << ")";
        oss << " | diff_llr_min=(";
        if (s.min_diff_x_rank >= 0) {
            oss << s.min_diff_x_rank;
        } else {
            oss << "NA";
        }
        oss << ",";
        if (s.min_diff_z_rank >= 0) {
            oss << s.min_diff_z_rank;
        } else {
            oss << "NA";
        }
        oss << ")";
        oss << " | flip=(";
        if (s.flip_x_count >= 0) {
            oss << s.flip_x_count;
        } else {
            oss << "NA";
        }
        oss << ",";
        if (s.flip_z_count >= 0) {
            oss << s.flip_z_count;
        } else {
            oss << "NA";
        }
        oss << ")";
        return oss.str();
    };

    auto print_status_en = [&](int iter_idx, const StatusMetrics &s) {
        std::cout << format_status_en(iter_idx, s) << "\n";
    };

    int min_us_x_nonzero = -1;
    int min_us_z_nonzero = -1;
    int min_diff_x_count = -1;
    int min_diff_z_count = -1;
    int min_diff_x_rank = -1;
    int min_diff_z_rank = -1;
    auto update_status_min = [&](StatusMetrics &s) {
        if (s.us_x == 0) {
            min_us_x_nonzero = 0;
        } else if (s.us_x > 0 && (min_us_x_nonzero < 0 || s.us_x < min_us_x_nonzero)) {
            min_us_x_nonzero = s.us_x;
        }
        if (s.us_z == 0) {
            min_us_z_nonzero = 0;
        } else if (s.us_z > 0 && (min_us_z_nonzero < 0 || s.us_z < min_us_z_nonzero)) {
            min_us_z_nonzero = s.us_z;
        }
        if (s.have_diff) {
            if (s.diff_x_count >= 0 &&
                (min_diff_x_count < 0 || s.diff_x_count < min_diff_x_count)) {
                min_diff_x_count = s.diff_x_count;
            }
            if (s.diff_z_count >= 0 &&
                (min_diff_z_count < 0 || s.diff_z_count < min_diff_z_count)) {
                min_diff_z_count = s.diff_z_count;
            }
        }
        if (s.have_diff && s.have_llr_rank) {
            if (s.diff_x_rank_max >= 0 &&
                (min_diff_x_rank < 0 || s.diff_x_rank_max < min_diff_x_rank)) {
                min_diff_x_rank = s.diff_x_rank_max;
            }
            if (s.diff_z_rank_max >= 0 &&
                (min_diff_z_rank < 0 || s.diff_z_rank_max < min_diff_z_rank)) {
                min_diff_z_rank = s.diff_z_rank_max;
            }
        }
        s.min_us_x_nonzero = min_us_x_nonzero;
        s.min_us_z_nonzero = min_us_z_nonzero;
        s.min_diff_x_count = min_diff_x_count;
        s.min_diff_z_count = min_diff_z_count;
        s.min_diff_x_rank = min_diff_x_rank;
        s.min_diff_z_rank = min_diff_z_rank;
    };

    bool want_llr_diag = (truth_err != nullptr && (verbose || (iter_logs != nullptr)));
    bool want_llr_pp = (pp_ctx && pp_ctx->enable);
    for (int iter = 0; iter < max_iter; ++iter) {
        bool collect_llr = want_llr_diag || want_llr_pp || (collect_llr_final && (iter + 1 == max_iter));
        if (collect_llr) {
            abs_llr_x.assign(nvars, 0.0);
            abs_llr_z.assign(nvars, 0.0);
        }
        if (!freeze_x_active) {
            for (int c = 0; c < mX; ++c) {
                int deg = static_cast<int>(x_checks[c].size());
                std::vector<double> q0(deg, 0.0), q1(deg, 0.0);
                for (int i = 0; i < deg; ++i) {
                    const Msg &m = x_v2c[c][i];
                    q0[i] = m[0] + m[2];
                    q1[i] = m[1] + m[3];
                }
                for (int i = 0; i < deg; ++i) {
                    double p_even = 1.0;
                    double p_odd = 0.0;
                    for (int j = 0; j < deg; ++j) {
                        if (j == i) continue;
                        double new_even = p_even * q0[j] + p_odd * q1[j];
                        double new_odd = p_even * q1[j] + p_odd * q0[j];
                        p_even = new_even;
                        p_odd = new_odd;
                    }
                    double val0 = (sx[c] == 0) ? p_even : p_odd;
                    double val1 = (sx[c] == 0) ? p_odd : p_even;
                    Msg out{val0, val1, val0, val1};
                    normalize_msg(out);
                    x_c2v[c][i] = out;
                }
            }
        }

        if (!freeze_z_active) {
            for (int c = 0; c < mZ; ++c) {
                int deg = static_cast<int>(z_checks[c].size());
                std::vector<double> q0(deg, 0.0), q1(deg, 0.0);
                for (int i = 0; i < deg; ++i) {
                    const Msg &m = z_v2c[c][i];
                    q0[i] = m[0] + m[1];
                    q1[i] = m[2] + m[3];
                }
                for (int i = 0; i < deg; ++i) {
                    double p_even = 1.0;
                    double p_odd = 0.0;
                    for (int j = 0; j < deg; ++j) {
                        if (j == i) continue;
                        double new_even = p_even * q0[j] + p_odd * q1[j];
                        double new_odd = p_even * q1[j] + p_odd * q0[j];
                        p_even = new_even;
                        p_odd = new_odd;
                    }
                    double val0 = (sz[c] == 0) ? p_even : p_odd;
                    double val1 = (sz[c] == 0) ? p_odd : p_even;
                    Msg out{val0, val0, val1, val1};
                    normalize_msg(out);
                    z_c2v[c][i] = out;
                }
            }
        }

        for (int v = 0; v < nvars; ++v) {
            Msg total = prior;
            for (const auto &e : var_to_x[v]) {
                total = multiply_msg(total, x_c2v[e.check][e.pos]);
            }
            for (const auto &e : var_to_z[v]) {
                total = multiply_msg(total, z_c2v[e.check][e.pos]);
            }
            normalize_msg(total);

            int best = 0;
            for (int s = 1; s < 4; ++s) {
                if (total[s] > total[best]) best = s;
            }
            est[v] = best;
            if (collect_llr) {
                const double eps = 1e-300;
                double px1 = total[1] + total[3];
                double px0 = total[0] + total[2];
                double pz1 = total[2] + total[3];
                double pz0 = total[0] + total[1];
                double llr_x = std::log(std::max(px1, eps)) - std::log(std::max(px0, eps));
                double llr_z = std::log(std::max(pz1, eps)) - std::log(std::max(pz0, eps));
                abs_llr_x[v] = std::abs(llr_x);
                abs_llr_z[v] = std::abs(llr_z);
            }

            if (!freeze_x_active) {
                for (const auto &e : var_to_x[v]) {
                    Msg out = prior;
                    for (const auto &e2 : var_to_x[v]) {
                        if (e2.check == e.check && e2.pos == e.pos) continue;
                        out = multiply_msg(out, x_c2v[e2.check][e2.pos]);
                    }
                    for (const auto &e2 : var_to_z[v]) {
                        out = multiply_msg(out, z_c2v[e2.check][e2.pos]);
                    }
                    normalize_msg(out);
                    x_v2c[e.check][e.pos] = scaled_mix(out, x_v2c[e.check][e.pos], damping);
                }
            }
            if (!freeze_z_active) {
                for (const auto &e : var_to_z[v]) {
                    Msg out = prior;
                    for (const auto &e2 : var_to_x[v]) {
                        out = multiply_msg(out, x_c2v[e2.check][e2.pos]);
                    }
                    for (const auto &e2 : var_to_z[v]) {
                        if (e2.check == e.check && e2.pos == e.pos) continue;
                        out = multiply_msg(out, z_c2v[e2.check][e2.pos]);
                    }
                    normalize_msg(out);
                    z_v2c[e.check][e.pos] = scaled_mix(out, z_v2c[e.check][e.pos], damping);
                }
            }
        }

        auto sx_hat = compute_syndrome(x_checks, est, true);
        auto sz_hat = compute_syndrome(z_checks, est, false);
        bool syn_x = (sx_hat == sx);
        bool syn_z = (sz_hat == sz);
        bool syn_all = syn_x && syn_z;
        if (freeze_syn) {
            if (syn_x && !freeze_x_active) {
                freeze_x_active = true;
                freeze_side_msgs(x_checks, x_v2c, x_c2v, est, true);
            }
            if (syn_z && !freeze_z_active) {
                freeze_z_active = true;
                freeze_side_msgs(z_checks, z_v2c, z_c2v, est, false);
            }
        }
        int us_x_count = 0;
        int us_z_count = 0;
        int us_signal = 0;
        bool pp_skip_due_us = false;
        if (!syn_all) {
            us_x_count = unsatisfied_count(sx_hat, sx);
            us_z_count = unsatisfied_count(sz_hat, sz);
            if (us_x_count == 0) {
                min_us_x_nonzero = 0;
            } else if (us_x_count > 0 && (min_us_x_nonzero < 0 || us_x_count < min_us_x_nonzero)) {
                min_us_x_nonzero = us_x_count;
            }
            if (us_z_count == 0) {
                min_us_z_nonzero = 0;
            } else if (us_z_count > 0 && (min_us_z_nonzero < 0 || us_z_count < min_us_z_nonzero)) {
                min_us_z_nonzero = us_z_count;
            }
            us_signal = std::max(us_x_count, us_z_count);
            pp_skip_due_us = (us_signal > pp_us_limit);
            us_signal_hist.push_back(us_signal);
            if (static_cast<int>(us_signal_hist.size()) > us_osc_window) {
                us_signal_hist.pop_front();
            }
        }
        bool same_est = (!prev_est.empty() && est == prev_est);
        bool oscillating = (!prev_prev_est.empty() && est == prev_prev_est && !same_est);
        if (!syn_all) {
            bool allow_pp_trigger = (iter + 1 >= pp_warmup_iters) || (iter + 1 == max_iter);
            if (same_est) {
                stall_same++;
            } else {
                stall_same = 0;
            }
            if (!stall_active && (stall_same >= stall_window || iter + 1 == max_iter) && allow_pp_trigger) {
                stall_active = true;
                if (!flip_hist_active && hist_window > 0) {
                    flip_hist_active = true;
                    flip_hist_x_window.clear();
                    flip_hist_z_window.clear();
                }
            }
            if (!oscillation_active && oscillating && allow_pp_trigger) {
                oscillation_active = true;
                if (!flip_hist_active && hist_window > 0) {
                    flip_hist_active = true;
                    flip_hist_x_window.clear();
                    flip_hist_z_window.clear();
                }
            }
            if (!us_oscillation_active && allow_pp_trigger) {
                int range = 0;
                int turns = 0;
                int trend = 0;
                double mean = 0.0;
                if (detect_us_oscillation(us_signal_hist, range, turns, trend, mean)) {
                    us_oscillation_active = true;
                    if (!flip_hist_active && hist_window > 0) {
                        flip_hist_active = true;
                        flip_hist_x_window.clear();
                        flip_hist_z_window.clear();
                    }
                    if (verbose) {
                        std::cout << "[PP] iter=" << (iter + 1)
                                  << " trigger=us_oscillation"
                                  << " us_range=" << range
                                  << " us_turns=" << turns
                                  << " us_trend=" << trend
                                  << "\n";
                    }
                }
            }
        }

        if (!prev_est.empty()) {
            size_t n = std::min(prev_est.size(), est.size());
            last_flip_x.clear();
            last_flip_z.clear();
            last_flip_x.reserve(n / 8);
            last_flip_z.reserve(n / 8);
            for (size_t i = 0; i < n; ++i) {
                if (xbit(prev_est[i]) != xbit(est[i])) {
                    last_flip_x.push_back(static_cast<int>(i));
                }
                if (zbit(prev_est[i]) != zbit(est[i])) {
                    last_flip_z.push_back(static_cast<int>(i));
                }
            }
            if (flip_hist_active) {
                flip_hist_x_window.push_back(last_flip_x);
                if (static_cast<int>(flip_hist_x_window.size()) > hist_window) {
                    flip_hist_x_window.pop_front();
                }
                flip_hist_z_window.push_back(last_flip_z);
                if (static_cast<int>(flip_hist_z_window.size()) > hist_window) {
                    flip_hist_z_window.pop_front();
                }
            }
        }
        bool need_detail = (verbose_all || stall_active || oscillation_active || us_oscillation_active) &&
                           (verbose || (iter_logs != nullptr));
        if (need_detail) {
            StatusMetrics status;
            compute_status_metrics(est, sx_hat, sz_hat,
                                   prev_est.empty() ? nullptr : &prev_est,
                                   collect_llr, status);
            update_status_min(status);
            if (verbose) {
                print_status_en(iter + 1, status);
            }
            if (iter_logs) {
                std::ostringstream oss;
                oss << "iter=" << (iter + 1)
                    << " syndrome_match=" << (status.syn_all ? "1" : "0")
                    << " us=(" << status.us_x << "," << status.us_z << ")";
                oss << " us_min=(";
                if (status.min_us_x_nonzero >= 0) {
                    oss << status.min_us_x_nonzero;
                } else {
                    oss << "NA";
                }
                oss << ",";
                if (status.min_us_z_nonzero >= 0) {
                    oss << status.min_us_z_nonzero;
                } else {
                    oss << "NA";
                }
                oss << ")";
                if (status.have_diff) {
                    oss << " diff=(" << status.diff_x_count << "," << status.diff_z_count << ")";
                    oss << " diff_min=(";
                    if (status.min_diff_x_count >= 0) {
                        oss << status.min_diff_x_count;
                    } else {
                        oss << "NA";
                    }
                    oss << ",";
                    if (status.min_diff_z_count >= 0) {
                        oss << status.min_diff_z_count;
                    } else {
                        oss << "NA";
                    }
                    oss << ")";
                    if (status.have_row) {
                        oss << " diff_row=(" << tf(status.diff_row_x)
                            << "," << tf(status.diff_row_z) << ")";
                    } else {
                        oss << " diff_row=(NA,NA)";
                    }
                } else {
                    oss << " diff=(NA,NA) diff_min=(NA,NA) diff_row=(NA,NA)";
                }
                oss << " diff_llr=(";
                if (status.have_diff && status.have_llr_rank) {
                    if (status.diff_x_rank_max >= 0) {
                        oss << status.diff_x_rank_max;
                    } else {
                        oss << "NA";
                    }
                    oss << ",";
                    if (status.diff_z_rank_max >= 0) {
                        oss << status.diff_z_rank_max;
                    } else {
                        oss << "NA";
                    }
                } else {
                    oss << "NA,NA";
                }
                oss << ")";
                oss << " diff_llr_min=(";
                if (status.min_diff_x_rank >= 0) {
                    oss << status.min_diff_x_rank;
                } else {
                    oss << "NA";
                }
                oss << ",";
                if (status.min_diff_z_rank >= 0) {
                    oss << status.min_diff_z_rank;
                } else {
                    oss << "NA";
                }
                oss << ")";
                oss << " flip=(";
                if (status.flip_x_count >= 0) {
                    oss << status.flip_x_count;
                } else {
                    oss << "NA";
                }
                oss << ",";
                if (status.flip_z_count >= 0) {
                    oss << status.flip_z_count;
                } else {
                    oss << "NA";
                }
                oss << ")";
                iter_logs->push_back(oss.str());
            }
        }
        if (syn_all) {
            return make_result(iter + 1, true, true);
        }
        if (pp_ctx && pp_ctx->enable &&
            (stall_active || oscillation_active || us_oscillation_active)) {
            if (pp_skip_due_us) {
                std::ostringstream oss;
                oss << "[PP] iter=" << (iter + 1)
                    << " skipped us=(" << us_x_count << "," << us_z_count << ")";
                std::string line = oss.str();
                if (pp_verbose) {
                    std::cout << line << "\n";
                }
                if (pp_logs) {
                    pp_logs->push_back(line);
                }
                prev_prev_est = prev_est;
                prev_est = est;
                continue;
            }
            auto pp_log_line = [&](const std::string &line) {
                if (pp_verbose) {
                    std::cout << line << "\n";
                }
                if (pp_logs) {
                    pp_logs->push_back(line);
                }
            };
            std::vector<int> trial_est = est;
            bool ok = false;
            bool ets_pp_ok = false;
            bool flip_pp_ok = false;
            pp_log_line("[PP] iter=" + std::to_string(iter + 1) + " start");
            auto print_ranked = [&](const std::string &label,
                                    const std::vector<int> &idxs,
                                    const std::vector<double> &abs_llr,
                                    const std::vector<double> &sorted_llr) {
                std::vector<std::pair<int, int>> ranked;
                ranked.reserve(idxs.size());
                for (int v : idxs) {
                    if (v < 0 || v >= static_cast<int>(abs_llr.size())) continue;
                    double val = abs_llr[v];
                    int rank = static_cast<int>(std::upper_bound(sorted_llr.begin(),
                                                                 sorted_llr.end(),
                                                                 val) - sorted_llr.begin());
                    ranked.push_back({rank, v});
                }
                std::sort(ranked.begin(), ranked.end(),
                          [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
                              if (a.first != b.first) return a.first < b.first;
                              return a.second < b.second;
                          });
                std::ostringstream oss;
                oss << "[PP] iter=" << (iter + 1) << " " << label << "=[";
                for (size_t i = 0; i < ranked.size(); ++i) {
                    if (i) oss << ",";
                    oss << ranked[i].second << ":" << ranked[i].first;
                }
                oss << "]";
                pp_log_line(oss.str());
            };
            bool have_llr_x = abs_llr_x.size() == static_cast<size_t>(nvars);
            bool have_llr_z = abs_llr_z.size() == static_cast<size_t>(nvars);
            std::vector<double> sorted_x;
            std::vector<double> sorted_z;
            if (have_llr_x) {
                sorted_x = abs_llr_x;
                std::sort(sorted_x.begin(), sorted_x.end());
            }
            if (have_llr_z) {
                sorted_z = abs_llr_z;
                std::sort(sorted_z.begin(), sorted_z.end());
            }
            auto us_x = unsatisfied_indices(sx_hat, sx);
            auto us_z = unsatisfied_indices(sz_hat, sz);
            {
                std::ostringstream us_oss;
                us_oss << "[PP-ETS] iter=" << (iter + 1) << " ";
                write_index_list(us_oss, "us_x", us_x);
                pp_log_line(us_oss.str());
            }
            {
                std::ostringstream us_oss;
                us_oss << "[PP-ETS] iter=" << (iter + 1) << " ";
                write_index_list(us_oss, "us_z", us_z);
                pp_log_line(us_oss.str());
            }
            std::vector<int> nus_x;
            std::vector<int> nus_z;
            if (!us_x.empty()) {
                std::vector<char> seen(nvars, 0);
                for (int c : us_x) {
                    if (c < 0 || c >= static_cast<int>(x_checks.size())) continue;
                    for (int v : x_checks[c]) {
                        if (v >= 0 && v < nvars) seen[v] = 1;
                    }
                }
                for (int i = 0; i < nvars; ++i) {
                    if (seen[i]) nus_x.push_back(i);
                }
            }
            if (!us_z.empty()) {
                std::vector<char> seen(nvars, 0);
                for (int c : us_z) {
                    if (c < 0 || c >= static_cast<int>(z_checks.size())) continue;
                    for (int v : z_checks[c]) {
                        if (v >= 0 && v < nvars) seen[v] = 1;
                    }
                }
                for (int i = 0; i < nvars; ++i) {
                    if (seen[i]) nus_z.push_back(i);
                }
            }
            if (have_llr_x && nus_x.size() <= 30) {
                print_ranked("nus_x_ranked", nus_x, abs_llr_x, sorted_x);
            }
            if (have_llr_z && nus_z.size() <= 30) {
                print_ranked("nus_z_ranked", nus_z, abs_llr_z, sorted_z);
            }
            std::vector<int> diff_x;
            std::vector<int> diff_z;
            if (truth_err != nullptr) {
                diff_x = diff_indices_bits(*truth_err, est, true);
                diff_z = diff_indices_bits(*truth_err, est, false);
                std::vector<char> flip_seen_x(nvars, 0);
                std::vector<char> flip_seen_z(nvars, 0);
                if (flip_hist_active && hist_window > 0) {
                    for (const auto &vec : flip_hist_x_window) {
                        for (int v : vec) {
                            if (v >= 0 && v < nvars) flip_seen_x[v] = 1;
                        }
                    }
                    for (const auto &vec : flip_hist_z_window) {
                        for (int v : vec) {
                            if (v >= 0 && v < nvars) flip_seen_z[v] = 1;
                        }
                    }
                }
                if (have_llr_x && diff_x.size() <= 30) {
                    print_ranked("diff_x_ranked", diff_x, abs_llr_x, sorted_x);
                    std::vector<int> diff_x_in_flip;
                    diff_x_in_flip.reserve(diff_x.size());
                    for (int v : diff_x) {
                        if (v >= 0 && v < nvars && flip_seen_x[v]) {
                            diff_x_in_flip.push_back(v);
                        }
                    }
                    print_ranked("diff_x_in_flip_ranked", diff_x_in_flip, abs_llr_x, sorted_x);
                }
                if (have_llr_z && diff_z.size() <= 30) {
                    print_ranked("diff_z_ranked", diff_z, abs_llr_z, sorted_z);
                    std::vector<int> diff_z_in_flip;
                    diff_z_in_flip.reserve(diff_z.size());
                    for (int v : diff_z) {
                        if (v >= 0 && v < nvars && flip_seen_z[v]) {
                            diff_z_in_flip.push_back(v);
                        }
                    }
                    print_ranked("diff_z_in_flip_ranked", diff_z_in_flip, abs_llr_z, sorted_z);
                }
                if (diff_x.size() <= 50) {
                    std::ostringstream diff_oss;
                    diff_oss << "[PP] iter=" << (iter + 1) << " ";
                    write_index_list(diff_oss, "diff_x", diff_x);
                    pp_log_line(diff_oss.str());
                }
                if (diff_z.size() <= 50) {
                    std::ostringstream diff_oss;
                    diff_oss << "[PP] iter=" << (iter + 1) << " ";
                    write_index_list(diff_oss, "diff_z", diff_z);
                    pp_log_line(diff_oss.str());
                }
                maybe_log_diff_trapset_cycle_union('X', diff_x, us_x, x_checks, var_to_x,
                                                   pp_ctx ? pp_ctx->cycles : nullptr, pp_log_line);
                maybe_log_diff_trapset_cycle_union('Z', diff_z, us_z, z_checks, var_to_z,
                                                   pp_ctx ? pp_ctx->cycles : nullptr, pp_log_line);
            }
            ETSFixResult ets_x = attempt_ets_fix_side_multi(
                "X", pp_ctx->ets, x_checks, var_to_x, sx, true, trial_est, pp_ctx->verbose
            );
            ETSFixResult ets_z = attempt_ets_fix_side_multi(
                "Z", pp_ctx->ets, z_checks, var_to_z, sz, false, trial_est, pp_ctx->verbose
            );
            bool fixed_x = ets_x.applied;
            bool fixed_z = ets_z.applied;
            bool ets_used = fixed_x || fixed_z;
            bool ets6 = (ets_x.size == 6) || (ets_z.size == 6);
            bool ets12 = (ets_x.size == 12) || (ets_z.size == 12);
            if (ets_used) {
                record_ets_label(ets_x);
                record_ets_label(ets_z);
            }
            if (truth_err != nullptr) {
                maybe_log_trapset('X', us_x, fixed_x, diff_x, x_checks, var_to_x,
                                  pp_ctx ? pp_ctx->cycles : nullptr, pp_log_line);
                maybe_log_trapset('Z', us_z, fixed_z, diff_z, z_checks, var_to_z,
                                  pp_ctx ? pp_ctx->cycles : nullptr, pp_log_line);
            }
            if (ets_used) {
                auto sx_pp = compute_syndrome(x_checks, trial_est, true);
                auto sz_pp = compute_syndrome(z_checks, trial_est, false);
                ok = (sx_pp == sx) && (sz_pp == sz);
                if (ok) {
                    ets_pp_ok = true;
                }
                std::ostringstream oss;
                oss << "[PP-ETS] iter=" << (iter + 1)
                    << " done pp-ets=" << tf(ets_used)
                    << " x=" << format_ets_label(ets_x)
                    << " z=" << format_ets_label(ets_z)
                    << " syndrome_ok=" << tf(ok);
                pp_log_line(oss.str());
                if (fixed_x && !us_x.empty()) {
                    std::ostringstream us_oss;
                    us_oss << "[PP-ETS] iter=" << (iter + 1) << " ";
                    write_index_list(us_oss, "us_x", us_x);
                    pp_log_line(us_oss.str());
                }
                if (fixed_z && !us_z.empty()) {
                    std::ostringstream us_oss;
                    us_oss << "[PP-ETS] iter=" << (iter + 1) << " ";
                    write_index_list(us_oss, "us_z", us_z);
                    pp_log_line(us_oss.str());
                }
            } else {
                pp_log_line("[PP-ETS] iter=" + std::to_string(iter + 1) +
                            " done pp-ets=false x=none z=none syndrome_ok=false");
            }
            if (!ok) {
                bool flip_x = false;
                bool flip_z = false;
                pp_log_line("[flip-PP] iter=" + std::to_string(iter + 1) + " start");
                std::vector<int> flip_x_union;
                std::vector<int> flip_z_union;
                std::vector<std::vector<int>> flip_x_window(flip_hist_x_window.begin(), flip_hist_x_window.end());
                std::vector<std::vector<int>> flip_z_window(flip_hist_z_window.begin(), flip_hist_z_window.end());
                if (flip_hist_active && hist_window > 0) {
                    std::vector<char> seen_x(nvars, 0);
                    std::vector<char> seen_z(nvars, 0);
                    for (const auto &vec : flip_hist_x_window) {
                        for (int v : vec) {
                            if (v >= 0 && v < nvars) seen_x[v] = 1;
                        }
                    }
                    for (const auto &vec : flip_hist_z_window) {
                        for (int v : vec) {
                            if (v >= 0 && v < nvars) seen_z[v] = 1;
                        }
                    }
                    flip_x_union.reserve(nvars / 8);
                    flip_z_union.reserve(nvars / 8);
                    for (int i = 0; i < nvars; ++i) {
                        if (seen_x[i]) flip_x_union.push_back(i);
                        if (seen_z[i]) flip_z_union.push_back(i);
                    }
                }
                auto x_sets = build_flip_sets(flip_x_union, flip_x_window, last_flip_x);
                for (const auto &set : x_sets) {
                    if (!set.vars || set.vars->empty()) continue;
                    if (attempt_flip_pp_side(
                            "X", *set.vars, x_checks, sx, true, trial_est, pp_ctx->verbose, nullptr, set.tag
                        )) {
                        flip_x = true;
                        break;
                    }
                }
                auto z_sets = build_flip_sets(flip_z_union, flip_z_window, last_flip_z);
                for (const auto &set : z_sets) {
                    if (!set.vars || set.vars->empty()) continue;
                    if (attempt_flip_pp_side(
                            "Z", *set.vars, z_checks, sz, false, trial_est, pp_ctx->verbose, nullptr, set.tag
                        )) {
                        flip_z = true;
                        break;
                    }
                }
                if (flip_x || flip_z) {
                    auto sx_pp = compute_syndrome(x_checks, trial_est, true);
                    auto sz_pp = compute_syndrome(z_checks, trial_est, false);
                    ok = (sx_pp == sx) && (sz_pp == sz);
                    if (ok) {
                        flip_pp_ok = true;
                    }
                }
                std::ostringstream oss;
                oss << "[flip-PP] iter=" << (iter + 1)
                    << " done x=" << tf(flip_x)
                    << " z=" << tf(flip_z)
                    << " syndrome_ok=" << tf(ok);
                pp_log_line(oss.str());
            }
            if (!ok) {
                auto sx_pp = compute_syndrome(x_checks, trial_est, true);
                auto sz_pp = compute_syndrome(z_checks, trial_est, false);
                int osd_us_x = unsatisfied_count(sx_pp, sx);
                int osd_us_z = unsatisfied_count(sz_pp, sz);
                int osd_us = std::max(osd_us_x, osd_us_z);
                if (osd_us > osd_us_limit) {
                    std::ostringstream oss;
                    oss << "[osd-PP] iter=" << (iter + 1)
                        << " skipped us=(" << osd_us_x << "," << osd_us_z << ")";
                    pp_log_line(oss.str());
                } else {
                    bool osd_x = false;
                    bool osd_z = false;
                    pp_log_line("[osd-PP] iter=" + std::to_string(iter + 1) + " start");
                    osd_x = attempt_osd_pp_side("X", abs_llr_x, x_checks, sx, true, trial_est,
                                                pp_ctx->verbose, truth_err);
                    osd_z = attempt_osd_pp_side("Z", abs_llr_z, z_checks, sz, false, trial_est,
                                                pp_ctx->verbose, truth_err);
                    if (osd_x || osd_z) {
                        auto sx_pp2 = compute_syndrome(x_checks, trial_est, true);
                        auto sz_pp2 = compute_syndrome(z_checks, trial_est, false);
                        ok = (sx_pp2 == sx) && (sz_pp2 == sz);
                    }
                    std::ostringstream oss;
                    oss << "[osd-PP] iter=" << (iter + 1)
                        << " done x=" << tf(osd_x)
                        << " z=" << tf(osd_z)
                        << " syndrome_ok=" << tf(ok);
                    pp_log_line(oss.str());
                }
            }
            if (ok) {
                est.swap(trial_est);
                pp_success = true;
                pp_success_ets = ets_pp_ok;
                pp_success_flip = flip_pp_ok;
                pp_used_ets_any = ets_used;
                pp_used_ets6 = ets6;
                pp_used_ets12 = ets12;
                if (verbose) {
                    auto sx_after = compute_syndrome(x_checks, est, true);
                    auto sz_after = compute_syndrome(z_checks, est, false);
                    StatusMetrics status_after;
                    compute_status_metrics(est, sx_after, sz_after, &trial_est, collect_llr, status_after);
                    update_status_min(status_after);
                    pp_log_line(format_status_en(iter + 1, status_after));
                }
                return make_result(iter + 1, false, true);
            }
        }
        prev_prev_est = prev_est;
        prev_est = est;
    }
    return make_result(max_iter, false, false);
}

static bool read_bit_vector(const std::string &path, std::vector<int> &out) {
    std::ifstream in(path);
    if (!in) return false;
    out.clear();
    int v = 0;
    while (in >> v) {
        out.push_back(v ? 1 : 0);
    }
    return true;
}

static bool read_int_vector(const std::string &path, std::vector<int> &out) {
    std::ifstream in(path);
    if (!in) return false;
    out.clear();
    int v = 0;
    while (in >> v) {
        out.push_back(v);
    }
    return true;
}

static bool write_error_vector(const std::string &path, const std::vector<int> &err) {
    std::ofstream ofs(path);
    if (!ofs) return false;
    for (int v : err) {
        ofs << v << "\n";
    }
    return true;
}

struct ETSEntry {
    std::vector<int> vars;
    std::array<int, 2> unsat_checks{};
};

struct ETSIndex {
    int P = 0;
    int var_count = 0;
    std::vector<ETSEntry> entries;
    std::unordered_map<uint64_t, std::vector<int>> by_pair;
    std::string source;
};

static std::string trim_ascii(const std::string &s) {
    size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

static char consume_side_prefix(std::string &line) {
    if (line.size() >= 2) {
        char ch = line[0];
        if ((ch == 'X' || ch == 'x' || ch == 'Z' || ch == 'z') &&
            std::isspace(static_cast<unsigned char>(line[1]))) {
            line = trim_ascii(line.substr(1));
            return static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        }
    }
    return '\0';
}

static std::vector<std::string> split_csv(const std::string &s) {
    std::vector<std::string> out;
    std::string cur;
    for (char ch : s) {
        if (ch == ',') {
            out.push_back(trim_ascii(cur));
            cur.clear();
        } else {
            cur.push_back(ch);
        }
    }
    if (!cur.empty()) out.push_back(trim_ascii(cur));
    return out;
}

static std::vector<std::string> split_ws(const std::string &s) {
    std::vector<std::string> out;
    std::string cur;
    for (char ch : s) {
        if (std::isspace(static_cast<unsigned char>(ch))) {
            if (!cur.empty()) {
                out.push_back(cur);
                cur.clear();
            }
        } else {
            cur.push_back(ch);
        }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

static bool parse_var_token(const std::string &tok, int P, int &var_id) {
    if (tok.size() < 3) return false;
    if (tok[0] != 'V' && tok[0] != 'v') return false;
    size_t pos = tok.find(':');
    if (pos == std::string::npos) return false;
    int col = std::stoi(tok.substr(1, pos - 1));
    int idx = std::stoi(tok.substr(pos + 1));
    var_id = col * P + idx;
    return true;
}

static bool parse_check_token(const std::string &tok, int P, int &check_id) {
    if (tok.size() < 3) return false;
    if (tok[0] != 'C' && tok[0] != 'c') return false;
    size_t pos = tok.find(':');
    if (pos == std::string::npos) return false;
    int row = std::stoi(tok.substr(1, pos - 1));
    int idx = std::stoi(tok.substr(pos + 1));
    check_id = row * P + idx;
    return true;
}

static bool read_ets_file(const std::string &path, int P, int expected_vars, char side_filter, ETSIndex &out) {
    std::ifstream in(path);
    if (!in) return false;
    out.entries.clear();
    out.by_pair.clear();
    out.var_count = expected_vars;
    out.P = P;
    out.source = path;
    if (expected_vars <= 0 || expected_vars > 32) return false;
    std::string line;
    while (std::getline(in, line)) {
        line = trim_ascii(line);
        if (line.empty()) continue;
        char line_side = consume_side_prefix(line);
        if ((side_filter == 'X' || side_filter == 'Z') &&
            (line_side != '\0' && line_side != side_filter)) {
            continue;
        }
        auto pos_vars = line.find("vars=");
        auto pos_unsat = line.find("unsat=");
        if (pos_vars == std::string::npos || pos_unsat == std::string::npos) continue;
        std::string vars_part = trim_ascii(line.substr(pos_vars + 5, pos_unsat - (pos_vars + 5)));
        std::string unsat_part = trim_ascii(line.substr(pos_unsat + 6));
        auto var_tokens = split_csv(vars_part);
        auto chk_tokens = split_csv(unsat_part);
        if (static_cast<int>(var_tokens.size()) != expected_vars || chk_tokens.size() != 2) continue;

        ETSEntry entry;
        entry.vars.resize(expected_vars);
        bool ok = true;
        for (size_t i = 0; i < var_tokens.size(); ++i) {
            int vid = -1;
            if (!parse_var_token(var_tokens[i], P, vid)) {
                ok = false;
                break;
            }
            entry.vars[i] = vid;
        }
        if (!ok) continue;
        for (size_t i = 0; i < chk_tokens.size(); ++i) {
            int cid = -1;
            if (!parse_check_token(chk_tokens[i], P, cid)) {
                ok = false;
                break;
            }
            entry.unsat_checks[i] = cid;
        }
        if (!ok) continue;
        std::sort(entry.vars.begin(), entry.vars.end());
        std::sort(entry.unsat_checks.begin(), entry.unsat_checks.end());
        int c0 = entry.unsat_checks[0];
        int c1 = entry.unsat_checks[1];
        uint64_t key = (static_cast<uint64_t>(static_cast<uint32_t>(c0)) << 32) |
                       static_cast<uint32_t>(c1);
        int idx = static_cast<int>(out.entries.size());
        out.entries.push_back(entry);
        out.by_pair[key].push_back(idx);
    }
    return true;
}

static void decode_var(int v, int P, int &col, int &idx) {
    col = v / P;
    idx = v % P;
}

static void decode_check(int c, int P, int &row, int &idx) {
    row = c / P;
    idx = c % P;
}

static bool read_cycles_file(const std::string &path,
                             int P,
                             int var_count,
                             std::vector<CycleEntry> &out_x,
                             std::vector<CycleEntry> &out_z) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    while (std::getline(in, line)) {
        line = trim_ascii(line);
        if (line.empty()) continue;
        char side = consume_side_prefix(line);
        auto tokens = split_ws(line);
        if (tokens.empty()) continue;
        size_t start = 0;
        if (!tokens.empty()) {
            char t0 = tokens[0].empty() ? '\0' : tokens[0][0];
            if (t0 != 'V' && t0 != 'v' && t0 != 'C' && t0 != 'c') {
                start = 1;
            }
        }
        std::vector<int> vars;
        std::vector<int> checks;
        vars.reserve(static_cast<size_t>(var_count));
        checks.reserve(static_cast<size_t>(var_count));
        for (size_t i = start; i < tokens.size(); ++i) {
            int id = -1;
            if (parse_var_token(tokens[i], P, id)) {
                vars.push_back(id);
                continue;
            }
            if (parse_check_token(tokens[i], P, id)) {
                checks.push_back(id);
                continue;
            }
        }
        if (static_cast<int>(vars.size()) != var_count ||
            static_cast<int>(checks.size()) != var_count) {
            continue;
        }
        CycleEntry entry;
        for (int i = 0; i < var_count; ++i) {
            entry.vars[i] = vars[i];
            entry.checks[i] = checks[i];
        }
        if (side == 'X') {
            out_x.push_back(entry);
        } else if (side == 'Z') {
            out_z.push_back(entry);
        } else {
            out_x.push_back(entry);
            out_z.push_back(entry);
        }
    }
    return true;
}

static std::string find_cycles_path(const std::string &log_dir,
                                    const std::vector<std::string> &names) {
    for (const auto &name : names) {
        if (!log_dir.empty()) {
            std::string path = log_dir + "/" + name;
            if (std::filesystem::exists(path)) return path;
        }
        if (std::filesystem::exists(name)) return name;
    }
    return "";
}

static void load_cycles_index(const std::string &log_dir, int P, CycleIndex &out) {
    std::string cycles8_path = find_cycles_path(log_dir, {"cycles8.txt"});
    if (!cycles8_path.empty()) {
        if (read_cycles_file(cycles8_path, P, 4, out.x.cycles8, out.z.cycles8)) {
            out.x.have8 = !out.x.cycles8.empty();
            out.z.have8 = !out.z.cycles8.empty();
        }
    }
    std::string cycle10_path = find_cycles_path(log_dir, {"cycle10.txt", "cycles10.txt"});
    if (!cycle10_path.empty()) {
        if (read_cycles_file(cycle10_path, P, 5, out.x.cycles10, out.z.cycles10)) {
            out.x.have10 = !out.x.cycles10.empty();
            out.z.have10 = !out.z.cycles10.empty();
        }
    }
    out.any_loaded = out.x.have8 || out.x.have10 || out.z.have8 || out.z.have10;
}

static bool solve_gf2(
    const std::vector<uint32_t> &rows,
    const std::vector<int> &rhs,
    int nvars,
    std::vector<int> &solution
) {
    int m = static_cast<int>(rows.size());
    std::vector<uint32_t> a = rows;
    std::vector<int> b = rhs;
    std::vector<int> where(nvars, -1);
    int r = 0;
    for (int c = 0; c < nvars && r < m; ++c) {
        int sel = r;
        while (sel < m && (((a[sel] >> c) & 1u) == 0u)) ++sel;
        if (sel == m) continue;
        std::swap(a[sel], a[r]);
        std::swap(b[sel], b[r]);
        where[c] = r;
        for (int i = 0; i < m; ++i) {
            if (i != r && (((a[i] >> c) & 1u) != 0u)) {
                a[i] ^= a[r];
                b[i] ^= b[r];
            }
        }
        ++r;
    }
    for (int i = r; i < m; ++i) {
        if (a[i] == 0u && b[i] != 0) return false;
    }
    solution.assign(nvars, 0);
    for (int c = 0; c < nvars; ++c) {
        if (where[c] != -1) solution[c] = b[where[c]] & 1;
    }
    return true;
}

struct GF2SolveInfo {
    bool solvable = false;
    int rank = 0;
    bool unique = false;
};

struct GF2BasisResult {
    bool solvable = false;
    int rank = 0;
    bool unique = false;
    std::vector<int> where;
    std::vector<std::vector<uint64_t>> rows;
    std::vector<int> rhs;
};

static GF2SolveInfo solve_gf2_bits(
    std::vector<std::vector<uint64_t>> rows,
    std::vector<int> rhs,
    int nvars,
    std::vector<int> &solution
) {
    GF2SolveInfo info;
    int m = static_cast<int>(rows.size());
    if (nvars <= 0) return info;
    int words = (nvars + 63) / 64;
    std::vector<int> where(nvars, -1);
    int r = 0;
    for (int c = 0; c < nvars && r < m; ++c) {
        int w = c >> 6;
        uint64_t mask = 1ULL << (c & 63);
        int sel = r;
        while (sel < m && (rows[sel][w] & mask) == 0) ++sel;
        if (sel == m) continue;
        std::swap(rows[sel], rows[r]);
        std::swap(rhs[sel], rhs[r]);
        where[c] = r;
        for (int i = 0; i < m; ++i) {
            if (i == r) continue;
            if ((rows[i][w] & mask) == 0) continue;
            for (int k = 0; k < words; ++k) {
                rows[i][k] ^= rows[r][k];
            }
            rhs[i] ^= rhs[r];
        }
        ++r;
    }
    for (int i = r; i < m; ++i) {
        bool zero = true;
        for (int k = 0; k < words; ++k) {
            if (rows[i][k] != 0) {
                zero = false;
                break;
            }
        }
        if (zero && (rhs[i] & 1)) {
            info.solvable = false;
            info.rank = r;
            info.unique = false;
            return info;
        }
    }
    info.solvable = true;
    info.rank = r;
    info.unique = (r == nvars);
    solution.assign(nvars, 0);
    for (int c = 0; c < nvars; ++c) {
        if (where[c] != -1) solution[c] = rhs[where[c]] & 1;
    }
    return info;
}

static GF2BasisResult solve_gf2_bits_basis(
    std::vector<std::vector<uint64_t>> rows,
    std::vector<int> rhs,
    int nvars
) {
    GF2BasisResult out;
    int m = static_cast<int>(rows.size());
    if (nvars <= 0) return out;
    int words = (nvars + 63) / 64;
    std::vector<int> where(nvars, -1);
    int r = 0;
    for (int c = 0; c < nvars && r < m; ++c) {
        int w = c >> 6;
        uint64_t mask = 1ULL << (c & 63);
        int sel = r;
        while (sel < m && (rows[sel][w] & mask) == 0) ++sel;
        if (sel == m) continue;
        std::swap(rows[sel], rows[r]);
        std::swap(rhs[sel], rhs[r]);
        where[c] = r;
        for (int i = 0; i < m; ++i) {
            if (i == r) continue;
            if ((rows[i][w] & mask) == 0) continue;
            for (int k = 0; k < words; ++k) {
                rows[i][k] ^= rows[r][k];
            }
            rhs[i] ^= rhs[r];
        }
        ++r;
    }
    for (int i = r; i < m; ++i) {
        bool zero = true;
        for (int k = 0; k < words; ++k) {
            if (rows[i][k] != 0) {
                zero = false;
                break;
            }
        }
        if (zero && (rhs[i] & 1)) {
            out.solvable = false;
            out.rank = r;
            out.unique = false;
            return out;
        }
    }
    out.solvable = true;
    out.rank = r;
    out.unique = (r == nvars);
    out.where = std::move(where);
    out.rows = std::move(rows);
    out.rhs = std::move(rhs);
    return out;
}

static int find_local_index(const std::vector<int> &vars, int v) {
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i] == v) return static_cast<int>(i);
    }
    return -1;
}

static bool attempt_ets_fix_side(
    const std::string &label,
    const ETSIndex &ets_index,
    const std::vector<std::vector<int>> &checks,
    const std::vector<std::vector<EdgeRef>> &var_to_checks,
    const std::vector<int> &target_syndrome,
    bool use_xbit,
    std::vector<int> &est,
    bool verbose
) {
    if (ets_index.entries.empty()) return false;
    int nvars = ets_index.var_count;
    if (nvars <= 0 || nvars > 32) return false;
    int P = ets_index.P;
    auto synd_hat = compute_syndrome(checks, est, use_xbit);
    std::vector<int> diff(synd_hat.size(), 0);
    std::vector<int> unsat;
    for (size_t i = 0; i < synd_hat.size(); ++i) {
        int d = synd_hat[i] ^ target_syndrome[i];
        diff[i] = d;
        if (d) unsat.push_back(static_cast<int>(i));
    }
    if (unsat.size() != 2) return false;
    int c0 = std::min(unsat[0], unsat[1]);
    int c1 = std::max(unsat[0], unsat[1]);
    uint64_t key = (static_cast<uint64_t>(static_cast<uint32_t>(c0)) << 32) |
                   static_cast<uint32_t>(c1);
    auto it = ets_index.by_pair.find(key);
    if (it == ets_index.by_pair.end()) return false;

    for (int idx : it->second) {
        const auto &ets = ets_index.entries[idx];
        if (static_cast<int>(ets.vars.size()) != nvars) continue;
        std::vector<int> check_list;
        check_list.reserve(3 * nvars);
        for (int v : ets.vars) {
            for (const auto &e : var_to_checks[v]) {
                bool seen = false;
                for (int c : check_list) {
                    if (c == e.check) {
                        seen = true;
                        break;
                    }
                }
                if (!seen) check_list.push_back(e.check);
            }
        }

        bool ok = true;
        for (int c : check_list) {
            if (diff[c] == 1 && c != c0 && c != c1) {
                ok = false;
                break;
            }
        }
        if (!ok) continue;

        std::vector<uint32_t> rows;
        std::vector<int> rhs;
        rows.reserve(check_list.size());
        rhs.reserve(check_list.size());

        for (int c : check_list) {
            uint32_t mask = 0;
            for (int v : checks[c]) {
                if (v < 0) continue;
                int local = find_local_index(ets.vars, v);
                if (local >= 0) {
                    mask |= (1u << local);
                }
            }
            if (mask == 0u) continue;
            rows.push_back(mask);
            rhs.push_back(diff[c] & 1);
        }

        std::vector<int> delta;
        if (!solve_gf2(rows, rhs, nvars, delta)) continue;

        std::vector<int> trial = est;
        for (int i = 0; i < nvars; ++i) {
            if (delta[i] == 0) continue;
            int v = ets.vars[i];
            trial[v] = use_xbit ? (trial[v] ^ 1) : (trial[v] ^ 2);
        }

        auto new_synd = compute_syndrome(checks, trial, use_xbit);
        if (new_synd == target_syndrome) {
            est.swap(trial);
            if (verbose) {
                std::ostringstream oss;
                oss << "[ETS-" << label << "] applied correction using ETS index " << idx;
                if (!ets_index.source.empty()) {
                    oss << " source=" << ets_index.source;
                }
                if (P > 0) {
                    oss << " vars=";
                    for (size_t i = 0; i < ets.vars.size(); ++i) {
                        int col = 0, vidx = 0;
                        decode_var(ets.vars[i], P, col, vidx);
                        if (i) oss << ",";
                        oss << "V" << col << ":" << vidx;
                    }
                    oss << " unsat=";
                    for (size_t i = 0; i < ets.unsat_checks.size(); ++i) {
                        int row = 0, cidx = 0;
                        decode_check(ets.unsat_checks[i], P, row, cidx);
                        if (i) oss << ",";
                        oss << "C" << row << ":" << cidx;
                    }
                }
                std::cout << oss.str() << "\n";
            }
            return true;
        }
    }
    return false;
}

static ETSFixResult attempt_ets_fix_side_multi(
    const std::string &label,
    const std::vector<ETSIndexRef> &ets_list,
    const std::vector<std::vector<int>> &checks,
    const std::vector<std::vector<EdgeRef>> &var_to_checks,
    const std::vector<int> &target_syndrome,
    bool use_xbit,
    std::vector<int> &est,
    bool verbose
) {
    ETSFixResult result;
    for (const auto &entry : ets_list) {
        const ETSIndex *ets = use_xbit ? entry.x : entry.z;
        if (!ets) continue;
        std::string tag = label + (entry.label.empty() ? std::to_string(entry.size) : entry.label);
        if (attempt_ets_fix_side(tag, *ets, checks, var_to_checks, target_syndrome,
                                 use_xbit, est, verbose)) {
            result.applied = true;
            result.size = entry.size;
            result.label = entry.label;
            break;
        }
    }
    return result;
}

static std::vector<FlipSetRef> build_flip_sets(
    const std::vector<int> &union_vars,
    const std::vector<std::vector<int>> &window,
    const std::vector<int> &last_vars
) {
    std::vector<FlipSetRef> sets;
    if (!union_vars.empty()) {
        sets.push_back({"union", &union_vars});
    }
    if (!window.empty()) {
        int idx = 0;
        for (int i = static_cast<int>(window.size()) - 1; i >= 0; --i, ++idx) {
            if (window[i].empty()) continue;
            sets.push_back({std::string("hist") + std::to_string(idx), &window[i]});
        }
    } else if (!last_vars.empty() && union_vars.empty()) {
        sets.push_back({"last", &last_vars});
    }
    return sets;
}

static bool apply_flip_pp_best(
    const std::string &label,
    const FlipPPBest &best,
    const std::vector<std::vector<int>> &checks,
    const std::vector<int> &target_syndrome,
    bool use_xbit,
    std::vector<int> &est,
    bool verbose
) {
    if (!best.has) return false;
    if (best.vars.size() != best.delta.size()) return false;
    std::vector<int> trial = est;
    for (size_t i = 0; i < best.vars.size(); ++i) {
        if ((best.delta[i] & 1) == 0) continue;
        int v = best.vars[i];
        trial[v] = use_xbit ? (trial[v] ^ 1) : (trial[v] ^ 2);
    }
    auto new_synd = compute_syndrome(checks, trial, use_xbit);
    bool ok = (new_synd == target_syndrome);
    if (verbose) {
        std::cout << "[flip-PP-" << label << "] fallback_min_weight=" << best.weight;
        if (!best.tag.empty()) {
            std::cout << " set=" << best.tag;
        }
        std::cout << " syndrome_ok=" << tf(ok) << "\n";
    }
    if (ok) {
        est.swap(trial);
        if (verbose) {
            std::cout << "[flip-PP-" << label << "] applied fallback using " << best.vars.size() << " flips\n";
        }
    }
    return ok;
}

static bool attempt_flip_pp_side(
    const std::string &label,
    const std::vector<int> &flip_vars_in,
    const std::vector<std::vector<int>> &checks,
    const std::vector<int> &target_syndrome,
    bool use_xbit,
    std::vector<int> &est,
    bool verbose,
    FlipPPBest *best,
    const std::string &set_tag
) {
    int flip_in = static_cast<int>(flip_vars_in.size());
    std::string tag = set_tag.empty() ? "" : (std::string(" set=") + set_tag);
    if (flip_vars_in.empty()) {
        if (verbose) {
            std::cout << "[flip-PP-" << label << "]" << tag << " flip_count=0 eq=H[N(flip),flip] rows=0 cols=0"
                      << " solvable=false dof=NA unique=false (no flips)\n";
        }
        return false;
    }
    int nvars_total = static_cast<int>(est.size());
    std::vector<int> vars;
    vars.reserve(flip_vars_in.size());
    for (int v : flip_vars_in) {
        if (v >= 0 && v < nvars_total) vars.push_back(v);
    }
    if (vars.empty()) {
        if (verbose) {
            std::cout << "[flip-PP-" << label << "]" << tag << " flip_count_in=" << flip_in
                      << " flip_count_used=0 eq=H[N(flip),flip] rows=0 cols=0"
                      << " solvable=false dof=NA unique=false (no valid flips)\n";
        }
        return false;
    }
    vars = sorted_unique(std::move(vars));
    int nvars = static_cast<int>(vars.size());
    if (nvars == 0) return false;
    if (nvars != flip_in) {
        if (verbose) {
            std::cout << "[flip-PP-" << label << "]" << tag << " flip_count_in=" << flip_in
                      << " flip_count_used=" << nvars << "\n";
        }
    }

    auto synd_hat = compute_syndrome(checks, est, use_xbit);
    std::vector<int> diff(synd_hat.size(), 0);
    for (size_t i = 0; i < synd_hat.size(); ++i) {
        diff[i] = synd_hat[i] ^ target_syndrome[i];
    }

    std::vector<int> pos(nvars_total, -1);
    for (int i = 0; i < nvars; ++i) {
        pos[vars[i]] = i;
    }

    int words = (nvars + 63) / 64;
    std::vector<std::vector<uint64_t>> rows;
    std::vector<int> rhs;
    rows.reserve(checks.size());
    rhs.reserve(checks.size());
    int rows_used = 0;
    for (size_t c = 0; c < checks.size(); ++c) {
        std::vector<uint64_t> row(words, 0);
        for (int v : checks[c]) {
            if (v < 0 || v >= nvars_total) continue;
            int p = pos[v];
            if (p < 0) continue;
            row[p >> 6] |= (1ULL << (p & 63));
        }
        bool empty = true;
        for (int k = 0; k < words; ++k) {
            if (row[k] != 0) {
                empty = false;
                break;
            }
        }
        if (empty) {
            if (diff[c] != 0) {
                if (verbose) {
                    std::cout << "[flip-PP-" << label << "]" << tag << " flip_count=" << nvars
                              << " eq=H[N(flip),flip] rows=0 cols=" << nvars
                              << " solvable=false dof=NA unique=false\n";
                }
                return false;
            }
            continue;
        }
        rows_used++;
        rows.push_back(std::move(row));
        rhs.push_back(diff[c] & 1);
    }

    if (rows.empty()) {
        int dof = nvars;
        if (verbose) {
            std::cout << "[flip-PP-" << label << "]" << tag << " flip_count=" << nvars
                      << " eq=H[N(flip),flip] rows=0 cols=" << nvars
                      << " solvable=true dof=" << dof << " unique=false\n";
        }
        return false;
    }

    std::vector<int> delta;
    GF2SolveInfo info = solve_gf2_bits(std::move(rows), std::move(rhs), nvars, delta);
    if (!info.solvable) {
        if (verbose) {
            std::cout << "[flip-PP-" << label << "]" << tag << " flip_count=" << nvars
                      << " eq=H[N(flip),flip] rows=" << rows_used
                      << " cols=" << nvars
                      << " solvable=false dof=NA unique=false\n";
        }
        return false;
    }
    int dof = nvars - info.rank;
    if (verbose) {
        std::cout << "[flip-PP-" << label << "]" << tag << " flip_count=" << nvars
                  << " eq=H[N(flip),flip] rows=" << rows_used
                  << " cols=" << nvars
                  << " solvable=true dof=" << dof
                  << " unique=" << tf(info.unique) << "\n";
    }

    std::vector<int> trial = est;
    for (int i = 0; i < nvars; ++i) {
        if (delta[i] == 0) continue;
        int v = vars[i];
        trial[v] = use_xbit ? (trial[v] ^ 1) : (trial[v] ^ 2);
    }
    auto new_synd = compute_syndrome(checks, trial, use_xbit);
    bool match = (new_synd == target_syndrome);
    if (info.unique && match) {
        est.swap(trial);
        if (verbose) {
            std::cout << "[flip-PP-" << label << "] applied using " << nvars << " flips\n";
        }
        return true;
    }
    if (match && best) {
        int weight = 0;
        for (int v : delta) {
            if (v & 1) weight++;
        }
        if (!best->has || weight < best->weight) {
            best->has = true;
            best->weight = weight;
            best->vars = vars;
            best->delta = delta;
            best->tag = set_tag;
        }
    }
    if (verbose && !match) {
        std::cout << "[flip-PP-" << label << "] solved but syndrome mismatch\n";
    }
    return false;
}

static OSDPPInfo solve_osd_pp_for_k(
    const std::vector<int> &order,
    int k,
    const std::vector<std::vector<int>> &checks,
    const std::vector<int> &diff,
    int nvars_total
) {
    OSDPPInfo info;
    info.cols = k;
    if (k <= 0 || static_cast<int>(order.size()) < k) {
        return info;
    }
    std::vector<int> pos(nvars_total, -1);
    for (int i = 0; i < k; ++i) {
        int v = order[i];
        if (v >= 0 && v < nvars_total) {
            pos[v] = i;
        }
    }
    int words = (k + 63) / 64;
    std::vector<std::vector<uint64_t>> rows;
    std::vector<int> rhs;
    rows.reserve(checks.size());
    rhs.reserve(checks.size());
    int rows_used = 0;
    for (size_t c = 0; c < checks.size(); ++c) {
        std::vector<uint64_t> row(words, 0);
        for (int v : checks[c]) {
            if (v < 0 || v >= nvars_total) continue;
            int p = pos[v];
            if (p < 0) continue;
            row[p >> 6] |= (1ULL << (p & 63));
        }
        bool empty = true;
        for (int w = 0; w < words; ++w) {
            if (row[w] != 0) {
                empty = false;
                break;
            }
        }
        if (empty) {
            if (diff[c] != 0) {
                info.solvable = false;
                info.rows = rows_used;
                return info;
            }
            continue;
        }
        rows_used++;
        rows.push_back(std::move(row));
        rhs.push_back(diff[c] & 1);
    }
    info.rows = rows_used;
    if (rows.empty()) {
        info.solvable = true;
        info.rank = 0;
        info.unique = (k == 0);
        info.delta.assign(k, 0);
        return info;
    }
    std::vector<int> delta;
    GF2SolveInfo gf = solve_gf2_bits(std::move(rows), std::move(rhs), k, delta);
    info.solvable = gf.solvable;
    info.rank = gf.rank;
    info.unique = gf.unique;
    info.delta = std::move(delta);
    return info;
}

static bool osd_min_weight_for_k(
    const std::vector<int> &order,
    int k,
    const std::vector<std::vector<int>> &checks,
    const std::vector<int> &diff,
    int nvars_total,
    OSDMinWeightCandidate &out
) {
    if (k <= 0 || static_cast<int>(order.size()) < k) {
        return false;
    }
    std::vector<int> pos(nvars_total, -1);
    for (int i = 0; i < k; ++i) {
        int v = order[i];
        if (v >= 0 && v < nvars_total) {
            pos[v] = i;
        }
    }
    int words = (k + 63) / 64;
    std::vector<std::vector<uint64_t>> rows;
    std::vector<int> rhs;
    rows.reserve(checks.size());
    rhs.reserve(checks.size());
    int rows_used = 0;
    for (size_t c = 0; c < checks.size(); ++c) {
        std::vector<uint64_t> row(words, 0);
        for (int v : checks[c]) {
            if (v < 0 || v >= nvars_total) continue;
            int p = pos[v];
            if (p < 0) continue;
            row[p >> 6] |= (1ULL << (p & 63));
        }
        bool empty = true;
        for (int w = 0; w < words; ++w) {
            if (row[w] != 0) {
                empty = false;
                break;
            }
        }
        if (empty) {
            if (diff[c] != 0) {
                return false;
            }
            continue;
        }
        rows_used++;
        rows.push_back(std::move(row));
        rhs.push_back(diff[c] & 1);
    }
    if (rows.empty()) {
        return false;
    }
    GF2BasisResult basis = solve_gf2_bits_basis(std::move(rows), std::move(rhs), k);
    if (!basis.solvable) {
        return false;
    }
    int dof = k - basis.rank;
    if (dof > 5) {
        return false;
    }

    std::vector<int> free_vars;
    free_vars.reserve(dof);
    for (int c = 0; c < k; ++c) {
        if (basis.where[c] == -1) free_vars.push_back(c);
    }
    if (static_cast<int>(free_vars.size()) != dof) {
        return false;
    }

    std::vector<uint64_t> particular_bits(words, 0);
    for (int c = 0; c < k; ++c) {
        int r = basis.where[c];
        if (r >= 0 && (basis.rhs[r] & 1)) {
            particular_bits[c >> 6] |= (1ULL << (c & 63));
        }
    }

    std::vector<std::vector<uint64_t>> basis_bits;
    basis_bits.reserve(dof);
    for (int fv : free_vars) {
        std::vector<uint64_t> vec(words, 0);
        vec[fv >> 6] |= (1ULL << (fv & 63));
        for (int c = 0; c < k; ++c) {
            int r = basis.where[c];
            if (r < 0) continue;
            if (basis.rows[r][fv >> 6] & (1ULL << (fv & 63))) {
                vec[c >> 6] |= (1ULL << (c & 63));
            }
        }
        basis_bits.push_back(std::move(vec));
    }

    int combos = 1 << dof;
    int best_weight = std::numeric_limits<int>::max();
    std::vector<uint64_t> best_bits;
    for (int mask = 0; mask < combos; ++mask) {
        std::vector<uint64_t> cand = particular_bits;
        for (int i = 0; i < dof; ++i) {
            if (mask & (1 << i)) {
                for (int w = 0; w < words; ++w) {
                    cand[w] ^= basis_bits[i][w];
                }
            }
        }
        int weight = 0;
        for (int w = 0; w < words; ++w) {
            weight += __builtin_popcountll(cand[w]);
        }
        if (weight < best_weight) {
            best_weight = weight;
            best_bits = std::move(cand);
        }
    }
    if (best_bits.empty()) {
        return false;
    }
    out.has = true;
    out.k = k;
    out.dof = dof;
    out.rows = rows_used;
    out.cols = k;
    out.weight = best_weight;
    out.delta_bits = std::move(best_bits);
    return true;
}

static bool attempt_osd_pp_side(
    const std::string &label,
    const std::vector<double> &abs_llr,
    const std::vector<std::vector<int>> &checks,
    const std::vector<int> &target_syndrome,
    bool use_xbit,
    std::vector<int> &est,
    bool verbose,
    const std::vector<int> *truth_err
) {
    int nvars_total = static_cast<int>(est.size());
    if (static_cast<int>(abs_llr.size()) != nvars_total) {
        if (verbose) {
            std::cout << "[osd-PP-" << label << "] skipped (llr unavailable)\n";
        }
        return false;
    }
    if (checks.empty()) {
        if (verbose) {
            std::cout << "[osd-PP-" << label << "] skipped (no checks)\n";
        }
        return false;
    }
    auto synd_hat = compute_syndrome(checks, est, use_xbit);
    std::vector<int> diff(synd_hat.size(), 0);
    bool any_diff = false;
    for (size_t i = 0; i < synd_hat.size(); ++i) {
        diff[i] = synd_hat[i] ^ target_syndrome[i];
        if (diff[i]) any_diff = true;
    }
    if (!any_diff) {
        if (verbose) {
            std::cout << "[osd-PP-" << label << "] skipped (no syndrome diff)\n";
        }
        return false;
    }

    std::vector<int> order(nvars_total, 0);
    for (int i = 0; i < nvars_total; ++i) {
        order[i] = i;
    }
    std::stable_sort(order.begin(), order.end(),
                     [&](int a, int b) {
                         if (abs_llr[a] == abs_llr[b]) return a < b;
                         return abs_llr[a] < abs_llr[b];
                     });

    bool have_truth = (truth_err != nullptr && truth_err->size() == est.size());
    std::vector<uint8_t> true_diff;
    if (have_truth) {
        true_diff.assign(nvars_total, 0);
        for (int i = 0; i < nvars_total; ++i) {
            int t = use_xbit ? xbit((*truth_err)[i]) : zbit((*truth_err)[i]);
            int e = use_xbit ? xbit(est[i]) : zbit(est[i]);
            if ((t ^ e) & 1) {
                true_diff[i] = 1;
            }
        }
    }
    auto true_diff_ok_for_k = [&](int k) -> bool {
        if (k <= 0) {
            for (int d : diff) {
                if ((d & 1) != 0) return false;
            }
            return true;
        }
        std::vector<char> in_subset(nvars_total, 0);
        for (int i = 0; i < k; ++i) {
            int v = order[i];
            if (v >= 0 && v < nvars_total) {
                in_subset[v] = 1;
            }
        }
        for (size_t c = 0; c < checks.size(); ++c) {
            int parity = 0;
            for (int v : checks[c]) {
                if (v < 0 || v >= nvars_total) continue;
                if (!in_subset[v]) continue;
                parity ^= (true_diff[v] & 1);
            }
            if (parity != (diff[c] & 1)) return false;
        }
        return true;
    };

    int k_max = nvars_total;
    if (k_max <= 0) {
        if (verbose) {
            std::cout << "[osd-PP-" << label << "] skipped (K_max=0)\n";
        }
        return false;
    }
    if (verbose) {
        std::cout << "[osd-PP-" << label << "] K_max=" << k_max << "\n";
    }

    std::unordered_map<int, OSDPPInfo> cache;
    auto get_info = [&](int k) -> OSDPPInfo & {
        auto &info = cache[k];
        if (!info.computed) {
            info = solve_osd_pp_for_k(order, k, checks, diff, nvars_total);
            info.computed = true;
            if (verbose) {
                if (info.solvable) {
                    int dof = info.cols - info.rank;
                    std::cout << "[osd-PP-" << label << "] K_test=" << k
                              << " rows=" << info.rows
                              << " cols=" << info.cols
                              << " solvable=true dof=" << dof
                              << " unique=" << tf(info.unique);
                    if (have_truth) {
                        std::cout << " true_diff_ok=" << tf(true_diff_ok_for_k(k));
                    } else {
                        std::cout << " true_diff_ok=NA";
                    }
                    std::cout << "\n";
                } else {
                    std::cout << "[osd-PP-" << label << "] K_test=" << k
                              << " rows=" << info.rows
                              << " cols=" << info.cols
                              << " solvable=false dof=NA"
                              << " unique=false";
                    if (have_truth) {
                        std::cout << " true_diff_ok=" << tf(true_diff_ok_for_k(k));
                    } else {
                        std::cout << " true_diff_ok=NA";
                    }
                    std::cout << "\n";
                }
            }
        }
        return info;
    };

    auto &info_max = get_info(k_max);
    if (!info_max.solvable) {
        if (verbose) {
            std::cout << "[osd-PP-" << label << "] K_max=" << k_max
                      << " solvable=false unique=false\n";
        }
        return false;
    }

    int lo = 1;
    int hi = k_max;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        bool solvable = get_info(mid).solvable;
        if (verbose) {
            std::cout << "[osd-PP-" << label << "] K_search_solve mid=" << mid
                      << " solvable=" << tf(solvable) << "\n";
        }
        if (solvable) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    int k_solve = lo;
    auto &info_solve = get_info(k_solve);
    int k_use = -1;
    if (!info_solve.unique) {
        int dof = info_solve.solvable ? (info_solve.cols - info_solve.rank) : -1;
        if (verbose) {
            std::cout << "[osd-PP-" << label << "] K_solve=" << k_solve
                      << " solvable=" << tf(info_solve.solvable)
                      << " dof=" << dof
                      << " unique=false\n";
        }
        return false;
    } else {
        int lo2 = k_solve;
        int hi2 = k_max;
        while (lo2 < hi2) {
            int mid = lo2 + (hi2 - lo2 + 1) / 2;
            bool unique = get_info(mid).unique;
            if (verbose) {
                std::cout << "[osd-PP-" << label << "] K_search_unique mid=" << mid
                          << " unique=" << tf(unique) << "\n";
            }
            if (unique) {
                lo2 = mid;
            } else {
                hi2 = mid - 1;
            }
        }
        k_use = lo2;
        auto &info_use = get_info(k_use);
        if (info_use.solvable && info_use.unique) {
            const int osd_weight_limit = 30;
            std::vector<int> trial = est;
            int weight = 0;
            for (int i = 0; i < k_use; ++i) {
                if (i >= static_cast<int>(info_use.delta.size())) break;
                if ((info_use.delta[i] & 1) == 0) continue;
                weight++;
                int v = order[i];
                trial[v] = use_xbit ? (trial[v] ^ 1) : (trial[v] ^ 2);
            }
            auto new_synd = compute_syndrome(checks, trial, use_xbit);
            bool ok = (new_synd == target_syndrome);
            int dof = info_use.cols - info_use.rank;
            bool apply = ok && (weight <= osd_weight_limit);
            if (verbose) {
                std::cout << "[osd-PP-" << label << "] K=" << k_use
                          << " rows=" << info_use.rows
                          << " cols=" << info_use.cols
                          << " solvable=" << tf(info_use.solvable)
                          << " dof=" << dof
                          << " unique=" << tf(info_use.unique)
                          << " weight=" << weight
                          << " syndrome_ok=" << tf(ok)
                          << " apply=" << tf(apply)
                          << "\n";
            }
            if (apply) {
                est.swap(trial);
                if (verbose) {
                    std::cout << "[osd-PP-" << label << "] applied using K=" << k_use << "\n";
                }
                return true;
            }
            return false;
        } else {
            if (verbose) {
                std::cout << "[osd-PP-" << label << "] K=" << k_use
                          << " solvable=" << tf(info_use.solvable)
                          << " unique=" << tf(info_use.unique) << "\n";
            }
        }
    }

    int min_dof = std::numeric_limits<int>::max();
    std::vector<int> min_dof_ks;
    for (const auto &kv : cache) {
        const auto &info = kv.second;
        if (!info.computed || !info.solvable) continue;
        int dof = info.cols - info.rank;
        if (dof < min_dof) {
            min_dof = dof;
            min_dof_ks.clear();
            min_dof_ks.push_back(kv.first);
        } else if (dof == min_dof) {
            min_dof_ks.push_back(kv.first);
        }
    }
    if (min_dof == std::numeric_limits<int>::max()) {
        if (verbose) {
            std::cout << "[osd-PP-" << label << "] fallback_min_dof=NA\n";
        }
        return false;
    }
    if (verbose) {
        std::cout << "[osd-PP-" << label << "] fallback_min_dof=" << min_dof
                  << " candidates=" << min_dof_ks.size() << "\n";
    }
    if (min_dof > 5) {
        return false;
    }

    OSDMinWeightCandidate best;
    for (int k : min_dof_ks) {
        OSDMinWeightCandidate cand;
        if (!osd_min_weight_for_k(order, k, checks, diff, nvars_total, cand)) {
            continue;
        }
        if (!best.has || cand.weight < best.weight ||
            (cand.weight == best.weight && cand.k < best.k)) {
            best = std::move(cand);
        }
    }
    if (!best.has) {
        if (verbose) {
            std::cout << "[osd-PP-" << label << "] fallback_min_weight=NA\n";
        }
        return false;
    }
    if (best.weight > 20) {
        if (verbose) {
            std::cout << "[osd-PP-" << label << "] fallback_min_weight=" << best.weight
                      << " threshold=20 skip\n";
        }
        return false;
    }

    std::vector<int> trial = est;
    for (int i = 0; i < best.k; ++i) {
        if (best.delta_bits[i >> 6] & (1ULL << (i & 63))) {
            int v = order[i];
            trial[v] = use_xbit ? (trial[v] ^ 1) : (trial[v] ^ 2);
        }
    }
    auto new_synd = compute_syndrome(checks, trial, use_xbit);
    bool ok = (new_synd == target_syndrome);
    if (verbose) {
        std::cout << "[osd-PP-" << label << "] fallback_min_weight=" << best.weight
                  << " K=" << best.k
                  << " dof=" << best.dof
                  << " syndrome_ok=" << tf(ok)
                  << "\n";
    }
    if (ok) {
        est.swap(trial);
        if (verbose) {
            std::cout << "[osd-PP-" << label << "] applied fallback using K=" << best.k << "\n";
        }
    }
    return ok;
}

static bool wilson_interval(long long n, long long k, double z, double &low, double &high) {
    if (n <= 0) return false;
    double nn = static_cast<double>(n);
    double phat = static_cast<double>(k) / nn;
    double z2 = z * z;
    double denom = 1.0 + z2 / nn;
    double center = (phat + z2 / (2.0 * nn)) / denom;
    double half = z * std::sqrt((phat * (1.0 - phat) + z2 / (4.0 * nn)) / nn) / denom;
    low = std::max(0.0, center - half);
    high = std::min(1.0, center + half);
    return true;
}

static uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31);
}

static uint64_t trial_seed_from_base(uint64_t base_seed, uint64_t trial_idx) {
    uint64_t mix = base_seed ^ (trial_idx * 0x9e3779b97f4a7c15ull);
    return splitmix64(mix);
}

static std::string shell_quote(const std::string &s) {
    std::string out = "'";
    for (char ch : s) {
        if (ch == '\'') {
            out += "'\\''";
        } else {
            out += ch;
        }
    }
    out += "'";
    return out;
}

static std::string strip_dir(const std::string &path) {
    size_t pos = path.find_last_of("/\\");
    if (pos == std::string::npos) return path;
    return path.substr(pos + 1);
}

static std::string basename_no_ext(const std::string &path) {
    std::string base = strip_dir(path);
    size_t dot = base.find_last_of('.');
    if (dot == std::string::npos) return base;
    return base.substr(0, dot);
}

static std::string resolve_existing_path(const std::string &path) {
    if (path.empty()) return "";
    std::error_code ec;
    if (std::filesystem::exists(path, ec)) return path;
    std::string base = strip_dir(path);
    if (base != path && std::filesystem::exists(base, ec)) return base;
    return "";
}

static std::string format_p_tag(double p) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << std::setprecision(6) << p;
    std::string s = oss.str();
    size_t dot = s.find('.');
    if (dot != std::string::npos) {
        while (!s.empty() && s.back() == '0') s.pop_back();
        if (!s.empty() && s.back() == '.') s.pop_back();
    }
    if (s.empty()) s = "0";
    return s;
}

static std::string make_run_id() {
    auto now = std::chrono::system_clock::now();
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count();
    long long sec = static_cast<long long>(us / 1000000);
    long long usec = static_cast<long long>(us % 1000000);
    std::ostringstream oss;
    oss << sec << "_" << std::setw(6) << std::setfill('0') << usec;
    return oss.str();
}

static bool ensure_dir(const std::string &path) {
    if (path.empty()) return false;
    std::error_code ec;
    if (std::filesystem::exists(path, ec)) {
        return std::filesystem::is_directory(path, ec);
    }
    return std::filesystem::create_directories(path, ec);
}

static bool compress_xz(const std::string &path) {
    std::ostringstream cmd;
    cmd << "xz -z -f " << shell_quote(path);
    int rc = std::system(cmd.str().c_str());
    return rc == 0;
}

static bool parse_int_after(const std::string &line, const std::string &key, int &out) {
    size_t pos = line.find(key);
    if (pos == std::string::npos) return false;
    pos += key.size();
    if (pos >= line.size()) return false;
    size_t end = pos;
    while (end < line.size() && std::isdigit(static_cast<unsigned char>(line[end]))) {
        ++end;
    }
    if (end == pos) return false;
    out = std::stoi(line.substr(pos, end - pos));
    return true;
}

static std::vector<int> diff_indices_bits(
    const std::vector<int> &truth,
    const std::vector<int> &est,
    bool use_xbit
) {
    std::vector<int> idx;
    size_t n = std::min(truth.size(), est.size());
    idx.reserve(n / 8);
    for (size_t i = 0; i < n; ++i) {
        int t = use_xbit ? xbit(truth[i]) : zbit(truth[i]);
        int e = use_xbit ? xbit(est[i]) : zbit(est[i]);
        if ((t ^ e) != 0) idx.push_back(static_cast<int>(i));
    }
    return idx;
}

static int unsatisfied_count(
    const std::vector<int> &hat,
    const std::vector<int> &target
) {
    int count = 0;
    size_t n = std::min(hat.size(), target.size());
    for (size_t i = 0; i < n; ++i) {
        if ((hat[i] ^ target[i]) != 0) count++;
    }
    return count;
}

static std::vector<int> unsatisfied_indices(
    const std::vector<int> &hat,
    const std::vector<int> &target
) {
    std::vector<int> idx;
    size_t n = std::min(hat.size(), target.size());
    idx.reserve(n / 8);
    for (size_t i = 0; i < n; ++i) {
        if ((hat[i] ^ target[i]) != 0) idx.push_back(static_cast<int>(i));
    }
    return idx;
}

static void print_index_list(const std::string &label, const std::vector<int> &idxs) {
    std::cout << label << "_count=" << idxs.size();
    std::cout << " " << label << "=[";
    for (size_t i = 0; i < idxs.size(); ++i) {
        if (i) std::cout << ",";
        std::cout << idxs[i];
    }
    std::cout << "]\n";
}

static void write_index_list(std::ostream &os, const std::string &label, const std::vector<int> &idxs) {
    os << label << "_count=" << idxs.size();
    os << " " << label << "=[";
    for (size_t i = 0; i < idxs.size(); ++i) {
        if (i) os << ",";
        os << idxs[i];
    }
    os << "]";
}

static std::vector<int> sorted_unique(std::vector<int> v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

static std::vector<std::string> sorted_unique_strings(std::vector<std::string> v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

static std::string format_progress_line(
    long long trials_done,
    long long failures,
    long long pp_success,
    long long ets_pp_success,
    long long flip_pp_success,
    long long stab_success,
    long long total_iters,
    double elapsed_sec
) {
    if (trials_done <= 0) return "";
    double fer = static_cast<double>(failures) / static_cast<double>(trials_done);
    double fps = (elapsed_sec > 0.0) ? (static_cast<double>(trials_done) / elapsed_sec) : 0.0;
    double avg_iter = static_cast<double>(total_iters) / static_cast<double>(trials_done);
    std::ostringstream oss;
    oss << "PROGRESS"
        << " | trials=" << trials_done
        << " | failures=" << failures
        << " | pp_success=" << pp_success
        << " | ets_success=" << ets_pp_success
        << " | flip_success=" << flip_pp_success
        << " | stab_success=" << stab_success
        << " | avg_iter=" << std::setprecision(2) << std::fixed << avg_iter
        << " | iters=" << total_iters
        << " | FER=" << std::setprecision(6) << std::fixed << fer
        << " | trials/sec=" << std::setprecision(2) << std::fixed << fps;
    return oss.str();
}

static std::string format_progress_line_compact(
    long long trials_done,
    long long failures,
    long long total_iters,
    double elapsed_sec
) {
    if (trials_done <= 0) return "";
    double fer = static_cast<double>(failures) / static_cast<double>(trials_done);
    double fps = (elapsed_sec > 0.0) ? (static_cast<double>(trials_done) / elapsed_sec) : 0.0;
    double avg_iter = static_cast<double>(total_iters) / static_cast<double>(trials_done);
    std::ostringstream oss;
    oss << "PROGRESS | t=" << trials_done
        << " | fail=" << failures
        << " | fer=" << std::setprecision(6) << std::fixed << fer
        << " | avg=" << std::setprecision(2) << std::fixed << avg_iter
        << " | fps=" << std::setprecision(2) << std::fixed << fps;
    return oss.str();
}

static void report_progress(
    long long trials_done,
    long long failures,
    long long bp_failures,
    long long pp_success,
    long long ets_pp_success,
    long long flip_pp_success,
    long long osd_pp_success,
    long long stab_success,
    long long ets_used,
    const std::vector<std::string> &ets_labels,
    const std::vector<long long> &ets_used_counts,
    long long total_iters,
    double elapsed_sec,
    long long k_value,
    double avg_latency_sec
) {
    if (trials_done <= 0) return;
    double fer = static_cast<double>(failures) / static_cast<double>(trials_done);
    double fps = (elapsed_sec > 0.0) ? (static_cast<double>(trials_done) / elapsed_sec) : 0.0;
    double qbps = fps * static_cast<double>(k_value);
    double ets_rate = (trials_done > 0) ? (static_cast<double>(ets_used) / trials_done) : 0.0;
    double avg_iter = static_cast<double>(total_iters) / static_cast<double>(trials_done);
    double pp_rate = (trials_done > 0) ? (static_cast<double>(pp_success) / trials_done) : 0.0;

    std::vector<std::string> lines;
    const int cols = 3;
    const int key_width = 12;
    const int value_width = 9;
    const int gap = 1;
    auto append_lines = [&](const std::vector<std::pair<std::string, std::string>> &items) {
        auto chunk = build_kv_columns(items, cols, key_width, value_width, gap);
        lines.insert(lines.end(), chunk.begin(), chunk.end());
    };

    append_lines({
        {"trials", std::to_string(trials_done)},
        {"failures", std::to_string(failures)},
        {"bp_fail", std::to_string(bp_failures)}
    });
    append_lines({
        {"FER", format_double_fixed(fer, 6)},
        {"elapsed_s", format_double_fixed(elapsed_sec, 2) + "s"},
        {"latency", format_latency(avg_latency_sec)}
    });
    append_lines({
        {"fps", format_double_fixed(fps, 2)},
        {"qbps", format_qbps(qbps)},
        {"avg_iter", format_double_fixed(avg_iter, 2)}
    });
    append_lines({
        {"iters", std::to_string(total_iters)},
        {"stab_success", std::to_string(stab_success)},
        {"pp_success", std::to_string(pp_success)},
    });
    append_lines({
        {"pp_rate", format_double_fixed(pp_rate, 4)},
        {"pp_ets", std::to_string(ets_pp_success)},
        {"pp_flip", std::to_string(flip_pp_success)}
    });
    append_lines({
        {"pp_osd", std::to_string(osd_pp_success)},
        {"ets_used", std::to_string(ets_used)},
        {"ets_rate", format_double_fixed(ets_rate, 4)}
    });

    std::vector<std::pair<std::string, std::string>> ets_items;
    ets_items.reserve(ets_labels.size());
    size_t count = std::min(ets_labels.size(), ets_used_counts.size());
    for (size_t i = 0; i < count; ++i) {
        ets_items.push_back({"ets" + ets_labels[i], std::to_string(ets_used_counts[i])});
    }
    append_lines(ets_items);

    print_progress_block(lines);
}

static void print_usage(const char *prog) {
    print_ascii_banner(std::cerr);
    print_tips(std::cerr, {
        "Provide --params FILE to load a code definition.",
        "Use --simulate for random trials or --sx/--sz/--err for single decode.",
        "Add --help to see all available options."
    });

    print_help_section("Usage");
    std::cerr << "  " << prog << " --params FILE [options]\n";

    print_help_section("Required");
    std::cerr << "  --params FILE\n";
    std::cerr << "    Code-parameter file. Expected format: ../apm_g8/...\n";

    print_help_section("Input Modes");
    std::cerr << "  --simulate\n";
    std::cerr << "    Generate random errors and run decoding.\n";
    std::cerr << "    Use with: --p, --trials, --min-fails, --seed, --trial-index\n";
    std::cerr << "  --sx FILE --sz FILE\n";
    std::cerr << "    Read syndrome from files and decode once.\n";
    std::cerr << "  --err FILE\n";
    std::cerr << "    Read the true error vector and decode once.\n";
    std::cerr << "    Mutually exclusive with --sx/--sz and --simulate.\n";

    print_help_section("Error Probability / Iterations / Damping");
    std::cerr << "  --p N           Error rate (default: 0.05)\n";
    std::cerr << "  --max-iter N    Maximum iterations (default: 50)\n";
    std::cerr << "  --damping R     Damping factor 0-1 (default: 0.0)\n";

    print_help_section("Randomness and Reproducibility");
    std::cerr << "  --seed N        Base RNG seed (default: current time)\n";
    std::cerr << "  --trial-index N Reproduce a specific trial (requires: --simulate --trials 1)\n";

    print_help_section("Trials and Early Stop");
    std::cerr << "  --trials N      Number of trials (simulate only, default: 1)\n";
    std::cerr << "  --min-fails N   Stop when failures reach N (simulate only)\n";
    std::cerr << "  --report-every N Progress print interval (in trials, default: 100)\n";

    print_help_section("Progress TSV");
    std::cerr << "  --progress-tsv FILE  Progress TSV output path (overwrites with latest line)\n";
    std::cerr << "                    Default: <prefix>/progress_p<...>_<run>.tsv\n";
    std::cerr << "  --progress-every N   TSV update interval (in trials, default: 100)\n";

    print_help_section("Output and Details");
    std::cerr << "  --verbose 0|1   Show per-iteration details (default: 0)\n";
    std::cerr << "    Show us/diff only after stagnation is detected.\n";
    std::cerr << "  --verbose-all   Show us/diff at every iteration (forces verbose)\n";
    std::cerr << "  --flip-hist N   History window after stagnation for flip-PP (default: 5)\n";
    std::cerr << "    If 0, do not use history union; use only the last flip.\n";
    std::cerr << "  --freeze-syn    Freeze BP on a side once its syndrome is satisfied.\n";
    std::cerr << "  --no-pp         Disable PP (ETS/flip).\n";
    std::cerr << "  --report-fail   Print summary on failure to stdout.\n";
    std::cerr << "  --report-ets    Print detailed ETS application status.\n";
    std::cerr << "  --est FILE      Write estimated error vector (trials=1 only)\n";

    print_help_section("ETS Files");
    std::cerr << "  --ets6-x FILE --ets6-z FILE\n";
    std::cerr << "  --ets12-x FILE --ets12-z FILE\n";
    std::cerr << "  --ets8-path4 FILE\n";
    std::cerr << "    Explicitly set ets_8_2_path4.txt (used for both X/Z)\n";
    std::cerr << "    By default, use files under the directory named by the H file prefix:\n";
    std::cerr << "    ets_6_2.txt / ets_8_2.txt / ets_8_2_path4.txt\n";
    std::cerr << "    ets_10_2.txt / ets_12_2.txt / ets_14_2.txt / ets_16_2.txt\n";
    std::cerr << "    ets_18_2.txt / ets_20_2.txt / ets_22_2.txt.\n";
    std::cerr << "    In the merged format, prefix lines with 'X ' or 'Z ' to distinguish.\n";
    std::cerr << "    If no prefix is given, it is used for both.\n";
    std::cerr << "  --ets-x FILE --ets-z FILE\n";
    std::cerr << "    For compatibility. Same as --ets12-x/--ets12-z.\n";

    print_help_section("Input File Formats");
    std::cerr << "  --sx/--sz: one 0/1 bit per line.\n";
    std::cerr << "  --err: one 0/1/2/3 per line (0=I,1=X,2=Z,3=Y).\n";

    print_help_section("Output Notes");
    std::cerr << "  syndrome_match=true means the estimated error's syndrome matches the input.\n";
    std::cerr << "  It does not necessarily match the true error.\n";
}

}  // namespace

int main(int argc, char **argv) {
    if (argc == 1) {
        print_usage(argv[0]);
        return 0;
    }
    std::string params_path;
    std::string sx_path;
    std::string sz_path;
    std::string est_path;
    std::string err_path;
    std::string ets6_x_path = "ets_6_2.txt";
    std::string ets6_z_path = "ets_6_2.txt";
    std::string ets12_x_path = "ets_12_2.txt";
    std::string ets12_z_path = "ets_12_2.txt";
    std::string ets8_path4_path = "ets_8_2_path4.txt";
    bool ets6_x_set = false;
    bool ets6_z_set = false;
    bool ets12_x_set = false;
    bool ets12_z_set = false;
    bool ets8_path4_set = false;
    std::string save_fail_prefix;
    double p_err = 0.05;
    int max_iter = 50;
    int flip_hist_window = 5;
    double damping = 0.0;
    bool freeze_syn = false;
    bool simulate = false;
    bool verbose = false;
    bool verbose_all = false;
    unsigned long long seed = 0;
    bool seed_set = false;
    long long trials = 1;
    long long min_fails = 0;
    long long report_every = 100;
    long long progress_every = 100;
    bool report_fail = false;
    bool report_ets = false;
    long long trial_index = -1;
    bool enable_pp = true;
    const bool enable_log_files = false;
    std::string progress_tsv_path;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto need = [&](int n) {
            if (i + n >= argc) {
                print_usage(argv[0]);
                std::exit(1);
            }
        };
        if (arg == "--params") {
            need(1);
            params_path = argv[++i];
        } else if (arg == "--sx") {
            need(1);
            sx_path = argv[++i];
        } else if (arg == "--sz") {
            need(1);
            sz_path = argv[++i];
        } else if (arg == "--err") {
            need(1);
            err_path = argv[++i];
        } else if (arg == "--est") {
            need(1);
            est_path = argv[++i];
        } else if (arg == "--p") {
            need(1);
            p_err = std::stod(argv[++i]);
        } else if (arg == "--max-iter") {
            need(1);
            max_iter = std::stoi(argv[++i]);
        } else if (arg == "--flip-hist") {
            need(1);
            flip_hist_window = std::stoi(argv[++i]);
            if (flip_hist_window < 0) flip_hist_window = 0;
        } else if (arg == "--freeze-syn") {
            freeze_syn = true;
        } else if (arg == "--damping") {
            need(1);
            damping = std::stod(argv[++i]);
        } else if (arg == "--simulate") {
            simulate = true;
        } else if (arg == "--verbose-all") {
            verbose = true;
            verbose_all = true;
        } else if (arg == "--verbose") {
            need(1);
            verbose = (std::stoi(argv[++i]) != 0);
        } else if (arg == "--seed") {
            need(1);
            seed = std::stoull(argv[++i]);
            seed_set = true;
        } else if (arg == "--ets6-x") {
            need(1);
            ets6_x_path = argv[++i];
            ets6_x_set = true;
        } else if (arg == "--ets6-z") {
            need(1);
            ets6_z_path = argv[++i];
            ets6_z_set = true;
        } else if (arg == "--ets12-x" || arg == "--ets-x") {
            need(1);
            ets12_x_path = argv[++i];
            ets12_x_set = true;
        } else if (arg == "--ets12-z" || arg == "--ets-z") {
            need(1);
            ets12_z_path = argv[++i];
            ets12_z_set = true;
        } else if (arg == "--ets8-path4") {
            need(1);
            ets8_path4_path = argv[++i];
            ets8_path4_set = true;
        } else if (arg == "--trials") {
            need(1);
            trials = std::stoll(argv[++i]);
        } else if (arg == "--min-fails") {
            need(1);
            min_fails = std::stoll(argv[++i]);
        } else if (arg == "--report-every") {
            need(1);
            report_every = std::stoll(argv[++i]);
        } else if (arg == "--progress-every") {
            need(1);
            progress_every = std::stoll(argv[++i]);
        } else if (arg == "--report-fail") {
            report_fail = true;
        } else if (arg == "--report-ets") {
            report_ets = true;
        } else if (arg == "--save-fail-prefix") {
            need(1);
            save_fail_prefix = argv[++i];
        } else if (arg == "--progress-tsv") {
            need(1);
            progress_tsv_path = argv[++i];
        } else if (arg == "--no-pp") {
            enable_pp = false;
        } else if (arg == "--trial-index") {
            need(1);
            trial_index = std::stoll(argv[++i]);
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else {
            print_usage(argv[0]);
            return 1;
        }
    }

    (void)save_fail_prefix;

    if (!seed_set) {
        auto now = std::chrono::system_clock::now();
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(
            now.time_since_epoch()).count();
        seed = static_cast<unsigned long long>(us);
    }

    if (params_path.empty()) {
        std::cerr << "--params is required\n";
        return 1;
    }

    Params params;
    if (!read_params(params_path, params)) {
        std::cerr << "Failed to read params from " << params_path << "\n";
        return 1;
    }
    g_inline_progress_enabled = is_tty_stdout() && !verbose && !verbose_all && !report_fail && !report_ets;
    print_stdout_header("JointBP ETS");
    print_tips(std::cout, {
        "Use --help to see all options.",
        "Use --simulate for randomized trials.",
        "Use --report-every to print progress periodically."
    });
    print_stdout_section("Input");
    std::cout << "params_path=" << params_path << "\n";
    std::cout << "fi=" << format_affine_pairs(params.a, params.b) << "\n";
    std::cout << "gi=" << format_affine_pairs(params.c, params.d) << "\n";
    std::string log_dir = basename_no_ext(params_path);
    std::string p_tag = format_p_tag(p_err);
    if (simulate && enable_log_files) {
        if (!ensure_dir(log_dir)) {
            std::cerr << "Failed to create log directory: " << log_dir << "\n";
            return 1;
        }
    }
    if (!ets6_x_set) {
        ets6_x_path = log_dir + "/ets_6_2.txt";
    }
    if (!ets6_z_set) {
        ets6_z_path = log_dir + "/ets_6_2.txt";
    }
    if (!ets12_x_set) {
        ets12_x_path = log_dir + "/ets_12_2.txt";
    }
    if (!ets12_z_set) {
        ets12_z_path = log_dir + "/ets_12_2.txt";
    }
    if (!ets8_path4_set) {
        ets8_path4_path = log_dir + "/ets_8_2_path4.txt";
    }

    bool progress_enabled = (enable_log_files && trials > 1);
    if (progress_enabled) {
        if (progress_every <= 0) {
            std::cerr << "--progress-every must be 1 or greater\n";
            return 1;
        }
        if (progress_tsv_path.empty()) {
            std::string run_id = make_run_id();
            progress_tsv_path = log_dir + "/progress_p" + p_tag + "_" + run_id + ".tsv";
        }
        {
            std::filesystem::path p(progress_tsv_path);
            std::filesystem::path parent = p.parent_path();
            if (!parent.empty()) {
                if (!ensure_dir(parent.string())) {
                    std::cerr << "Failed to create progress TSV directory: " << parent.string() << "\n";
                    return 1;
                }
            }
        }
        {
            std::ofstream progress_tsv(progress_tsv_path, std::ios::trunc);
            if (!progress_tsv) {
                std::cerr << "Failed to write progress TSV: " << progress_tsv_path << "\n";
                return 1;
            }
        }
    }

    bool progress_tsv_failed = false;
    auto write_progress_tsv = [&](long long trials_done,
                                  long long failures,
                                  long long pp_success,
                                  long long ets_pp_success,
                                  long long flip_pp_success,
                                  long long stab_success,
                                  long long total_iters,
                                  double elapsed_sec) {
        if (!progress_enabled || trials_done <= 0 || progress_tsv_failed) return;
        std::string line = format_progress_line(trials_done, failures, pp_success, ets_pp_success,
                                                flip_pp_success, stab_success, total_iters, elapsed_sec);
        if (line.empty()) return;
        std::ofstream progress_tsv(progress_tsv_path, std::ios::trunc);
        if (!progress_tsv) {
            std::cerr << "Failed to write progress TSV: " << progress_tsv_path << "\n";
            progress_tsv_failed = true;
            return;
        }
        double fer = static_cast<double>(failures) / static_cast<double>(trials_done);
        double fps = (elapsed_sec > 0.0) ? (static_cast<double>(trials_done) / elapsed_sec) : 0.0;
        double avg_iter = static_cast<double>(total_iters) / static_cast<double>(trials_done);
        auto ts_epoch = std::chrono::duration_cast<std::chrono::seconds>(
                            std::chrono::system_clock::now().time_since_epoch())
                            .count();
        progress_tsv << ts_epoch
                     << "\t" << trials_done
                     << "\t" << failures
                     << "\t" << pp_success
                     << "\t" << ets_pp_success
                     << "\t" << flip_pp_success
                     << "\t" << stab_success
                     << "\t" << std::setprecision(2) << std::fixed << avg_iter
                     << "\t" << total_iters
                     << "\t" << std::setprecision(6) << std::fixed << fer
                     << "\t" << std::setprecision(2) << std::fixed << fps
                     << "\t" << line
                     << "\n";
    };

    auto hx_blocks = build_block_matrix_hx(params.P, params.J, params.L2, params.a, params.b, params.c, params.d);
    auto hz_blocks = build_block_matrix_hz(params.P, params.J, params.L2, params.a, params.b, params.c, params.d);
    auto x_checks = build_check_neighbors(hx_blocks, params.P);
    auto z_checks = build_check_neighbors(hz_blocks, params.P);
    int nvars = params.L * static_cast<int>(params.P);

    auto var_to_x = build_var_adjacency(x_checks, nvars);
    auto var_to_z = build_var_adjacency(z_checks, nvars);

    bool need_basis = simulate || !err_path.empty();
    RowBasis hx_basis;
    RowBasis hz_basis;
    RowBasis rank_basis_x_storage;
    RowBasis rank_basis_z_storage;
    const RowBasis *rank_basis_x = nullptr;
    const RowBasis *rank_basis_z = nullptr;
    if (need_basis) {
        hx_basis = build_row_basis(x_checks, nvars);
        hz_basis = build_row_basis(z_checks, nvars);
        rank_basis_x = &hx_basis;
        rank_basis_z = &hz_basis;
    } else {
        rank_basis_x_storage = build_row_basis(x_checks, nvars);
        rank_basis_z_storage = build_row_basis(z_checks, nvars);
        rank_basis_x = &rank_basis_x_storage;
        rank_basis_z = &rank_basis_z_storage;
    }

    DegreeStats x_var_stats = compute_var_degree_stats(var_to_x);
    DegreeStats z_var_stats = compute_var_degree_stats(var_to_z);
    DegreeStats x_check_stats = compute_check_degree_stats(x_checks);
    DegreeStats z_check_stats = compute_check_degree_stats(z_checks);
    OrthCheckResult orth = check_hx_hz_orthogonality(x_checks, z_checks, nvars);
    int girth_x = compute_girth(x_checks, var_to_x);
    int girth_z = compute_girth(z_checks, var_to_z);
    long long mx = static_cast<long long>(x_checks.size());
    long long mz = static_cast<long long>(z_checks.size());
    long long rank_x = static_cast<long long>(rank_basis_x->rows.size());
    long long rank_z = static_cast<long long>(rank_basis_z->rows.size());
    long long stab_rank = rank_x + rank_z;
    long long k = static_cast<long long>(nvars) - stab_rank;
    double rate = (nvars > 0) ? (static_cast<double>(k) / static_cast<double>(nvars)) : 0.0;

    print_stdout_section("Code Summary");
    std::vector<std::pair<std::string, std::string>> summary_items = {
        {"P", std::to_string(params.P)},
        {"J", std::to_string(params.J)},
        {"L", std::to_string(params.L)},
        {"L2", std::to_string(params.L2)},
        {"p", format_double_fixed(p_err, 6)},
        {"n", std::to_string(nvars)},
        {"mx", std::to_string(mx)},
        {"mz", std::to_string(mz)},
        {"rank_x", std::to_string(rank_x)},
        {"rank_z", std::to_string(rank_z)},
        {"stab_rank", std::to_string(stab_rank)},
        {"k", std::to_string(k)},
        {"rate", format_double_fixed(rate, 6)}
    };
    print_kv_compact(std::cout, summary_items, kBoxWidth);
    std::vector<std::pair<std::string, std::string>> x_items = {
        {"edges", std::to_string(x_var_stats.edges)},
        {"var_deg", format_deg_triplet(x_var_stats.min, x_var_stats.avg, x_var_stats.max)},
        {"check_deg", format_deg_triplet(x_check_stats.min, x_check_stats.avg, x_check_stats.max)}
    };
    std::vector<std::pair<std::string, std::string>> z_items = {
        {"edges", std::to_string(z_var_stats.edges)},
        {"var_deg", format_deg_triplet(z_var_stats.min, z_var_stats.avg, z_var_stats.max)},
        {"check_deg", format_deg_triplet(z_check_stats.min, z_check_stats.avg, z_check_stats.max)}
    };
    print_kv_compact_prefixed(std::cout, "X: ", x_items, kBoxWidth);
    print_kv_compact_prefixed(std::cout, "Z: ", z_items, kBoxWidth);

    print_stdout_section("Graph Checks");
    std::vector<std::pair<std::string, std::string>> graph_items = {
        {"hx_hz_orthogonal", tf(orth.orthogonal)},
        {"hx_hz_odd_pairs", std::to_string(orth.odd_pairs)},
        {"girth_x", format_girth(girth_x)},
        {"girth_z", format_girth(girth_z)}
    };
    if (!orth.orthogonal && orth.first_x >= 0 && orth.first_z >= 0) {
        graph_items.push_back({"hx_hz_first_odd_pair",
                               std::to_string(orth.first_x) + "," +
                                   std::to_string(orth.first_z)});
        graph_items.push_back({"hx_hz_odd_overlap", std::to_string(orth.first_overlap)});
    }
    print_kv_compact(std::cout, graph_items, kBoxWidth);

    struct ETSGroup {
        int size = 0;
        std::string label;
        ETSIndex x;
        ETSIndex z;
        bool x_loaded = false;
        bool z_loaded = false;
        int x_count = 0;
        int z_count = 0;
    };
    enum class ETSLoadState { Missing, Empty, Loaded };
    struct ETSLoadInfo {
        std::string label;
        int x_count = 0;
        int z_count = 0;
        ETSLoadState x_state = ETSLoadState::Missing;
        ETSLoadState z_state = ETSLoadState::Missing;
    };
    const std::vector<int> ets_sizes = {6, 8, 10, 12, 14, 16, 18, 20, 22};
    const std::vector<std::string> ets_labels = {"6", "8", "8p4", "10", "12", "14", "16", "18", "20", "22"};
    std::unordered_map<std::string, size_t> ets_label_index;
    for (size_t i = 0; i < ets_labels.size(); ++i) {
        ets_label_index[ets_labels[i]] = i;
    }
    std::vector<ETSGroup> ets_groups;
    ets_groups.reserve(ets_sizes.size() + 1);
    auto load_ets = [&](int size,
                        const std::string &path,
                        char side,
                        ETSIndex &out,
                        bool &loaded,
                        int &count_out,
                        ETSLoadState &state_out,
                        bool warn_missing) {
        bool warn = warn_missing && (report_ets || verbose);
        count_out = 0;
        state_out = ETSLoadState::Missing;
        std::string resolved = resolve_existing_path(path);
        if (resolved.empty()) {
            if (warn) {
                std::cerr << "Warning: failed to read " << path << " (ETS" << size
                          << " " << side << " disabled)\n";
            }
            loaded = false;
            return;
        }
        if (report_ets && resolved != path) {
            std::cout << "[ETS] fallback path=" << path << " -> " << resolved << "\n";
        }
        std::error_code ec;
        auto file_size = std::filesystem::file_size(resolved, ec);
        if (!ec && file_size == 0) {
            state_out = ETSLoadState::Empty;
            loaded = false;
            return;
        }
        if (!read_ets_file(resolved, static_cast<int>(params.P), size, side, out)) {
            if (warn) {
                std::cerr << "Warning: failed to read " << resolved << " (ETS" << size
                          << " " << side << " disabled)\n";
            }
            loaded = false;
            return;
        }
        if (out.entries.empty()) {
            if (warn) {
                std::cerr << "Warning: no ETS" << size << " " << side
                          << " entries loaded from " << resolved << "\n";
            }
            state_out = ETSLoadState::Empty;
            loaded = false;
            return;
        }
        loaded = true;
        state_out = ETSLoadState::Loaded;
        count_out = static_cast<int>(out.entries.size());
        if (report_ets) {
            std::cout << "[ETS] loaded size=" << size << " side=" << side
                      << " entries=" << out.entries.size()
                      << " path=" << resolved << "\n";
        }
    };
    std::vector<ETSLoadInfo> ets_load_infos;
    auto add_group = [&](int size,
                         const std::string &label,
                         const std::string &path_x,
                         const std::string &path_z,
                         bool warn_missing) {
        ETSGroup group;
        group.size = size;
        group.label = label;
        ETSLoadState x_state = ETSLoadState::Missing;
        ETSLoadState z_state = ETSLoadState::Missing;
        load_ets(size, path_x, 'X', group.x, group.x_loaded, group.x_count, x_state, warn_missing);
        load_ets(size, path_z, 'Z', group.z, group.z_loaded, group.z_count, z_state, warn_missing);
        ETSLoadInfo info;
        info.label = label;
        info.x_count = group.x_count;
        info.z_count = group.z_count;
        info.x_state = x_state;
        info.z_state = z_state;
        ets_load_infos.push_back(info);
        if (group.x_loaded || group.z_loaded) {
            ets_groups.push_back(std::move(group));
        }
    };

    for (int size : ets_sizes) {
        std::string path_x;
        std::string path_z;
        if (size == 6) {
            path_x = ets6_x_path;
            path_z = ets6_z_path;
        } else if (size == 12) {
            path_x = ets12_x_path;
            path_z = ets12_z_path;
        } else {
            path_x = log_dir + "/ets_" + std::to_string(size) + "_2.txt";
            path_z = path_x;
        }
        bool warn_missing = (size == 6 || size == 12);
        add_group(size, std::to_string(size), path_x, path_z, warn_missing);
        if (size == 8) {
            add_group(size, "8p4", ets8_path4_path, ets8_path4_path, ets8_path4_set);
        }
    }
    print_stdout_section("ETS Load");
    if (ets_load_infos.empty()) {
        std::cout << "ETS=none\n";
    } else {
        auto is_loaded = [&](ETSLoadState state) {
            return state == ETSLoadState::Loaded;
        };
        auto is_empty = [&](ETSLoadState state) {
            return state == ETSLoadState::Empty;
        };
        auto row_status = [&](const ETSLoadInfo &info) {
            bool any_loaded = is_loaded(info.x_state) || is_loaded(info.z_state);
            bool any_empty = is_empty(info.x_state) || is_empty(info.z_state);
            if (any_loaded) return std::string("ok");
            if (any_empty) return std::string("empty");
            return std::string("missing");
        };
        size_t label_width = 0;
        size_t x_val_width = 1;
        size_t z_val_width = 1;
        size_t status_width = 0;
        for (const auto &info : ets_load_infos) {
            std::string label = "ETS" + info.label + ":";
            label_width = std::max(label_width, label.size());
            std::string x_val = (info.x_state == ETSLoadState::Missing)
                                    ? "-"
                                    : std::to_string(info.x_count);
            std::string z_val = (info.z_state == ETSLoadState::Missing)
                                    ? "-"
                                    : std::to_string(info.z_count);
            x_val_width = std::max(x_val_width, x_val.size());
            z_val_width = std::max(z_val_width, z_val.size());
            status_width = std::max(status_width, row_status(info).size());
        }
        for (const auto &info : ets_load_infos) {
            std::string label = "ETS" + info.label + ":";
            std::string x_val = (info.x_state == ETSLoadState::Missing)
                                    ? "-"
                                    : std::to_string(info.x_count);
            std::string z_val = (info.z_state == ETSLoadState::Missing)
                                    ? "-"
                                    : std::to_string(info.z_count);
            std::string status = row_status(info);
            std::ostringstream line;
            line << std::left << std::setw(static_cast<int>(label_width)) << label << " ";
            line << "X=" << std::right << std::setw(static_cast<int>(x_val_width)) << x_val << " ";
            line << "Z=" << std::right << std::setw(static_cast<int>(z_val_width)) << z_val << " ";
            line << "status=" << std::left << std::setw(static_cast<int>(status_width)) << status;
            std::cout << line.str() << "\n";
        }
    }
    std::vector<ETSIndexRef> ets_refs;
    ets_refs.reserve(ets_groups.size());
    for (const auto &group : ets_groups) {
        ETSIndexRef ref;
        ref.size = group.size;
        ref.label = group.label;
        ref.x = group.x_loaded ? &group.x : nullptr;
        ref.z = group.z_loaded ? &group.z : nullptr;
        if (ref.x || ref.z) {
            ets_refs.push_back(ref);
        }
    }

    CycleIndex cycles_index;
    if (enable_pp && (simulate || !err_path.empty())) {
        load_cycles_index(log_dir, static_cast<int>(params.P), cycles_index);
    }
    const CycleIndex *cycles_ptr = cycles_index.any_loaded ? &cycles_index : nullptr;

    std::vector<int> sx;
    std::vector<int> sz;
    std::vector<int> err;

    if (simulate && !err_path.empty()) {
        std::cerr << "--err cannot be used with --simulate\n";
        return 1;
    }
    if (simulate && trial_index >= 0) {
        if (trials != 1 || min_fails > 0) {
            std::cerr << "--trial-index requires --trials 1 and no --min-fails\n";
            return 1;
        }
    }

    if (!simulate) {
        if (trials != 1 || min_fails > 0) {
            std::cerr << "--trials/--min-fails are only valid with --simulate\n";
            return 1;
        }
        if (!err_path.empty() && (!sx_path.empty() || !sz_path.empty())) {
            std::cerr << "--err cannot be combined with --sx/--sz\n";
            return 1;
        }
        if (!err_path.empty()) {
            if (!read_int_vector(err_path, err)) {
                std::cerr << "Failed to read error vector from " << err_path << "\n";
                return 1;
            }
            if (static_cast<int>(err.size()) != nvars) {
                std::cerr << "Error vector size mismatch\n";
                return 1;
            }
            for (int v : err) {
                if (v < 0 || v > 3) {
                    std::cerr << "Error vector has invalid value\n";
                    return 1;
                }
            }
            sx = compute_syndrome(x_checks, err, true);
            sz = compute_syndrome(z_checks, err, false);
        } else {
            if (sx_path.empty() || sz_path.empty()) {
                std::cerr << "--sx and --sz are required unless --simulate is set\n";
                return 1;
            }
            if (!read_bit_vector(sx_path, sx) || !read_bit_vector(sz_path, sz)) {
                std::cerr << "Failed to read syndrome files\n";
                return 1;
            }
            if (static_cast<int>(sx.size()) != static_cast<int>(x_checks.size()) ||
                static_cast<int>(sz.size()) != static_cast<int>(z_checks.size())) {
                std::cerr << "Syndrome size mismatch\n";
                return 1;
            }
        }
    }

    Msg prior{1.0 - p_err, p_err / 3.0, p_err / 3.0, p_err / 3.0};
    normalize_msg(prior);

    if (!simulate) {
        const std::vector<int> *truth_ptr = err_path.empty() ? nullptr : &err;
        bool pp_verbose = report_ets || verbose;
        PPEarlyContext pp_ctx{enable_pp, pp_verbose, ets_refs, cycles_ptr};
        std::vector<std::string> pp_log_lines;
        auto latency_start = std::chrono::steady_clock::now();
        auto t_start = latency_start;
        auto res = joint_bp_decode(
            x_checks, z_checks, var_to_x, var_to_z,
            sx, sz, prior, max_iter, flip_hist_window, freeze_syn, damping,
            verbose, verbose_all, enable_pp, enable_pp ? &pp_ctx : nullptr, truth_ptr,
            need_basis ? &hx_basis : nullptr, need_basis ? &hz_basis : nullptr, nullptr,
            &pp_log_lines
        );
        double elapsed_sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();

        bool ok = res.syndrome_match;
        bool ets_pp_ok = res.pp_success_ets;
        bool flip_pp_ok = res.pp_success_flip;
        bool osd_pp_ok = false;
        bool used6 = res.pp_used_ets6;
        bool used12 = res.pp_used_ets12;
        bool used_any = res.pp_used_ets_any || used6 || used12;
        auto pp_log_line = [&](const std::string &line) {
            if (pp_verbose) {
                std::cout << line << "\n";
            }
            pp_log_lines.push_back(line);
        };
        if (!ok && enable_pp) {
            pp_log_line("[PP] start");
            auto sx_hat_before = compute_syndrome(x_checks, res.est, true);
            auto sz_hat_before = compute_syndrome(z_checks, res.est, false);
            auto us_x = unsatisfied_indices(sx_hat_before, sx);
            auto us_z = unsatisfied_indices(sz_hat_before, sz);
            {
                std::ostringstream us_oss;
                us_oss << "[PP-ETS] ";
                write_index_list(us_oss, "us_x", us_x);
                pp_log_line(us_oss.str());
            }
            {
                std::ostringstream us_oss;
                us_oss << "[PP-ETS] ";
                write_index_list(us_oss, "us_z", us_z);
                pp_log_line(us_oss.str());
            }
            std::vector<int> diff_x;
            std::vector<int> diff_z;
            if (truth_ptr) {
                diff_x = diff_indices_bits(*truth_ptr, res.est, true);
                diff_z = diff_indices_bits(*truth_ptr, res.est, false);
                if (diff_x.size() <= 50) {
                    std::ostringstream diff_oss;
                    diff_oss << "[PP] ";
                    write_index_list(diff_oss, "diff_x", diff_x);
                    pp_log_line(diff_oss.str());
                }
                if (diff_z.size() <= 50) {
                    std::ostringstream diff_oss;
                    diff_oss << "[PP] ";
                    write_index_list(diff_oss, "diff_z", diff_z);
                    pp_log_line(diff_oss.str());
                }
                maybe_log_diff_trapset_cycle_union('X', diff_x, us_x, x_checks, var_to_x,
                                                   pp_ctx.cycles, pp_log_line);
                maybe_log_diff_trapset_cycle_union('Z', diff_z, us_z, z_checks, var_to_z,
                                                   pp_ctx.cycles, pp_log_line);
            }
            ETSFixResult ets_x = attempt_ets_fix_side_multi(
                "X", pp_ctx.ets, x_checks, var_to_x, sx, true, res.est, pp_verbose
            );
            ETSFixResult ets_z = attempt_ets_fix_side_multi(
                "Z", pp_ctx.ets, z_checks, var_to_z, sz, false, res.est, pp_verbose
            );
            bool fixed_x = ets_x.applied;
            bool fixed_z = ets_z.applied;
            bool ets_used = fixed_x || fixed_z;
            used_any = used_any || ets_used;
            used6 = used6 || (ets_x.size == 6) || (ets_z.size == 6);
            used12 = used12 || (ets_x.size == 12) || (ets_z.size == 12);
            if (truth_ptr) {
                maybe_log_trapset('X', us_x, fixed_x, diff_x, x_checks, var_to_x,
                                  pp_ctx.cycles, pp_log_line);
                maybe_log_trapset('Z', us_z, fixed_z, diff_z, z_checks, var_to_z,
                                  pp_ctx.cycles, pp_log_line);
            }
            if (ets_used) {
                auto sx_hat = compute_syndrome(x_checks, res.est, true);
                auto sz_hat = compute_syndrome(z_checks, res.est, false);
                ok = (sx_hat == sx) && (sz_hat == sz);
                if (ok) {
                    ets_pp_ok = true;
                }
                std::ostringstream oss;
                oss << "[PP-ETS] done pp-ets=" << tf(ets_used)
                    << " x=" << format_ets_label(ets_x)
                    << " z=" << format_ets_label(ets_z)
                    << " syndrome_ok=" << tf(ok);
                pp_log_line(oss.str());
                if (fixed_x && !us_x.empty()) {
                    std::ostringstream us_oss;
                    us_oss << "[PP-ETS] ";
                    write_index_list(us_oss, "us_x", us_x);
                    pp_log_line(us_oss.str());
                }
                if (fixed_z && !us_z.empty()) {
                    std::ostringstream us_oss;
                    us_oss << "[PP-ETS] ";
                    write_index_list(us_oss, "us_z", us_z);
                    pp_log_line(us_oss.str());
                }
            } else {
                pp_log_line("[PP-ETS] done pp-ets=false x=none z=none syndrome_ok=false");
            }
        }
        if (!ok && enable_pp) {
            bool flip_x = false;
            bool flip_z = false;
            pp_log_line("[flip-PP] start");
            auto x_sets = build_flip_sets(res.flip_x_history, res.flip_x_window, res.flip_x);
            for (const auto &set : x_sets) {
                if (!set.vars || set.vars->empty()) continue;
                if (attempt_flip_pp_side("X", *set.vars, x_checks, sx, true, res.est,
                                         pp_verbose, nullptr, set.tag)) {
                    flip_x = true;
                    break;
                }
            }
            auto z_sets = build_flip_sets(res.flip_z_history, res.flip_z_window, res.flip_z);
            for (const auto &set : z_sets) {
                if (!set.vars || set.vars->empty()) continue;
                if (attempt_flip_pp_side("Z", *set.vars, z_checks, sz, false, res.est,
                                         pp_verbose, nullptr, set.tag)) {
                    flip_z = true;
                    break;
                }
            }
            if (flip_x || flip_z) {
                auto sx_hat = compute_syndrome(x_checks, res.est, true);
                auto sz_hat = compute_syndrome(z_checks, res.est, false);
                ok = (sx_hat == sx) && (sz_hat == sz);
                if (ok) {
                    flip_pp_ok = true;
                }
            }
            std::ostringstream oss;
            oss << "[flip-PP] done applied_x=" << tf(flip_x)
                << " applied_z=" << tf(flip_z)
                << " syndrome_ok=" << tf(ok);
            pp_log_line(oss.str());
        }
        if (!ok && enable_pp) {
            const int osd_us_limit = 20;
            bool osd_x = false;
            bool osd_z = false;
            auto sx_hat = compute_syndrome(x_checks, res.est, true);
            auto sz_hat = compute_syndrome(z_checks, res.est, false);
            int osd_us_x = unsatisfied_count(sx_hat, sx);
            int osd_us_z = unsatisfied_count(sz_hat, sz);
            int osd_us = std::max(osd_us_x, osd_us_z);
            if (osd_us > osd_us_limit) {
                std::ostringstream oss;
                oss << "[osd-PP] skipped us=(" << osd_us_x << "," << osd_us_z << ")";
                pp_log_line(oss.str());
            } else {
                pp_log_line("[osd-PP] start");
                osd_x = attempt_osd_pp_side("X", res.abs_llr_x, x_checks, sx, true, res.est,
                                            pp_verbose, truth_ptr);
                osd_z = attempt_osd_pp_side("Z", res.abs_llr_z, z_checks, sz, false, res.est,
                                            pp_verbose, truth_ptr);
                if (osd_x || osd_z) {
                    auto sx_hat2 = compute_syndrome(x_checks, res.est, true);
                    auto sz_hat2 = compute_syndrome(z_checks, res.est, false);
                    ok = (sx_hat2 == sx) && (sz_hat2 == sz);
                    if (ok) {
                        osd_pp_ok = true;
                    }
                }
                std::ostringstream oss;
                oss << "[osd-PP] done applied_x=" << tf(osd_x)
                    << " applied_z=" << tf(osd_z)
                    << " syndrome_ok=" << tf(ok);
                pp_log_line(oss.str());
            }
        }

        auto sx_hat = compute_syndrome(x_checks, res.est, true);
        auto sz_hat = compute_syndrome(z_checks, res.est, false);
        auto us_x = unsatisfied_indices(sx_hat, sx);
        auto us_z = unsatisfied_indices(sz_hat, sz);

        bool exact = false;
        bool stab_equiv = false;
        bool success = ok;
        if (!err_path.empty()) {
            exact = (res.est == err);
            if (ok && need_basis) {
                std::vector<uint64_t> diff_x_bits(hz_basis.words, 0);
                std::vector<uint64_t> diff_z_bits(hx_basis.words, 0);
                fill_diff_bits(diff_x_bits, err, res.est, true);
                fill_diff_bits(diff_z_bits, err, res.est, false);
                stab_equiv = in_row_space(hz_basis, diff_x_bits) && in_row_space(hx_basis, diff_z_bits);
            }
            if (need_basis) {
                success = ok && stab_equiv;
            }
        }
        bool pp_success = success && (res.pp_success || ets_pp_ok || flip_pp_ok || osd_pp_ok);
        bool pp_success_ets = success && ets_pp_ok;
        bool pp_success_flip = success && flip_pp_ok;
        double latency_sec = std::chrono::duration<double>(
                                 std::chrono::steady_clock::now() - latency_start)
                                 .count();
        if (pp_success && enable_log_files) {
            if (!ensure_dir(log_dir)) {
                std::cerr << "Failed to create log directory: " << log_dir << "\n";
            } else {
                std::ostringstream pp_path;
                pp_path << log_dir << "/pp-success_p" << p_tag << "_single.txt";
                std::ofstream pp_out(pp_path.str());
                if (!pp_out) {
                    std::cerr << "Failed to write pp-success file: " << pp_path.str() << "\n";
                } else {
                    bool osd_success = osd_pp_ok ||
                                       (res.pp_success && !res.pp_success_ets && !res.pp_success_flip);
                    std::string pp_method;
                    if (pp_success_ets) pp_method += (pp_method.empty() ? "ets" : ",ets");
                    if (pp_success_flip) pp_method += (pp_method.empty() ? "flip" : ",flip");
                    if (osd_success) pp_method += (pp_method.empty() ? "osd" : ",osd");
                    pp_out << "p=" << std::setprecision(8) << std::fixed << p_err << "\n";
                    pp_out << "bp_ok=" << (res.bp_syndrome_match ? "1" : "0") << "\n";
                    pp_out << "syndrome_match=" << (ok ? "1" : "0") << "\n";
                    pp_out << "decode_success=" << (success ? "1" : "0") << "\n";
                    pp_out << "pp_success=" << (pp_success ? "1" : "0") << "\n";
                    pp_out << "pp_success_in_decode=" << (res.pp_success ? "1" : "0") << "\n";
                    pp_out << "pp_success_ets=" << (pp_success_ets ? "1" : "0") << "\n";
                    pp_out << "pp_success_flip=" << (pp_success_flip ? "1" : "0") << "\n";
                    pp_out << "pp_success_osd=" << (osd_success ? "1" : "0") << "\n";
                    pp_out << "pp_method=" << (pp_method.empty() ? "unknown" : pp_method) << "\n";
                    pp_out << "pp_used_ets_any=" << (used_any ? "1" : "0") << "\n";
                    pp_out << "pp_used_ets6=" << (used6 ? "1" : "0") << "\n";
                    pp_out << "pp_used_ets12=" << (used12 ? "1" : "0") << "\n";
                    pp_out << "pp_log_count=" << pp_log_lines.size() << "\n";
                    pp_out << "pp_log_begin\n";
                    for (const auto &line : pp_log_lines) {
                        pp_out << line << "\n";
                    }
                    pp_out << "pp_log_end\n";
                    pp_out << "repro_cmd=./jointbp_ets --params " << params_path
                           << " --p " << std::setprecision(8) << std::fixed << p_err;
                    if (!err_path.empty()) {
                        pp_out << " --err " << err_path;
                    } else if (!sx_path.empty() && !sz_path.empty()) {
                        pp_out << " --sx " << sx_path << " --sz " << sz_path;
                    }
                    pp_out << " --max-iter " << max_iter
                           << " --flip-hist " << flip_hist_window
                           << " --damping " << std::setprecision(6) << std::fixed << damping
                           << " --verbose 1 --verbose-all";
                    if (freeze_syn) {
                        pp_out << " --freeze-syn";
                    }
                    pp_out << "\n";
                }
            }
        }
        long long stab_success = 0;
        if (stab_equiv && success && !exact) {
            stab_success = 1;
        }

        print_stdout_section("Decode Result");
        std::cout << "iterations=" << res.iterations << " syndrome_match=" << tf(ok) << "\n";
        std::cout << "decode_success=" << tf(success) << "\n";
        std::cout << "pp_enabled=" << tf(enable_pp) << "\n";
        std::cout << "pp_success=" << tf(pp_success) << "\n";
        std::cout << "pp_success_ets=" << tf(pp_success_ets) << "\n";
        std::cout << "pp_success_flip=" << tf(pp_success_flip) << "\n";
        std::cout << "avg_latency=" << format_latency(latency_sec) << "\n";
        write_progress_tsv(1, success ? 0 : 1, pp_success ? 1 : 0, pp_success_ets ? 1 : 0,
                           pp_success_flip ? 1 : 0, stab_success, res.iterations, elapsed_sec);
        bool syn_x = (sx_hat == sx);
        bool syn_z = (sz_hat == sz);
        if (verbose) {
            if (!syn_x) {
                print_index_list("us_x", us_x);
            }
            if (!syn_z) {
                print_index_list("us_z", us_z);
            }
        }
        if (!err_path.empty()) {
            std::cout << "exact_match=" << tf(exact) << "\n";
            std::cout << "stabilizer_equiv=" << tf(stab_equiv) << "\n";
            auto diff_x = diff_indices_bits(err, res.est, true);
            auto diff_z = diff_indices_bits(err, res.est, false);
            if (verbose) {
                if (!syn_x) {
                    print_index_list("diff_x", diff_x);
                }
                if (!syn_z) {
                    print_index_list("diff_z", diff_z);
                }
            }
        }
        if (!est_path.empty()) {
            if (!write_error_vector(est_path, res.est)) {
                std::cerr << "Failed to write estimate to " << est_path << "\n";
                return 1;
            }
            std::cout << "Saved estimate: " << est_path << "\n";
        }
        return success ? 0 : 2;
    }

    if (trials <= 0) {
        std::cerr << "--trials must be positive in --simulate mode\n";
        return 1;
    }
    if (trials != 1 && !est_path.empty()) {
        std::cerr << "--est is only supported with --trials 1\n";
        return 1;
    }

    bool iter_verbose = verbose;
    bool ets_verbose = report_ets || iter_verbose;
    PPEarlyContext pp_ctx{enable_pp, ets_verbose, ets_refs, cycles_ptr};
    const PPEarlyContext *pp_ctx_ptr = enable_pp ? &pp_ctx : nullptr;
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    long long failures = 0;
    long long bp_failures = 0;
    long long pp_success = 0;
    long long flip_pp_success = 0;
    long long osd_pp_success = 0;
    long long ets_used = 0;
    long long ets_saves = 0;
    long long ets6_used = 0;
    long long ets6_saves = 0;
    long long ets12_used = 0;
    long long ets12_saves = 0;
    long long exact_matches = 0;
    long long stab_success = 0;
    long long total_iters = 0;
    long long done = 0;
    long long last_progress_logged = 0;
    double latency_sum_sec = 0.0;
    long long latency_samples = 0;
    std::vector<long long> ets_used_counts(ets_labels.size(), 0);
    auto record_ets_labels = [&](const std::vector<std::string> &labels) {
        for (const auto &label : labels) {
            auto it = ets_label_index.find(label);
            if (it != ets_label_index.end()) {
                ets_used_counts[it->second]++;
            }
        }
    };
    auto add_ets_label = [&](const ETSFixResult &res, std::vector<std::string> &labels) {
        if (!res.applied) return;
        if (!res.label.empty()) {
            labels.push_back(res.label);
        } else if (res.size > 0) {
            labels.push_back(std::to_string(res.size));
        }
    };
    auto start = std::chrono::steady_clock::now();
    std::vector<uint64_t> diff_x_bits;
    std::vector<uint64_t> diff_z_bits;
    if (need_basis) {
        diff_x_bits.assign(hz_basis.words, 0);
        diff_z_bits.assign(hx_basis.words, 0);
    }

    for (long long t = 0; t < trials; ++t) {
        long long trial_idx = (trial_index >= 0) ? trial_index : (t + 1);
        uint64_t trial_seed = trial_seed_from_base(seed, static_cast<uint64_t>(trial_idx));
        std::mt19937_64 rng(trial_seed);
        err.assign(nvars, 0);
        for (int i = 0; i < nvars; ++i) {
            double u = dist(rng);
            if (u < (1.0 - p_err)) {
                err[i] = 0;
            } else {
                double v = (u - (1.0 - p_err)) / p_err;
                if (v < 1.0 / 3.0) err[i] = 1;
                else if (v < 2.0 / 3.0) err[i] = 2;
                else err[i] = 3;
            }
        }
        sx = compute_syndrome(x_checks, err, true);
        sz = compute_syndrome(z_checks, err, false);
        auto latency_start = std::chrono::steady_clock::now();

        std::vector<std::string> pp_log_lines;
        auto res = joint_bp_decode(
            x_checks, z_checks, var_to_x, var_to_z,
            sx, sz, prior, max_iter, flip_hist_window, freeze_syn, damping,
            iter_verbose, verbose_all, enable_pp, pp_ctx_ptr, &err,
            need_basis ? &hx_basis : nullptr, need_basis ? &hz_basis : nullptr, nullptr,
            &pp_log_lines
        );

        bool bp_ok = res.bp_syndrome_match;
        bool used6 = res.pp_used_ets6;
        bool used12 = res.pp_used_ets12;
        bool used_any = res.pp_used_ets_any || used6 || used12;
        std::vector<std::string> res_ets_labels = res.pp_used_ets_labels;
        bool ok = res.syndrome_match;
        bool ets_pp_ok = res.pp_success_ets;
        bool flip_pp_ok = res.pp_success_flip;
        bool osd_pp_ok = false;
        auto pp_log_line = [&](const std::string &line) {
            if (ets_verbose) {
                std::cout << line << "\n";
            }
            pp_log_lines.push_back(line);
        };
        if (!bp_ok) {
            bp_failures++;
        }
        if (res.pp_success && used_any) {
            ets_used++;
            record_ets_labels(res_ets_labels);
            if (used6) ets6_used++;
            if (used12) ets12_used++;
        }
        if (!ok && enable_pp) {
            std::ostringstream oss;
            oss << "[PP] trial=" << trial_idx << " start";
            pp_log_line(oss.str());
            auto sx_hat_before = compute_syndrome(x_checks, res.est, true);
            auto sz_hat_before = compute_syndrome(z_checks, res.est, false);
            auto us_x = unsatisfied_indices(sx_hat_before, sx);
            auto us_z = unsatisfied_indices(sz_hat_before, sz);
            {
                std::ostringstream us_oss;
                us_oss << "[PP-ETS] trial=" << trial_idx << " ";
                write_index_list(us_oss, "us_x", us_x);
                pp_log_line(us_oss.str());
            }
            {
                std::ostringstream us_oss;
                us_oss << "[PP-ETS] trial=" << trial_idx << " ";
                write_index_list(us_oss, "us_z", us_z);
                pp_log_line(us_oss.str());
            }
            std::vector<int> diff_x = diff_indices_bits(err, res.est, true);
            std::vector<int> diff_z = diff_indices_bits(err, res.est, false);
            if (diff_x.size() <= 50) {
                std::ostringstream diff_oss;
                diff_oss << "[PP] trial=" << trial_idx << " ";
                write_index_list(diff_oss, "diff_x", diff_x);
                pp_log_line(diff_oss.str());
            }
            if (diff_z.size() <= 50) {
                std::ostringstream diff_oss;
                diff_oss << "[PP] trial=" << trial_idx << " ";
                write_index_list(diff_oss, "diff_z", diff_z);
                pp_log_line(diff_oss.str());
            }
            maybe_log_diff_trapset_cycle_union('X', diff_x, us_x, x_checks, var_to_x,
                                               pp_ctx.cycles, pp_log_line);
            maybe_log_diff_trapset_cycle_union('Z', diff_z, us_z, z_checks, var_to_z,
                                               pp_ctx.cycles, pp_log_line);
            ETSFixResult ets_x = attempt_ets_fix_side_multi(
                "X", pp_ctx.ets, x_checks, var_to_x, sx, true, res.est, ets_verbose
            );
            ETSFixResult ets_z = attempt_ets_fix_side_multi(
                "Z", pp_ctx.ets, z_checks, var_to_z, sz, false, res.est, ets_verbose
            );
            bool fixed_x = ets_x.applied;
            bool fixed_z = ets_z.applied;
            bool ets_applied = fixed_x || fixed_z;
            used_any = used_any || ets_applied;
            used6 = used6 || (ets_x.size == 6) || (ets_z.size == 6);
            used12 = used12 || (ets_x.size == 12) || (ets_z.size == 12);
            std::vector<std::string> trial_ets_labels;
            add_ets_label(ets_x, trial_ets_labels);
            add_ets_label(ets_z, trial_ets_labels);
            if (!trial_ets_labels.empty()) {
                trial_ets_labels = sorted_unique_strings(std::move(trial_ets_labels));
                record_ets_labels(trial_ets_labels);
            }
            maybe_log_trapset('X', us_x, fixed_x, diff_x, x_checks, var_to_x,
                              pp_ctx.cycles, pp_log_line);
            maybe_log_trapset('Z', us_z, fixed_z, diff_z, z_checks, var_to_z,
                              pp_ctx.cycles, pp_log_line);
            if (ets_applied) {
                ets_used++;
            }
            if (used6) {
                ets6_used++;
            }
            if (used12) {
                ets12_used++;
            }
            if (ets_applied) {
                auto sx_hat = compute_syndrome(x_checks, res.est, true);
                auto sz_hat = compute_syndrome(z_checks, res.est, false);
                ok = (sx_hat == sx) && (sz_hat == sz);
                if (ok) {
                    ets_pp_ok = true;
                }
                std::ostringstream oss2;
                oss2 << "[PP-ETS] trial=" << trial_idx
                     << " done pp-ets=" << tf(ets_applied)
                     << " x=" << format_ets_label(ets_x)
                     << " z=" << format_ets_label(ets_z)
                     << " syndrome_ok=" << tf(ok);
                pp_log_line(oss2.str());
                if (fixed_x && !us_x.empty()) {
                    std::ostringstream us_oss;
                    us_oss << "[PP-ETS] trial=" << trial_idx << " ";
                    write_index_list(us_oss, "us_x", us_x);
                    pp_log_line(us_oss.str());
                }
                if (fixed_z && !us_z.empty()) {
                    std::ostringstream us_oss;
                    us_oss << "[PP-ETS] trial=" << trial_idx << " ";
                    write_index_list(us_oss, "us_z", us_z);
                    pp_log_line(us_oss.str());
                }
            } else {
                std::ostringstream oss2;
                oss2 << "[PP-ETS] trial=" << trial_idx
                     << " done pp-ets=false x=none z=none syndrome_ok=false";
                pp_log_line(oss2.str());
            }
        }
        if (!ok && enable_pp) {
            bool flip_x = false;
            bool flip_z = false;
            std::ostringstream oss;
            oss << "[flip-PP] trial=" << trial_idx << " start";
            pp_log_line(oss.str());
            auto x_sets = build_flip_sets(res.flip_x_history, res.flip_x_window, res.flip_x);
            for (const auto &set : x_sets) {
                if (!set.vars || set.vars->empty()) continue;
                if (attempt_flip_pp_side("X", *set.vars, x_checks, sx, true, res.est, ets_verbose, nullptr, set.tag)) {
                    flip_x = true;
                    break;
                }
            }
            auto z_sets = build_flip_sets(res.flip_z_history, res.flip_z_window, res.flip_z);
            for (const auto &set : z_sets) {
                if (!set.vars || set.vars->empty()) continue;
                if (attempt_flip_pp_side("Z", *set.vars, z_checks, sz, false, res.est, ets_verbose, nullptr, set.tag)) {
                    flip_z = true;
                    break;
                }
            }
            if (flip_x || flip_z) {
                auto sx_hat = compute_syndrome(x_checks, res.est, true);
                auto sz_hat = compute_syndrome(z_checks, res.est, false);
                ok = (sx_hat == sx) && (sz_hat == sz);
                if (ok) {
                    flip_pp_ok = true;
                }
            }
            std::ostringstream oss2;
            oss2 << "[flip-PP] trial=" << trial_idx
                 << " done x=" << tf(flip_x)
                 << " z=" << tf(flip_z)
                 << " syndrome_ok=" << tf(ok);
            pp_log_line(oss2.str());
        }
        if (!ok && enable_pp) {
            const int osd_us_limit = 20;
            bool osd_x = false;
            bool osd_z = false;
            auto sx_hat = compute_syndrome(x_checks, res.est, true);
            auto sz_hat = compute_syndrome(z_checks, res.est, false);
            int osd_us_x = unsatisfied_count(sx_hat, sx);
            int osd_us_z = unsatisfied_count(sz_hat, sz);
            int osd_us = std::max(osd_us_x, osd_us_z);
            if (osd_us > osd_us_limit) {
                std::ostringstream oss;
                oss << "[osd-PP] trial=" << trial_idx
                    << " skipped us=(" << osd_us_x << "," << osd_us_z << ")";
                pp_log_line(oss.str());
            } else {
                std::ostringstream oss;
                oss << "[osd-PP] trial=" << trial_idx << " start";
                pp_log_line(oss.str());
                osd_x = attempt_osd_pp_side("X", res.abs_llr_x, x_checks, sx, true, res.est,
                                            ets_verbose, &err);
                osd_z = attempt_osd_pp_side("Z", res.abs_llr_z, z_checks, sz, false, res.est,
                                            ets_verbose, &err);
                if (osd_x || osd_z) {
                    auto sx_hat2 = compute_syndrome(x_checks, res.est, true);
                    auto sz_hat2 = compute_syndrome(z_checks, res.est, false);
                    ok = (sx_hat2 == sx) && (sz_hat2 == sz);
                    if (ok) {
                        osd_pp_ok = true;
                    }
                }
                std::ostringstream oss2;
                oss2 << "[osd-PP] trial=" << trial_idx
                     << " done x=" << tf(osd_x)
                     << " z=" << tf(osd_z)
                     << " syndrome_ok=" << tf(ok);
                pp_log_line(oss2.str());
            }
        }

        bool exact = (res.est == err);
        bool stab_equiv = false;
        bool success = ok;
        if (ok && need_basis) {
            fill_diff_bits(diff_x_bits, err, res.est, true);
            fill_diff_bits(diff_z_bits, err, res.est, false);
            stab_equiv = in_row_space(hz_basis, diff_x_bits) && in_row_space(hx_basis, diff_z_bits);
            success = ok && stab_equiv;
        }
        if (exact) exact_matches++;
        if (success && !exact && stab_equiv) {
            stab_success++;
        }
        if (success && ets_pp_ok) {
            ets_saves++;
            if (used6) ets6_saves++;
            if (used12) ets12_saves++;
        }
        if (success && flip_pp_ok) {
            flip_pp_success++;
        }
        bool osd_success = osd_pp_ok ||
                           (res.pp_success && !res.pp_success_ets && !res.pp_success_flip);
        bool pp_success_this = success && (res.pp_success || ets_pp_ok || flip_pp_ok || osd_pp_ok);
        if (pp_success_this) {
            pp_success++;
        }
        if (success && osd_success) {
            osd_pp_success++;
        }
        double latency_sec = std::chrono::duration<double>(
                                 std::chrono::steady_clock::now() - latency_start)
                                 .count();
        latency_sum_sec += latency_sec;
        latency_samples++;
        if (pp_success_this && enable_log_files) {
            if (!ensure_dir(log_dir)) {
                std::cerr << "Failed to create log directory: " << log_dir << "\n";
            } else {
                std::ostringstream pp_path;
                pp_path << log_dir << "/pp-success_p" << p_tag
                        << "_seed" << trial_seed
                        << "_trial" << trial_idx << ".txt";
                std::ofstream pp_out(pp_path.str());
                if (!pp_out) {
                    std::cerr << "Failed to write pp-success file: " << pp_path.str() << "\n";
                } else {
                    std::string pp_method;
                    if (ets_pp_ok) pp_method += (pp_method.empty() ? "ets" : ",ets");
                    if (flip_pp_ok) pp_method += (pp_method.empty() ? "flip" : ",flip");
                    if (osd_success) pp_method += (pp_method.empty() ? "osd" : ",osd");
                    pp_out << "trial=" << trial_idx << "\n";
                    pp_out << "trial_seed=" << trial_seed << "\n";
                    pp_out << "p=" << std::setprecision(8) << std::fixed << p_err << "\n";
                    pp_out << "bp_ok=" << (bp_ok ? "1" : "0") << "\n";
                    pp_out << "syndrome_match=" << (ok ? "1" : "0") << "\n";
                    pp_out << "decode_success=" << (success ? "1" : "0") << "\n";
                    pp_out << "pp_success=" << (pp_success_this ? "1" : "0") << "\n";
                    pp_out << "pp_success_in_decode=" << (res.pp_success ? "1" : "0") << "\n";
                    pp_out << "pp_success_ets=" << (ets_pp_ok ? "1" : "0") << "\n";
                    pp_out << "pp_success_flip=" << (flip_pp_ok ? "1" : "0") << "\n";
                    pp_out << "pp_success_osd=" << (osd_success ? "1" : "0") << "\n";
                    pp_out << "pp_method=" << (pp_method.empty() ? "unknown" : pp_method) << "\n";
                    pp_out << "pp_used_ets_any=" << (used_any ? "1" : "0") << "\n";
                    pp_out << "pp_used_ets6=" << (used6 ? "1" : "0") << "\n";
                    pp_out << "pp_used_ets12=" << (used12 ? "1" : "0") << "\n";
                    pp_out << "pp_log_count=" << pp_log_lines.size() << "\n";
                    pp_out << "pp_log_begin\n";
                    for (const auto &line : pp_log_lines) {
                        pp_out << line << "\n";
                    }
                    pp_out << "pp_log_end\n";
                    pp_out << "repro_cmd=./jointbp_ets --params " << params_path
                           << " --simulate --p " << std::setprecision(8) << std::fixed << p_err
                           << " --seed " << seed
                           << " --trial-index " << trial_idx
                           << " --trials 1"
                           << " --max-iter " << max_iter
                           << " --flip-hist " << flip_hist_window
                           << " --damping " << std::setprecision(6) << std::fixed << damping
                           << " --verbose 1 --verbose-all";
                    if (freeze_syn) {
                        pp_out << " --freeze-syn";
                    }
                    pp_out << "\n";
                }
            }
        }

        total_iters += res.iterations;
        if (!success) failures++;
        if (!success && enable_log_files) {
            std::ostringstream fail_path;
            fail_path << log_dir << "/fail_p" << p_tag
                      << "_seed" << trial_seed
                      << "_trial" << trial_idx << ".txt";
            std::ofstream fail_out(fail_path.str());
            if (!fail_out) {
                std::cerr << "Failed to write fail file: " << fail_path.str() << "\n";
            } else {
                bool pp_done = (!bp_ok && enable_pp);
                bool pp_ok = (!bp_ok && success && enable_pp);
                auto diff_x = diff_indices_bits(err, res.est, true);
                auto diff_z = diff_indices_bits(err, res.est, false);
                fail_out << "trial=" << trial_idx << "\n";
                fail_out << "trial_seed=" << trial_seed << "\n";
                fail_out << "p=" << std::setprecision(8) << std::fixed << p_err << "\n";
                fail_out << "bp_ok=" << (bp_ok ? "1" : "0") << "\n";
                fail_out << "syndrome_match=" << (ok ? "1" : "0") << "\n";
                fail_out << "decode_success=" << (success ? "1" : "0") << "\n";
                fail_out << "pp_performed=" << (pp_done ? "1" : "0") << "\n";
                fail_out << "pp_success=" << (pp_ok ? "1" : "0") << "\n";
                fail_out << "diff_x_size=" << diff_x.size() << "\n";
                fail_out << "diff_z_size=" << diff_z.size() << "\n";
                fail_out << "repro_cmd=./jointbp_ets --params " << params_path
                         << " --simulate --p " << std::setprecision(8) << std::fixed << p_err
                         << " --seed " << seed
                         << " --trial-index " << trial_idx
                         << " --trials 1"
                         << " --max-iter " << max_iter
                         << " --flip-hist " << flip_hist_window
                         << " --damping " << std::setprecision(6) << std::fixed << damping
                         << " --verbose 1 --verbose-all";
                if (freeze_syn) {
                    fail_out << " --freeze-syn";
                }
                fail_out << "\n";
            }
        }

        done = t + 1;
        if (report_fail && !success) {
            bool ets6 = used6;
            bool ets12 = used12;
            std::cout << "[fail] trial=" << trial_idx
                      << " iter=" << res.iterations
                      << " bp_ok=" << tf(bp_ok)
                      << " syndrome_match=" << tf(ok)
                      << " ets6=" << tf(ets6)
                      << " ets12=" << tf(ets12)
                      << " trial_seed=" << trial_seed
                      << "\n";
            auto sx_hat = compute_syndrome(x_checks, res.est, true);
            auto sz_hat = compute_syndrome(z_checks, res.est, false);
            auto us_x = unsatisfied_indices(sx_hat, sx);
            auto us_z = unsatisfied_indices(sz_hat, sz);
            auto diff_x = diff_indices_bits(err, res.est, true);
            auto diff_z = diff_indices_bits(err, res.est, false);
            bool syn_x = (sx_hat == sx);
            bool syn_z = (sz_hat == sz);
            if (!syn_x) {
                print_index_list("us_x", us_x);
                print_index_list("diff_x", diff_x);
            }
            if (!syn_z) {
                print_index_list("us_z", us_z);
                print_index_list("diff_z", diff_z);
            }
        }
        if (!success) {
            bool redo_print = report_fail && !iter_verbose;
            if (redo_print) {
                std::cout << "[redo] trial=" << trial_idx << " start\n";
                auto redo = joint_bp_decode(
                    x_checks, z_checks, var_to_x, var_to_z,
                    sx, sz, prior, max_iter, flip_hist_window, freeze_syn, damping,
                    redo_print, verbose_all, false, nullptr, &err,
                    need_basis ? &hx_basis : nullptr, need_basis ? &hz_basis : nullptr, nullptr,
                    nullptr
                );
                std::cout << "[redo] trial=" << trial_idx
                          << " done iter=" << redo.iterations
                          << " syndrome_match=" << tf(redo.syndrome_match)
                          << "\n";
            }
        }
        if (report_every > 0 && (done % report_every == 0)) {
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - start).count();
            double avg_latency_sec = (latency_samples > 0) ? (latency_sum_sec / latency_samples) : 0.0;
            report_progress(done, failures, bp_failures, pp_success, ets_saves, flip_pp_success,
                            osd_pp_success, stab_success, ets_used, ets_labels, ets_used_counts,
                            total_iters, elapsed, k, avg_latency_sec);
        }
        if (progress_every > 0 && (done % progress_every == 0)) {
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - start).count();
            write_progress_tsv(done, failures, pp_success, ets_saves, flip_pp_success, stab_success, total_iters, elapsed);
            last_progress_logged = done;
        }
        if (min_fails > 0 && failures >= min_fails) {
            break;
        }
    }

    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();
    if (done > 0 && done != last_progress_logged) {
        write_progress_tsv(done, failures, pp_success, ets_saves, flip_pp_success, stab_success, total_iters, elapsed);
    }
    double fer = (done > 0) ? (static_cast<double>(failures) / static_cast<double>(done)) : 0.0;
    double ci_lo = 0.0;
    double ci_hi = 0.0;
    bool has_ci = wilson_interval(done, failures, 1.96, ci_lo, ci_hi);
    double avg_iter = (done > 0) ? (static_cast<double>(total_iters) / static_cast<double>(done)) : 0.0;
    double avg_latency_sec = (latency_samples > 0) ? (latency_sum_sec / latency_samples) : 0.0;
    double exact_rate = (done > 0) ? (static_cast<double>(exact_matches) / static_cast<double>(done)) : 0.0;

    print_stdout_section("Run Summary");
    std::cout << "trials=" << done << " failures=" << failures
              << " bp_fail=" << bp_failures
              << " pp_success=" << pp_success
              << " flip_pp_success=" << flip_pp_success
              << " stab_success=" << stab_success
              << " ets_used=" << ets_used
              << " ets_saved=" << ets_saves
              << " ets6_used=" << ets6_used
              << " ets6_saved=" << ets6_saves
              << " ets12_used=" << ets12_used
              << " ets12_saved=" << ets12_saves
              << " fer=" << std::setprecision(8) << std::fixed << fer;
    if (has_ci) {
        std::cout << " ci95=[" << std::setprecision(8) << std::fixed << ci_lo
                  << "," << std::setprecision(8) << std::fixed << ci_hi << "]";
    }
    std::cout << " iters=" << total_iters
              << " avg_iter=" << std::setprecision(2) << std::fixed << avg_iter
              << " avg_latency=" << format_latency(avg_latency_sec)
              << " exact_rate=" << std::setprecision(6) << std::fixed << exact_rate
              << " elapsed_s=" << std::setprecision(2) << std::fixed << elapsed << "s"
              << "\n";
    return 0;
}
