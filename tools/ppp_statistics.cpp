// Standalone C++17 statistics for RTKLIB calendar-time ECEF / decimal-LLH .pos.
// All time differences use the time system printed in the file (normally GPST).
// Build: cmake --build . --config Release --target ppp_statistics
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif
namespace fs = std::filesystem;
using Vec3 = std::array<double, 3>;
constexpr double PI = 3.14159265358979323846;
constexpr double A = 6378137.0;
constexpr double E2 = 6.6943799901413165e-3;
constexpr double EPS = 1e-6;
const double NA = std::numeric_limits<double>::quiet_NaN();
struct Epoch {
    double time{};
    Vec3 xyz{}, enu{};
    int quality{}, satellites{};
    double h2{}, u2{};
};
struct Options {
    fs::path input, output;
    std::string label = "solution", format = "auto", start;
    double referenceSeconds = 180, windowSeconds = 300;
    double horizontal = 0.10, vertical = 0.20, maxGap = 0;
    int quality = 6;
};
struct Reference { Vec3 xyz{}, llh{}; size_t begin{}, count{}; };
struct Window {
    size_t begin{}, end{}; // end is exclusive
    double stop{}, hrms{}, urms{};
    bool valid{}, pass{};
};
struct Convergence {
    std::vector<Window> windows;
    size_t index = std::numeric_limits<size_t>::max();
    bool found() const { return index != std::numeric_limits<size_t>::max(); }
};
struct Metrics {
    size_t count{}, speedPairs{};
    double hrms = NA, urms = NA, speed = NA;
    std::array<double, 3> hq{NA, NA, NA}, uq{NA, NA, NA};
};
// Gregorian calendar arithmetic, independent of local timezone or daylight saving.
long long daysFromCivil(int y, int m, int d) {
    y -= m <= 2;
    const int era = (y >= 0 ? y : y - 399) / 400;
    const unsigned yo = static_cast<unsigned>(y - era * 400);
    const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
    const unsigned doe = yo * 365 + yo / 4 - yo / 100 + doy;
    return static_cast<long long>(era) * 146097 + doe - 719468;
}
std::array<int, 3> civilFromDays(long long z) {
    z += 719468;
    const long long era = (z >= 0 ? z : z - 146096) / 146097;
    const unsigned doe = static_cast<unsigned>(z - era * 146097);
    const unsigned yo = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    int y = static_cast<int>(yo + era * 400);
    const unsigned doy = doe - (365 * yo + yo / 4 - yo / 100);
    const unsigned mp = (5 * doy + 2) / 153;
    const int d = static_cast<int>(doy - (153 * mp + 2) / 5 + 1);
    const int m = static_cast<int>(mp) + (mp < 10 ? 3 : -9);
    y += m <= 2;
    return {y, m, d};
}
double parseTime(const std::string& date, const std::string& clock) {
    int y, m, d, hh, mm, consumed = 0;
    double ss;
    if (std::sscanf(date.c_str(), "%d/%d/%d%n", &y, &m, &d, &consumed) != 3 ||
        consumed != static_cast<int>(date.size()) || y < 1970 || y > 2200 ||
        m < 1 || m > 12 || d < 1 || d > 31)
        throw std::runtime_error("Invalid calendar date: " + date);
    consumed = 0;
    if (std::sscanf(clock.c_str(), "%d:%d:%lf%n", &hh, &mm, &ss, &consumed) != 3 ||
        consumed != static_cast<int>(clock.size()) || hh < 0 || hh > 23 ||
        mm < 0 || mm > 59 || !std::isfinite(ss) || ss < 0 || ss >= 60)
        throw std::runtime_error("Invalid clock: " + clock);
    const long long days = daysFromCivil(y, m, d);
    if (civilFromDays(days) != std::array<int, 3>{y, m, d})
        throw std::runtime_error("Invalid calendar date: " + date);
    return days * 86400.0 + hh * 3600.0 + mm * 60.0 + ss;
}
double parseTime(const std::string& str) {
    std::istringstream in(str);
    std::string date, clock, extra;
    if (!(in >> date >> clock) || (in >> extra))
        throw std::runtime_error("Use --start \"YYYY/MM/DD hh:mm:ss.sss\"");
    return parseTime(date, clock);
}
std::string timeText(double time) {
    const long long ms = std::llround(time * 1000), day = ms / 86400000, rem = ms % 86400000;
    const auto c = civilFromDays(day);
    char b[64];
    std::snprintf(b, sizeof(b), "%04d/%02d/%02d %02d:%02d:%06.3f",
                  c[0], c[1], c[2], static_cast<int>(rem / 3600000),
                  static_cast<int>(rem / 60000 % 60), (rem % 60000) / 1000.0);
    return b;
}
Vec3 toXyz(const Vec3& llh) {
    const double lat = llh[0] * PI / 180, lon = llh[1] * PI / 180;
    const double n = A / std::sqrt(1 - E2 * std::pow(std::sin(lat), 2));
    return {(n + llh[2]) * std::cos(lat) * std::cos(lon),
            (n + llh[2]) * std::cos(lat) * std::sin(lon),
            (n * (1 - E2) + llh[2]) * std::sin(lat)};
}
Vec3 toLlh(const Vec3& xyz) {
    const double p = std::hypot(xyz[0], xyz[1]);
    if (p < 1e-8)
        return {std::copysign(90.0, xyz[2]), 0, std::abs(xyz[2]) - A * std::sqrt(1 - E2)};
    double lat = std::atan2(xyz[2], p * (1 - E2)), n = A;
    for (int i = 0; i < 20; ++i) {
        n = A / std::sqrt(1 - E2 * std::pow(std::sin(lat), 2));
        const double next = std::atan2(xyz[2] + E2 * n * std::sin(lat), p);
        if (std::abs(next - lat) < 1e-14) { lat = next; break; }
        lat = next;
    }
    n = A / std::sqrt(1 - E2 * std::pow(std::sin(lat), 2));
    const double height = p * std::cos(lat) + xyz[2] * std::sin(lat) -
                          n * (1 - E2 * std::pow(std::sin(lat), 2));
    return {lat * 180 / PI, std::atan2(xyz[1], xyz[0]) * 180 / PI, height};
}
Vec3 toEnu(const Vec3& xyz, const Reference& ref) {
    const double lat = ref.llh[0] * PI / 180, lon = ref.llh[1] * PI / 180;
    const double sl = std::sin(lat), cl = std::cos(lat), so = std::sin(lon), co = std::cos(lon);
    const double x = xyz[0] - ref.xyz[0], y = xyz[1] - ref.xyz[1], z = xyz[2] - ref.xyz[2];
    return {-so * x + co * y, -sl * co * x - sl * so * y + cl * z,
            cl * co * x + cl * so * y + sl * z};
}
std::vector<Epoch> readPos(std::istream& in, std::string format, std::string& timeSystem) {
    const bool automatic = format == "auto";
    std::vector<Epoch> epochs;
    std::string line;
    size_t lineNumber = 0;
    while (std::getline(in, line)) {
        ++lineNumber;
        if (lineNumber == 1 && line.compare(0, 3, "\xEF\xBB\xBF") == 0) line.erase(0, 3);
        const auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos) continue;
        if (line[first] == '%' || line[first] == '#') {
            if (automatic) {
                if (line.find("x-ecef") != std::string::npos) format = "ecef";
                if (line.find("latitude(deg)") != std::string::npos) format = "llh";
            }
            for (const std::string sys : {"GPST", "UTC", "JST"})
                if (line.find(sys) != std::string::npos) timeSystem = sys;
            if (line.find("latitude(d") != std::string::npos &&
                line.find("latitude(deg)") == std::string::npos)
                throw std::runtime_error("DMS coordinates are unsupported; export decimal degrees or ECEF.");
            continue;
        }
        if (format != "ecef" && format != "llh")
            throw std::runtime_error("No coordinate header; specify --format ecef or llh.");
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        std::string date, clock;
        Epoch e;
        if (!(row >> date >> clock >> e.xyz[0] >> e.xyz[1] >> e.xyz[2] >> e.quality >> e.satellites))
            throw std::runtime_error("Malformed .pos record at line " + std::to_string(lineNumber));
        e.time = parseTime(date, clock);
        if (e.quality < 0 || e.quality > 7 || e.satellites < 0 ||
            !std::all_of(e.xyz.begin(), e.xyz.end(), [](double x) { return std::isfinite(x); }))
            throw std::runtime_error("Invalid solution values at line " + std::to_string(lineNumber));
        if (!epochs.empty() && e.time <= epochs.back().time)
            throw std::runtime_error("Timestamps must be strictly increasing (no duplicate epochs).");
        if (format == "llh") {
            if (std::abs(e.xyz[0]) > 90 || std::abs(e.xyz[1]) > 180)
                throw std::runtime_error("Invalid decimal latitude or longitude.");
            e.xyz = toXyz(e.xyz);
        }
        if (e.quality > 0) {
            const double radius = std::hypot(std::hypot(e.xyz[0], e.xyz[1]), e.xyz[2]);
            if (radius < 6e6 || radius > 7e6)
                throw std::runtime_error("Position is outside the supported near-Earth surface range.");
        }
        epochs.push_back(e);
    }
    if (epochs.empty()) throw std::runtime_error("No solution epochs in input.");
    return epochs;
}
double quantile(std::vector<double> values, double p) {
    if (values.empty()) return NA;
    std::sort(values.begin(), values.end());
    const double k = p * (values.size() - 1);
    const size_t i = static_cast<size_t>(k), j = std::min(i + 1, values.size() - 1);
    return values[i] + (values[j] - values[i]) * (k - i);
}
double sampleInterval(const std::vector<Epoch>& rows) {
    std::vector<double> dt;
    for (size_t i = 1; i < rows.size(); ++i) dt.push_back(rows[i].time - rows[i - 1].time);
    if (dt.empty()) throw std::runtime_error("At least two target-quality epochs are required.");
    return quantile(dt, 0.5);
}
Reference referenceMean(const std::vector<Epoch>& rows, double seconds, double dt, double maxGap) {
    Reference r;
    const double boundary = rows.back().time - seconds;
    while (r.begin < rows.size() && rows[r.begin].time <= boundary + EPS) ++r.begin;
    r.count = rows.size() - r.begin;
    if (r.count < 2 || rows[r.begin].time > boundary + dt + EPS)
        throw std::runtime_error("Insufficient data coverage in the final reference interval.");
    for (size_t i = r.begin + 1; i < rows.size(); ++i)
        if (rows[i].time - rows[i - 1].time > maxGap + EPS)
            throw std::runtime_error("Gap in final reference interval; cannot form the requested truth mean.");
    for (size_t k = 0; k < 3; ++k) {
        long double sum = 0;
        for (size_t i = r.begin; i < rows.size(); ++i) sum += rows[i].xyz[k];
        r.xyz[k] = static_cast<double>(sum / r.count);
    }
    r.llh = toLlh(r.xyz);
    return r;
}
void setErrors(std::vector<Epoch>& rows, const Reference& ref) {
    for (auto& e : rows) {
        e.enu = toEnu(e.xyz, ref);
        e.h2 = e.enu[0] * e.enu[0] + e.enu[1] * e.enu[1];
        e.u2 = e.enu[2] * e.enu[2];
    }
}
Convergence convergence(const std::vector<Epoch>& rows, const Options& o, double dt, double maxGap) {
    Convergence out;
    const size_t n = rows.size();
    std::vector<long double> hs(n + 1), us(n + 1);
    std::vector<size_t> gaps(n);
    for (size_t i = 0; i < n; ++i) {
        hs[i + 1] = hs[i] + rows[i].h2;
        us[i + 1] = us[i] + rows[i].u2;
        if (i) gaps[i] = gaps[i - 1] + (rows[i].time - rows[i - 1].time > maxGap + EPS);
    }
    size_t j = 0;
    for (size_t i = 0; i < n; ++i) {
        const double end = rows[i].time + o.windowSeconds;
        // 300 regular one-second epochs cover [t, t+300), through t+299.
        if (end > rows.back().time + dt + EPS) break;
        while (j < n && rows[j].time < end - EPS) ++j;
        const size_t count = j - i;
        Window w;
        w.begin = i; w.end = j; w.stop = end;
        w.hrms = std::sqrt(static_cast<double>((hs[j] - hs[i]) / count));
        w.urms = std::sqrt(static_cast<double>((us[j] - us[i]) / count));
        w.valid = count >= 2 && gaps[j - 1] == gaps[i] &&
                  rows[j - 1].time + maxGap >= end - EPS;
        // Roundoff guard only (1 picometre), preserving inclusive <= at the boundary.
        w.pass = w.valid && w.hrms <= o.horizontal + 1e-12 && w.urms <= o.vertical + 1e-12;
        out.windows.push_back(w);
    }
    // A temporary crossing is not convergence: EVERY later full window must pass.
    bool sustained = true;
    for (size_t k = out.windows.size(); k-- > 0;) {
        sustained = sustained && out.windows[k].pass;
        if (sustained) out.index = out.windows[k].begin;
    }
    return out;
}
Metrics metrics(const std::vector<Epoch>& all, int quality, double begin, double maxGap) {
    Metrics m;
    long double h2 = 0, u2 = 0, v2 = 0;
    std::vector<double> h, u;
    for (size_t i = 0; i < all.size(); ++i) {
        const auto& e = all[i];
        if (e.quality != quality || e.time < begin - EPS) continue;
        ++m.count;
        if (!quality) continue;
        h2 += e.h2; u2 += e.u2;
        h.push_back(std::sqrt(e.h2)); u.push_back(std::sqrt(e.u2));
        if (i && all[i - 1].quality == quality && all[i - 1].time >= begin - EPS) {
            const double dt = e.time - all[i - 1].time;
            if (dt <= maxGap + EPS) {
                for (size_t k = 0; k < 3; ++k)
                    v2 += std::pow((e.xyz[k] - all[i - 1].xyz[k]) / dt, 2);
                ++m.speedPairs;
            }
        }
    }
    if (!m.count || !quality) return m;
    m.hrms = std::sqrt(static_cast<double>(h2 / m.count));
    m.urms = std::sqrt(static_cast<double>(u2 / m.count));
    if (m.speedPairs) m.speed = std::sqrt(static_cast<double>(v2 / m.speedPairs));
    const double p[] = {0.68, 0.95, 0.997};
    for (size_t k = 0; k < 3; ++k) { m.hq[k] = quantile(h, p[k]); m.uq[k] = quantile(u, p[k]); }
    return m;
}
std::string number(double x, int digits = 6) {
    if (!std::isfinite(x)) return "N/A";
    std::ostringstream out;
    out << std::fixed << std::setprecision(digits) << x;
    return out.str();
}
std::string csvCell(const std::string& s) {
    std::string out = "\"";
    for (const char c : s) { out += c; if (c == '"') out += '"'; }
    return out + '"';
}
std::string qualityName(int q) {
    const char* names[] = {"无效解", "固定解", "浮点解", "SBAS", "差分解", "单点定位", "精密单点 PPP", "DR"};
    return names[q];
}
std::ofstream openOutput(const fs::path& p, bool csv = false) {
    std::ofstream out(p, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot write output: " + p.u8string());
    out.exceptions(std::ios::badbit | std::ios::failbit);
    if (csv) out << "\xEF\xBB\xBF"; // Excel recognizes Chinese UTF-8 headers.
    return out;
}

void outputs(const Options& o, const std::vector<Epoch>& all, const std::vector<Epoch>& rows,
             const Reference& ref, const Convergence& cv, double dt, double maxGap,
             const std::string& timeSystem) {
    if (!o.output.parent_path().empty()) fs::create_directories(o.output.parent_path());
    const auto path = [&](const std::string& suffix) { fs::path p = o.output; p += suffix; return p; };
    auto md = openOutput(path(".md"));
    auto csv = openOutput(path(".csv"), true);
    auto windows = openOutput(path("_windows.csv"), true);
    auto epochs = openOutput(path("_epochs.csv"), true);
    const double start = o.start.empty() ? all.front().time : parseTime(o.start);
    const double tc = cv.found() ? rows[cv.index].time : NA;
    const double elapsed = cv.found() ? tc - start : NA;
    std::ostringstream report;
    report << "# " << o.label << " 定位精度与收敛统计\n\n"
           << "- 输入文件：" << o.input.u8string() << "\n"
           << "- 时间系统：" << timeSystem << "；统计解类型：Q=" << o.quality << "。\n"
           << "- 解算文件：" << timeText(all.front().time) << " 至 " << timeText(all.back().time)
           << "，共 " << all.size() << " 个历元。\n"
           << "- 参考真值：最后 " << number(o.referenceSeconds, 0) << " 秒目标解的 ECEF 算术平均值，区间 ("
           << timeText(rows.back().time - o.referenceSeconds) << ", " << timeText(rows.back().time)
           << "]，实际 " << ref.count << " 个历元（" << timeText(rows[ref.begin].time)
           << " 至 " << timeText(rows.back().time) << "）。\n"
           << "- 参考 ECEF (m)：X=" << number(ref.xyz[0], 9) << "，Y=" << number(ref.xyz[1], 9)
           << "，Z=" << number(ref.xyz[2], 9) << "。\n"
           << "- 参考纬度/经度/椭球高：" << number(ref.llh[0], 10) << "° / "
           << number(ref.llh[1], 10) << "° / " << number(ref.llh[2], 6) << " m。\n"
           << "- 收敛判据：每个目标解历元 t 检查 [t, t+" << number(o.windowSeconds, 0)
           << " s) 窗口，水平 RMS≤" << number(o.horizontal * 100, 2)
           << " cm 且高程 RMS≤" << number(o.vertical * 100, 2)
           << " cm；首个此后所有完整窗口持续达标的窗口起点为收敛时刻。\n"
           << "- 中位采样间隔 " << number(dt, 3) << " s；间隔超过 " << number(maxGap, 3)
           << " s 的窗口无效。尾部不足一个窗口不参与判定；收敛后总体 RMS 仍使用直到末历元的全部目标解。\n"
           << "- 指定计时起点：" << timeText(start)
           << (o.start.empty() ? "（解算文件首历元，未提供采集起点）。\n" : "（采集首个观测历元）。\n")
           << "\n";
    if (cv.found()) {
        const auto& w = cv.windows.at(cv.index); // one candidate per target epoch
        report << "**收敛时刻：" << timeText(tc) << " " << timeSystem << "。**\n\n"
               << "- 从指定计时起点：**" << number(elapsed, 3) << " s（"
               << number(elapsed / 60, 6) << " min）**。\n"
               << "- 从解算文件首历元：" << number(tc - all.front().time, 3) << " s。\n"
               << "- 从首个目标质量解：" << number(tc - rows.front().time, 3) << " s。\n"
               << "- 首个达标窗口：" << timeText(tc) << " 至 " << timeText(w.stop)
               << "（右端不含），" << w.end - w.begin << " 个历元；水平 RMS="
               << number(w.hrms, 9) << " m，高程 RMS=" << number(w.urms, 9) << " m。\n"
               << "- 首窗口完成边界：" << timeText(w.stop)
               << "；上述持续达标结论由整段数据事后确认。\n";
        if (cv.index) {
            const auto& prev = cv.windows[cv.index - 1];
            report << "- 前一个候选窗口：" << timeText(rows[prev.begin].time) << "，水平 RMS="
                   << number(prev.hrms, 9) << " m，高程 RMS=" << number(prev.urms, 9)
                   << " m，" << (prev.valid ? "至少一项超过阈值" : "数据不连续") << "。\n";
        }
    } else report << "**未收敛：没有满足持续达标要求的完整窗口。**\n";
    report << "\n统计公式：将 ECEF 坐标差旋转到参考点固定 ENU 坐标系，"
           << "水平 RMS = sqrt(mean(E²+N²))，高程 RMS = sqrt(mean(U²))；"
           << "相对固定真值计算，不再减去各窗口均值。\n\n"
           << "水平 CEP 为水平误差模长的经验分位数，高程 P 为 |U| 的经验分位数，"
           << "采用排序后 (n−1)p 位置线性插值。分位数列依次为 68%、95%、99.7%。"
           << "坐标差分速度 RMS 为相邻同质量解的三维坐标差/时间差的模长 RMS（非多普勒速度），"
           << "不跨越数据间断。各行占比的分母为对应统计区间内解算文件全部历元，含无效解；"
           << "没有数据的指标记为 N/A。参考真值按本次指定规则构造，属于内部一致性统计。\n\n";
    csv << "子目录,统计区间,定位质量,Q,数据数量,占区间总数据(%),水平RMS(m),高程RMS(m),"
        << "水平CEP68(m),水平CEP95(m),水平CEP99.7(m),高程P68(m),高程P95(m),高程P99.7(m),"
        << "坐标差分速度RMS(m/s),速度样本数,收敛时刻(" << timeSystem
        << "),收敛时间_指定起点(s),收敛时间_首解(s),窗口(s),参考区间(s)\r\n";
    std::map<int, bool> qualities{{5, true}, {o.quality, true}};
    for (const auto& e : all) qualities[e.quality] = true;
    for (int scope = 0; scope < 2; ++scope) {
        if (scope && !cv.found()) continue;
        const double begin = scope ? tc : all.front().time;
        const std::string scopeName = scope ? "收敛后" : "全时段";
        size_t denominator = 0;
        for (const auto& e : all) if (e.time >= begin - EPS) ++denominator;
        report << "## " << scopeName << "\n\n"
               << "| 子目录 | 定位质量 | 数据数量 | 占比 | 水平RMS(m) | 高程RMS(m) | "
               << "水平CEP68/95/99.7(m) | 高程P68/95/99.7(m) | 坐标差分速度RMS(m/s) | 收敛时间(s) |\n"
               << "|---|---|---:|---:|---:|---:|---|---|---:|---:|\n";
        for (const auto& q : qualities) {
            const Metrics m = metrics(all, q.first, begin, maxGap);
            const double percent = denominator ? m.count * 100.0 / denominator : 0;
            const bool target = q.first == o.quality;
            csv << csvCell(o.label) << ',' << scopeName << ',' << qualityName(q.first) << ',' << q.first
                << ',' << m.count << ',' << number(percent, 2) << ',' << number(m.hrms) << ',' << number(m.urms);
            for (const double x : m.hq) csv << ',' << number(x);
            for (const double x : m.uq) csv << ',' << number(x);
            csv << ',' << number(m.speed) << ',' << m.speedPairs << ','
                << (target && cv.found() ? timeText(tc) : "N/A") << ','
                << number(target ? elapsed : NA, 3) << ','
                << number(target && cv.found() ? tc - all.front().time : NA, 3)
                << ',' << number(o.windowSeconds, 3) << ',' << number(o.referenceSeconds, 3) << "\r\n";
            report << "| " << o.label << " | " << qualityName(q.first) << " | " << m.count << " | "
                   << number(percent, 2) << "% | " << number(m.hrms) << " | " << number(m.urms) << " | "
                   << number(m.hq[0]) << " / " << number(m.hq[1]) << " / " << number(m.hq[2]) << " | "
                   << number(m.uq[0]) << " / " << number(m.uq[1]) << " / " << number(m.uq[2]) << " | "
                   << number(m.speed) << " | " << number(target ? elapsed : NA, 3) << " |\n";
        }
        report << '\n';
    }
    windows << "窗口起点(" << timeSystem << "),窗口终点_不含(" << timeSystem
            << "),历元数,窗口完整且连续,水平RMS(m),高程RMS(m),当前窗口达标,从此后所有完整窗口达标\r\n";
    for (const auto& w : cv.windows)
        windows << timeText(rows[w.begin].time) << ',' << timeText(w.stop) << ',' << w.end - w.begin
                << ',' << w.valid << ',' << number(w.hrms, 9) << ',' << number(w.urms, 9)
                << ',' << w.pass << ',' << (cv.found() && w.begin >= cv.index) << "\r\n";
    epochs << "时间(" << timeSystem << "),Q,卫星数,E(m),N(m),U(m),水平误差(m),高程绝对误差(m),收敛后,用于参考均值\r\n";
    for (const auto& e : all)
        epochs << timeText(e.time) << ',' << e.quality << ',' << e.satellites << ','
               << number(e.quality ? e.enu[0] : NA, 9) << ',' << number(e.quality ? e.enu[1] : NA, 9) << ','
               << number(e.quality ? e.enu[2] : NA, 9) << ',' << number(e.quality ? std::sqrt(e.h2) : NA, 9) << ','
               << number(e.quality ? std::sqrt(e.u2) : NA, 9) << ',' << (cv.found() && e.time >= tc - EPS)
               << ',' << (e.quality == o.quality && e.time > rows.back().time - o.referenceSeconds + EPS) << "\r\n";
    md << report.str();
    md.close(); csv.close(); windows.close(); epochs.close();
    std::cout << report.str() << "\nOutputs: " << path(".md").u8string() << ", "
              << path(".csv").u8string() << ", *_windows.csv, *_epochs.csv\n";
}
int selfTest() {
    int checks = 0;
    const auto check = [&](bool ok, const char* name) {
        ++checks;
        if (!ok) throw std::runtime_error(std::string("Self-test failed: ") + name);
    };
    const auto rejects = [&](auto fn, const char* name) {
        bool rejected = false;
        try { fn(); } catch (const std::exception&) { rejected = true; }
        check(rejected, name);
    };
    const double t = parseTime("2026/09/09", "23:59:59.500");
    check(timeText(t + 1) == "2026/09/10 00:00:00.500", "midnight and milliseconds");
    check(parseTime("2024/03/01 00:00:00") - parseTime("2024/02/28 00:00:00") == 172800, "leap year");
    rejects([] { parseTime("2026/02/29 00:00:00"); }, "invalid calendar date");
    rejects([] { parseTime("2026/09/09 03:40:60"); }, "invalid seconds");
    const Vec3 original{31.474066, 120.264656, 31.2};
    const Vec3 back = toLlh(toXyz(original));
    for (int k = 0; k < 3; ++k) check(std::abs(back[k] - original[k]) < 1e-7, "LLH round trip");
    check(std::abs(toLlh(toXyz({90, 0, 7}))[2] - 7) < 1e-7, "polar height");
    Reference ref; ref.xyz = {A, 0, 0}; ref.llh = {0, 0, 0};
    check(toEnu({A + 3, 1, 2}, ref) == Vec3{1, 2, 3}, "known ENU axes");
    std::vector<Epoch> a(601);
    for (size_t i = 0; i < a.size(); ++i) {
        a[i].time = static_cast<double>(i);
        a[i].quality = 6; a[i].xyz = {A, static_cast<double>(i), 0};
    }
    const Reference mean = referenceMean(a, 180, 1, 1.5);
    check(mean.count == 180 && mean.begin == 421, "reference interval excludes left endpoint");
    check(mean.xyz[1] == 510.5, "reference arithmetic mean");
    auto gap = a; gap.erase(gap.begin() + 500);
    rejects([&] { referenceMean(gap, 180, 1, 1.5); }, "reference gap");
    rejects([&] { referenceMean(std::vector<Epoch>(a.begin(), a.begin() + 30), 180, 1, 1.5); },
            "short reference interval");
    std::vector<Epoch> pair(2);
    pair[0].time = 0; pair[1].time = 1;
    pair[0].xyz = {A, 0.06, 0.08}; pair[1].xyz = {A + 0.1, 0.06, 0.08};
    pair[0].quality = pair[1].quality = 6;
    setErrors(pair, ref);
    const auto mm = metrics(pair, 6, 0, 1.5);
    check(std::abs(mm.hrms - 0.1) < 1e-12, "horizontal vector RMS");
    check(std::abs(mm.urms - std::sqrt(0.005)) < 1e-8, "vertical RMS retains bias");
    check(std::abs(mm.speed - 0.1) < 1e-8, "coordinate difference speed");
    check(!std::isfinite(metrics(pair, 5, 0, 1.5).hrms), "empty quality is N/A");
    check(quantile({1, 2, 3, 4}, 0.5) == 2.5, "linear empirical quantile");
    check(std::abs(quantile({1, 2, 3, 4}, 0.997) - 3.991) < 1e-12, "99.7 percentile");
    Options o; o.windowSeconds = 3; o.horizontal = 0.1; o.vertical = 0.2;
    std::vector<Epoch> synthetic(15);
    for (size_t i = 0; i < synthetic.size(); ++i) {
        synthetic[i].time = static_cast<double>(i);
        synthetic[i].h2 = i == 5 ? 1 : 0;
    }
    auto cv = convergence(synthetic, o, 1, 1.5);
    check(cv.windows.size() == 13 && cv.windows.back().end == 15, "full half-open window counting");
    check(cv.windows.front().pass && cv.index == 6, "later failure cancels temporary convergence");
    check(cv.windows[6].stop == 9, "window start differs from confirmation boundary");
    for (auto& e : synthetic) { e.h2 = 0; e.u2 = 1; }
    check(!convergence(synthetic, o, 1, 1.5).found(), "vertical threshold also required");
    for (auto& e : synthetic) { e.h2 = o.horizontal * o.horizontal; e.u2 = o.vertical * o.vertical; }
    check(convergence(synthetic, o, 1, 1.5).index == 0, "inclusive threshold equality");
    for (auto& e : synthetic) e.h2 = std::pow(o.horizontal + 1e-8, 2);
    check(!convergence(synthetic, o, 1, 1.5).found(), "values above limit still fail");
    synthetic.resize(2);
    check(!convergence(synthetic, o, 1, 1.5).found(), "insufficient full window");
    synthetic.resize(10);
    for (size_t i = 0; i < synthetic.size(); ++i) {
        synthetic[i].time = static_cast<double>(i + (i >= 5 ? 5 : 0));
        synthetic[i].h2 = synthetic[i].u2 = 0;
    }
    check(convergence(synthetic, o, 1, 1.5).index == 5, "gap breaks sustained convergence");
    std::string system = "unknown";
    std::istringstream pos("% GPST x-ecef(m)\n2026/09/09 03:41:17.000 6378137 0 0 6 10\n"
                           "2026/09/09 03:41:18.000 6378137 1 0 5 8\n");
    const auto parsed = readPos(pos, "auto", system);
    check(parsed.size() == 2 && parsed[1].quality == 5 && system == "GPST", "mixed-quality ECEF parser");
    std::istringstream llh("% GPST latitude(deg)\n2026/09/09 03:41:17.000,0,0,0,6,10\n");
    check(readPos(llh, "auto", system)[0].xyz == Vec3{A, 0, 0}, "decimal LLH and comma separator");
    rejects([&] {
        std::istringstream duplicate("2026/09/09 00:00:00 6378137 0 0 6 8\n"
                                     "2026/09/09 00:00:00 6378137 0 0 6 8\n");
        readPos(duplicate, "ecef", system);
    }, "duplicate timestamps rejected");
    rejects([&] {
        std::istringstream bad("2026/09/09 00:00:00 0 0 0 6 8\n");
        readPos(bad, "ecef", system);
    }, "invalid Earth position rejected");
    std::cout << "ppp_statistics: " << checks << " checks passed.\n";
    return 0;
}
void usage() {
    std::cout << "ppp_statistics --input solution.pos --output result_prefix [options]\n"
              << "  --label NAME             Table dataset label\n"
              << "  --start \"YYYY/MM/DD hh:mm:ss.sss\"  Acquisition start, same time system as .pos\n"
              << "  --reference-seconds 180  Final ECEF mean interval (t_end-R, t_end]\n"
              << "  --window-seconds 300     Sustained forward RMS windows [t, t+W)\n"
              << "  --horizontal 0.10        Horizontal RMS limit (m)\n"
              << "  --vertical 0.20          Vertical RMS limit (m)\n"
              << "  --quality 6              Reference/convergence quality (RTKLIB Q=6 PPP)\n"
              << "  --format auto|ecef|llh    Decimal LLH or ECEF, calendar timestamps\n"
              << "  --max-gap SECONDS        Default: 1.5 * median target sampling interval\n"
              << "  --self-test              Run analytic regression cases\n"
              << "Outputs: .md, .csv, _windows.csv, _epochs.csv (CSV UTF-8 with BOM).\n";
}
double positiveNumber(const std::string& s) {
    size_t count = 0;
    const double x = std::stod(s, &count);
    if (count != s.size() || !std::isfinite(x) || x <= 0)
        throw std::runtime_error("Expected positive finite number: " + s);
    return x;
}
int run(const std::vector<std::string>& args) {
    if (args.size() == 1 || (args.size() == 2 && args[1] == "--help")) { usage(); return 0; }
    if (args.size() == 2 && args[1] == "--self-test") return selfTest();
    Options o;
    for (size_t i = 1; i < args.size(); ++i) {
        const std::string key = args[i];
        if (++i == args.size()) throw std::runtime_error("Missing value for " + key);
        const std::string value = args[i];
        if (key == "--input") o.input = fs::u8path(value);
        else if (key == "--output") o.output = fs::u8path(value);
        else if (key == "--label") o.label = value;
        else if (key == "--start") o.start = value;
        else if (key == "--format") o.format = value;
        else if (key == "--reference-seconds") o.referenceSeconds = positiveNumber(value);
        else if (key == "--window-seconds") o.windowSeconds = positiveNumber(value);
        else if (key == "--horizontal") o.horizontal = positiveNumber(value);
        else if (key == "--vertical") o.vertical = positiveNumber(value);
        else if (key == "--max-gap") o.maxGap = positiveNumber(value);
        else if (key == "--quality") {
            const double q = positiveNumber(value);
            if (q != std::floor(q) || q > 7) throw std::runtime_error("Quality must be an integer 1..7.");
            o.quality = static_cast<int>(q);
        } else throw std::runtime_error("Unknown option: " + key);
    }
    if (o.input.empty() || o.output.empty()) throw std::runtime_error("--input and --output are required.");
    if (o.format != "auto" && o.format != "ecef" && o.format != "llh")
        throw std::runtime_error("Format must be auto, ecef, or llh.");
    for (const std::string suffix : {".md", ".csv", "_windows.csv", "_epochs.csv"}) {
        fs::path out = o.output; out += suffix;
        if (fs::exists(out) && fs::equivalent(o.input, out))
            throw std::runtime_error("Output must not overwrite the input solution.");
    }
    std::ifstream input(o.input);
    if (!input) throw std::runtime_error("Cannot read input: " + o.input.u8string());
    std::string timeSystem = "file time (unspecified)";
    auto all = readPos(input, o.format, timeSystem);
    if (!o.start.empty() && parseTime(o.start) > all.front().time + EPS)
        throw std::runtime_error("Acquisition start must be at or before the first solution epoch.");
    std::vector<Epoch> rows;
    for (const auto& e : all) if (e.quality == o.quality) rows.push_back(e);
    const double dt = sampleInterval(rows);
    if (o.windowSeconds < 2 * dt || o.referenceSeconds < 2 * dt)
        throw std::runtime_error("Reference and RMS intervals must cover at least two sampling intervals.");
    const double maxGap = o.maxGap > 0 ? o.maxGap : 1.5 * dt;
    if (maxGap < dt) throw std::runtime_error("--max-gap is smaller than median sampling interval.");
    const auto ref = referenceMean(rows, o.referenceSeconds, dt, maxGap);
    setErrors(all, ref);
    setErrors(rows, ref);
    const auto cv = convergence(rows, o, dt, maxGap);
    outputs(o, all, rows, ref, cv, dt, maxGap, timeSystem);
    return 0;
}
#ifdef _WIN32
int wmain(int argc, wchar_t* argv[]) {
    SetConsoleOutputCP(CP_UTF8);
    std::vector<std::string> args;
    for (int i = 0; i < argc; ++i) args.push_back(fs::path(argv[i]).u8string());
#else
int main(int argc, char* argv[]) {
    std::vector<std::string> args(argv, argv + argc);
#endif
    try { return run(args); }
    catch (const std::exception& e) { std::cerr << "ppp_statistics: " << e.what() << '\n'; return 1; }
}
