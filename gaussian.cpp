#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

template <typename T>
void require_finite(T value) {
    if (!std::isfinite(value)) {
        throw std::runtime_error("Arithmetic overflow or non-finite value encountered.");
    }
}

template <typename T>
void solve(const std::filesystem::path& input_path, bool spp) {
    std::ifstream input(input_path);
    if (!input) throw std::runtime_error("Cannot open input file.");

    // Parse the dimension separately so fractional or malformed sizes are rejected.
    std::string size_token;
    if (!(input >> size_token) || size_token.empty() ||
        size_token.find_first_not_of("0123456789") != std::string::npos) {
        throw std::runtime_error("The first value must be a positive integer n.");
    }
    const auto dimension = std::stoull(size_token);
    if (dimension == 0 || dimension > std::numeric_limits<std::size_t>::max()) {
        throw std::runtime_error("Invalid system dimension.");
    }
    const auto n = static_cast<std::size_t>(dimension);
    std::vector<std::vector<T>> a(n, std::vector<T>(n));
    std::vector<T> b(n), scale(n), x(n);

    auto read_value = [&input]() {
        T value;
        if (!(input >> value)) throw std::runtime_error("Missing or invalid numeric input.");
        require_finite(value);
        return value;
    };
    for (auto& row : a) {
        for (auto& value : row) value = read_value();
    }
    for (auto& value : b) value = read_value();
    std::string extra;
    if (input >> extra) throw std::runtime_error("Unexpected data after the right-hand side.");

    // Scales are computed from the original coefficient rows, excluding b.
    if (spp) {
        for (std::size_t i = 0; i < n; ++i) {
            for (T value : a[i]) scale[i] = std::max(scale[i], std::abs(value));
            if (scale[i] == T(0)) {
                throw std::runtime_error("A zero coefficient row prevents a unique solution.");
            }
        }
    }

    for (std::size_t k = 0; k < n; ++k) {
        if (spp) {
            std::size_t pivot = k;
            T best_ratio = std::abs(a[k][k]) / scale[k];
            require_finite(best_ratio);
            for (std::size_t i = k + 1; i < n; ++i) {
                const T ratio = std::abs(a[i][k]) / scale[i];
                require_finite(ratio);
                if (ratio > best_ratio) {
                    best_ratio = ratio;
                    pivot = i;
                }
            }
            std::swap(a[k], a[pivot]);
            std::swap(b[k], b[pivot]);
            std::swap(scale[k], scale[pivot]);
        }
        if (a[k][k] == T(0)) {
            throw std::runtime_error(spp
                ? "Zero pivot: no unique solution can be computed in this precision."
                : "Zero pivot in naive elimination. Try --spp; the system may be singular.");
        }
        for (std::size_t i = k + 1; i < n; ++i) {
            const T factor = a[i][k] / a[k][k];
            require_finite(factor);
            a[i][k] = T(0);
            for (std::size_t j = k + 1; j < n; ++j) {
                a[i][j] -= factor * a[k][j];
                require_finite(a[i][j]);
            }
            b[i] -= factor * b[k];
            require_finite(b[i]);
        }
    }

    // Back substitution on the upper triangular system.
    for (std::size_t i = n; i-- > 0;) {
        T sum = b[i];
        for (std::size_t j = i + 1; j < n; ++j) {
            sum -= a[i][j] * x[j];
            require_finite(sum);
        }
        x[i] = sum / a[i][i];
        require_finite(x[i]);
    }

    auto output_path = input_path;
    output_path.replace_extension(".sol");
    std::ofstream output(output_path);
    if (!output) throw std::runtime_error("Cannot open solution file for writing.");
    output << std::setprecision(std::numeric_limits<T>::max_digits10);
    for (std::size_t i = 0; i < n; ++i) {
        if (i != 0) output << ' ';
        output << x[i];
    }
    output << '\n';
    output.close();
    if (!output) throw std::runtime_error("Failed to write solution file.");
}

int main(int argc, char* argv[]) {
    bool spp = false, use_double = false;
    std::filesystem::path input_path;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--spp") spp = true;
        else if (arg == "--double") use_double = true;
        else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: gaussian [--spp] [--double] system.lin\n"
                         "Default: naive Gaussian elimination, single precision.\n";
            return 0;
        } else if (arg.empty() || arg[0] == '-' || !input_path.empty()) {
            std::cerr << "Error: unknown option or extra argument: " << arg << '\n';
            return 1;
        } else input_path = arg;
    }
    if (input_path.empty() || input_path.extension() != ".lin") {
        std::cerr << "Usage: gaussian [--spp] [--double] system.lin\n"
                     "Input must have the .lin extension.\n";
        return 1;
    }
    try {
        if (use_double) solve<double>(input_path, spp);
        else solve<float>(input_path, spp);
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
