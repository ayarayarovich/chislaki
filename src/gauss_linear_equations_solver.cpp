#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>
#include <optional>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <unordered_set>
#include <unordered_map>

#include <CLI/CLI.hpp>

void noop(...)
{}

#define ALGORITHMS(algo) \
    algo(COLUMN_SEARCH)  \
    algo(ROW_SEARCH)     \
    algo(FULL_SEARCH)

#define _CREATE_ALGORITHMS_ENUM(name) ALG_##name
#define CREATE_ALGORITHMS_ENUM(name) _CREATE_ALGORITHMS_ENUM(name),
#define CREATE_STRINGS(name) #name,
#define CREATE_STRING_TO_ENUM_MAP(name) {#name, Algorithms::_CREATE_ALGORITHMS_ENUM(name)},

enum class Algorithms
{
    ALGORITHMS(CREATE_ALGORITHMS_ENUM)
};
std::unordered_set<std::string> AVAILABLE_ALGORITHMS{
    ALGORITHMS(CREATE_STRINGS)
};
std::unordered_map<std::string, Algorithms> MAP_ALGORITHMS_STR_TO_ENUM{
    ALGORITHMS(CREATE_STRING_TO_ENUM_MAP)
};

struct
{
    std::optional<std::filesystem::path> filename_option;
    Algorithms algorithm = Algorithms::ALG_FULL_SEARCH;
} config;

namespace matrix
{
    using element_type = double;
    using row_type = std::vector<element_type>;
    using type = std::vector<row_type>;

    void swap_rows(type &mtx, int y1, int y2)
    {
        if (y1 == y2) return;
        std::swap(mtx[y1], mtx[y2]);
    }

    void swap_columns(type &mtx, int x1, int x2)
    {
        if (x1 == x2) return;
        for (row_type &y: mtx)
        {
            std::swap(y[x1], y[x2]);
        }
    }
}

std::ostream &operator <<(std::ostream &os, const matrix::type &mtx)
{
    size_t min_field_width = 0;
    std::stringstream ss;
    for (const matrix::row_type &row: mtx)
    {
        for (const matrix::element_type &value: row)
        {
            ss << std::setprecision(3) << std::fixed << value;
            min_field_width = std::max(min_field_width, ss.str().size());
            ss.str("");
        }
    }

    for (const matrix::row_type &row: mtx)
    {
        for (const matrix::element_type &value: row)
        {
            os << std::setw(min_field_width) << std::setprecision(3) << std::fixed << value << ' ';
        }
        os << std::endl;
    }

    return os;
}

matrix::type gauss_per_column_algo(matrix::type mtx)
{
    const int MATRIX_SIZE = static_cast<int>(mtx.size());

    // Direct pass
    for (int k = 0; k < MATRIX_SIZE; ++k)
    {
        int mxy = 0;
        matrix::element_type mx = 0;
        for (int yi = k; yi < MATRIX_SIZE; ++yi)
        {
            const matrix::element_type abs_val = abs(mtx[yi][k]);
            if (abs_val > mx)
            {
                mx = abs_val;
                mxy = yi;
            }
        }

        matrix::swap_rows(mtx, k, mxy);

        for (int yi = k + 1; yi < MATRIX_SIZE; ++yi)
        {
            const matrix::element_type scale = -mtx[yi][k] / mtx[k][k];
            for (int xi = k; xi < MATRIX_SIZE + 1; ++xi)
            {
                mtx[yi][xi] += mtx[k][xi] * scale;
            }
        }
    }

    // Inverse pass
    for (int k = MATRIX_SIZE - 1; k >= 0; --k)
    {
        for (int yi = k - 1; yi >= 0; --yi)
        {
            const matrix::element_type scale = -mtx[yi][k] / mtx[k][k];
            for (int xi = k; xi < MATRIX_SIZE + 1; ++xi)
            {
                mtx[yi][xi] += mtx[k][xi] * scale;
            }
        }
    }
    return mtx;
}

matrix::type gauss_per_row_algo(matrix::type mtx)
{
    const int MATRIX_SIZE = static_cast<int>(mtx.size());

    // Used to keep track of column permutations during horizontal search for max element
    std::vector<int> columns(MATRIX_SIZE);
    std::iota(columns.begin(), columns.end(), 0);

    for (int k = 0; k < MATRIX_SIZE; ++k)
    {
        // Direct pass
        int mxx = 0;
        matrix::element_type mx = 0;
        for (int xi = k; xi < MATRIX_SIZE; ++xi)
        {
            const matrix::element_type abs_val = abs(mtx[k][xi]);
            if (abs_val > mx)
            {
                mx = abs_val;
                mxx = xi;
            }
        }

        matrix::swap_columns(mtx, k, mxx);
        std::swap(columns[k], columns[mxx]);

        for (int yi = k + 1; yi < MATRIX_SIZE; ++yi)
        {
            const matrix::element_type scale = -mtx[yi][k] / mtx[k][k];
            for (int xi = k; xi < MATRIX_SIZE + 1; ++xi)
            {
                mtx[yi][xi] += mtx[k][xi] * scale;
            }
        }
    }

    for (int k = MATRIX_SIZE - 1; k >= 0; --k)
    {
        // Inverse pass
        for (int yi = k - 1; yi >= 0; --yi)
        {
            const matrix::element_type scale = -mtx[yi][k] / mtx[k][k];
            for (int xi = k; xi < MATRIX_SIZE + 1; ++xi)
            {
                mtx[yi][xi] += mtx[k][xi] * scale;
            }
        }
    }

    // Placing rearranged columns to their right place
    for (int i = 0; i < MATRIX_SIZE;)
    {
        const int desired_position = columns[i];
        const int actual_position = i;

        if (desired_position == actual_position)
        {
            ++i;
        }
        else
        {
            matrix::swap_rows(mtx, desired_position, actual_position);
            std::swap(columns[desired_position], columns[actual_position]);
        }
    }

    return mtx;
}

matrix::type gauss_full_search_algo(matrix::type mtx)
{
    const int MATRIX_SIZE = static_cast<int>(mtx.size());

    // Used to keep track of column permutations during horizontal search for max element
    std::vector<int> columns(MATRIX_SIZE);
    std::iota(columns.begin(), columns.end(), 0);

    for (int k = 0; k < MATRIX_SIZE; ++k)
    {
        // Direct pass
        int mxy = 0;
        int mxx = 0;
        matrix::element_type mx = 0;

        // Search for max element vertically
        for (int yi = k; yi < MATRIX_SIZE; ++yi)
        {
            const matrix::element_type abs_val = abs(mtx[yi][k]);
            if (abs_val > mx)
            {
                mx = abs_val;
                mxy = yi;
                mxx = k;
            }
        }

        // Search for max element horizontally
        for (int xi = k; xi < MATRIX_SIZE; ++xi)
        {
            const matrix::element_type abs_val = abs(mtx[k][xi]);
            if (abs_val > mx)
            {
                mx = abs_val;
                mxy = k;
                mxx = xi;
            }
        }

        matrix::swap_rows(mtx, k, mxy);
        matrix::swap_columns(mtx, k, mxx);
        std::swap(columns[k], columns[mxx]);

        for (int yi = k + 1; yi < MATRIX_SIZE; ++yi)
        {
            const matrix::element_type scale = -mtx[yi][k] / mtx[k][k];
            for (int xi = k; xi < MATRIX_SIZE + 1; ++xi)
            {
                mtx[yi][xi] += mtx[k][xi] * scale;
            }
        }
    }

    for (int k = MATRIX_SIZE - 1; k >= 0; --k)
    {
        // Inverse pass
        for (int yi = k - 1; yi >= 0; --yi)
        {
            const matrix::element_type scale = -mtx[yi][k] / mtx[k][k];
            for (int xi = k; xi < MATRIX_SIZE + 1; ++xi)
            {
                mtx[yi][xi] += mtx[k][xi] * scale;
            }
        }
    }

    // Placing rearranged columns to their right place
    for (int i = 0; i < MATRIX_SIZE;)
    {
        const int desired_position = columns[i];
        const int actual_position = i;

        if (desired_position == actual_position)
        {
            ++i;
        }
        else
        {
            matrix::swap_rows(mtx, desired_position, actual_position);
            std::swap(columns[desired_position], columns[actual_position]);
        }
    }

    return mtx;
}

// Just a placeholder instead of actual variable. So CLI11 doesn't assign anything to an actual config.algorithm field,
// and we can control assignment/validation.
// Since there is some #define magic which generates enum with available algorithms and creates MAP_ALGORITHMS_STR_TO_ENUM std::unordered_map,
// so we can map string representation of algorithm to its enum representation.
std::string placeholder;

int main(int argc, char **argv)
{
    CLI::App app;

    app.add_option("-f,--file", config.filename_option, "Read linear equations system from file (extended matrix format)")->check(
            CLI::ExistingFile);
    app.add_option("-a,--algorithm", placeholder, "Algorithm to use for Gaussian elimination")->transform([](const std::string &algorithm) {
        if (AVAILABLE_ALGORITHMS.contains(algorithm))
        {
            config.algorithm = MAP_ALGORITHMS_STR_TO_ENUM[algorithm];
            return "";
        }
        else
        {
            return "Unknown algorithm specified :(";
        }
    });

    try
    {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError &e)
    {
        return app.exit(e);
    }

    std::shared_ptr<std::istream> in_stream;

    if (config.filename_option)
    {
        std::cerr << "Filename provided: " << config.filename_option.value() << " - reading equations from there." << std::endl;
        in_stream = std::make_shared<std::ifstream>(config.filename_option.value());
    }
    else
    {
        std::cerr << "Filename is missing. Using stdin." << std::endl;
        in_stream.reset(&std::cin, noop);
    }

    matrix::type mtx;
    std::optional<size_t> first_row_size;

    while (!in_stream->eof())
    {
        std::string line;
        std::getline(*in_stream, line);
        if (line.empty()) continue;
        std::istringstream iss(line);

        mtx.emplace_back();

        double number;
        while (iss >> number)
        {
            mtx.back().push_back(number);
        }

        if (first_row_size)
        {
            if (mtx.back().size() != first_row_size.value())
            {
                std::cerr << "Invalid matrix: different sizes of rows. Exit." << std::endl;
                std::cerr << mtx << std::endl;
                std::cerr << mtx.back().size() << ' ' << first_row_size.value() << std::endl;
                return 1;
            }
        }
        else
        {
            first_row_size = mtx.back().size();
        }
    }

    if (mtx.size() + 1 != mtx.begin()->size())
    {
        std::cerr << "Invalid matrix: height + 1 != width. Exit.";
        return 1;
    }

//    TODO: добавить проверку на совместность системы, а пока ненадо...

    matrix::type result;

    switch (config.algorithm)
    {
        case Algorithms::ALG_COLUMN_SEARCH:
        {
            std::cerr << "Using column search algorithm..." << std::endl;
            result = gauss_per_column_algo(mtx);
            break;
        }
        case Algorithms::ALG_ROW_SEARCH:
        {
            std::cerr << "Using row search algorithm..." << std::endl;
            result = gauss_per_row_algo(mtx);
            break;
        }
        case Algorithms::ALG_FULL_SEARCH:
        {
            std::cerr << "Using full search algorithm..." << std::endl;
            result = gauss_per_row_algo(mtx);
            break;
        }
    }

    std::cout << result << std::endl;

    for (int i = 0; i < result.size(); ++i) {
        result[i][result.size()] /= result[2 - i][i];
        std::cout << "X" << i << ": " << result[i][result.size()] << std::endl;
    }

    std::cout << result << std::endl;

    return 0;
}