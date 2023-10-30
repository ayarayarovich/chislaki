#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <string>

namespace matrix {
    using element_type = double;
    using row_type = std::vector<element_type>;
    using type = std::vector<row_type>;
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

std::pair<matrix::type, matrix::type> lu_decompose(const matrix::type& a) {
    const int n = static_cast<int>(a.size());

    matrix::type l(n, matrix::row_type(n, 0));
    matrix::type u(n, matrix::row_type(n, 0));

    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (j < i)
                l[j][i] = 0;
            else
            {
                l[j][i] = a[j][i];
                for (k = 0; k < i; k++)
                {
                    l[j][i] = l[j][i] - l[j][k] * u[k][i];
                }
            }
        }
        for (j = 0; j < n; j++)
        {
            if (j < i)
                u[i][j] = 0;
            else if (j == i)
                u[i][j] = 1;
            else
            {
                u[i][j] = a[i][j] / l[i][i];
                for (k = 0; k < i; k++)
                {
                    u[i][j] = u[i][j] - ((l[i][k] * u[k][j]) / l[i][i]);
                }
            }
        }
    }

    return {l, u};
}

int main() {
    std::cout << "\nMake sure that the input.txt file that I will be reading the matrix from is in this folder: " << std::filesystem::current_path().string() << "\n\n";

    std::ifstream ifs("input.txt");

    matrix::type A;
    std::vector<matrix::element_type> B;

    while (!ifs.eof()) {
        std::string line;
        std::getline(ifs, line);
        if (line.empty()) continue;
        std::istringstream iss(line);

        A.emplace_back();

        std::string token;
        while (iss >> token)
        {
            if (token == "|") {
                iss >> token;
                B.push_back(std::stod(token));
            }
            else {
                A.back().push_back(std::stod(token));
            }
        }
    }

    if (A.size() != A.front().size()) {
        std::cerr << "It is not A square matrix. Quit." << std::endl;
    }

    std::cout << "Matrix A:\n"
              << A << '\n';
    std::cout << "Vector B:\n";
    for (double bi : B) {
        std::cout << bi << '\n';
    }

    auto [L, U] = lu_decompose(A);

    std::cout << "\n"
              << "Matrix L:\n"
              << L << '\n'
              << "Matrix U:\n"
              << U << '\n';

    const int variables_count = static_cast<int>(A.size()); // since A is square matrix

    // Ly = B -> solve for vector y
    std::vector<matrix::element_type> y(variables_count);
    for (int i = 0; i < variables_count; ++i) {
        matrix::element_type rhs = B[i];
        for (int k = 0; k < i; ++k) {
            rhs -= y[k] * L[i][k];
        }
        y[i] = rhs / L[i][i];
    }

    // Ux = y -> solve for vector x
    std::vector<matrix::element_type> x(variables_count);
    for (int i = variables_count - 1; i >= 0; --i) {
        matrix::element_type rhs = y[i];
        for (int k = i + 1; k < variables_count; ++k) {
            rhs -= x[k] * U[i][k];
        }
        x[i] = rhs;
    }

    std::cout << "Solution: \n";
    for (int i = 0; i < variables_count; ++i) {
        std::cout << "x" << i + 1 << " = " << x[i] << '\n';
    }


    return 0;
}