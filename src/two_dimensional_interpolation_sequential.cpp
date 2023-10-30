#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>
#include <optional>
#include <filesystem>
#include <sstream>

#include <CLI/CLI.hpp>

void noop(...)
{}

struct
{
    std::optional<std::filesystem::path> filename_option;
} config;

class BasePolynomial {
  public:
    BasePolynomial(int i, int j, const std::shared_ptr<std::vector<double>> &xValues, const std::shared_ptr<std::vector<double>> &yValues)
            : i(i), j(j), x_values(xValues), y_values(yValues)
    {}

    double evaluate(double x, double y) {
        double divider = 1, result = 1;
        for (int p = 0; p < x_values->size(); ++p) {
            for (int q = 0; q < y_values->size(); ++q) {
                if (p != i && q != j) {
                    result *= (x - x_values->at(i)) * (y - y_values->at(j));
                    divider *= (x_values->at(p) - x_values->at(i)) * (y_values->at(q) - y_values->at(j));
                }
            }
        }
        return result / divider;
    }

  private:
    int i, j;
    std::shared_ptr<std::vector<double>> x_values, y_values;
};

//std::function<double(double, double)> create_base_polynomial(const std::vector<double>& x_values, const std::vector<double>& y_values, int i, int j) {
//    return [i, j, x_values, y_values](double x, double y) -> double {
//        double divider = 1, result = 1;
//        for (int p = 0; p < x_values.size(); ++p) {
//            for (int q = 0; q < y_values.size(); ++q) {
//                if (p != i && q != j) {
//                    result *= (x - x_values[i]) * (y - y_values[j]);
//                    divider *= (x_values[p] - x_values[i]) * (y_values[q] - y_values[j]);
//                }
//            }
//        }
//        return result / divider;
//    };
//}

class LagrangePolynomial {
  public:
    LagrangePolynomial(const std::shared_ptr<std::vector<double>> &xValues, const std::shared_ptr<std::vector<double>> &yValues,
                       const std::shared_ptr<std::vector<std::vector<double>>> &zValues) : x_values(xValues),
                                                                                           y_values(yValues),
                                                                                           z_values(zValues),
                                                                                           base_polynomials(std::make_shared<std::vector<std::vector<BasePolynomial>>>())
    {
        for (int j = 0; j < y_values->size(); ++j) {
            base_polynomials->emplace_back();
            for (int i = 0; i < x_values->size(); ++i) {
                base_polynomials->back().emplace_back(i, j, x_values, y_values);
            }
        }
    }

    double evaluate(double x, double y) {
        double result = 0;
        for (int i = 0; i < (*x_values).size(); ++i) {
            for (int j = 0; j < (*y_values).size(); ++j) {
                result += (*z_values)[j][i] * (*base_polynomials)[j][i].evaluate(x, y);
            }
        }
        return result;
    }

  private:
    const std::shared_ptr<std::vector<double>> x_values;
    const std::shared_ptr<std::vector<double>> y_values;
    const std::shared_ptr<std::vector<std::vector<double>>> z_values;
    const std::shared_ptr<std::vector<std::vector<BasePolynomial>>> base_polynomials;
};

//std::function<double(double, double)> create_lagrange_polynomial(
//        const std::shared_ptr<std::vector<double>> x_values,
//        const std::shared_ptr<std::vector<double>> y_values,
//        const std::shared_ptr<std::vector<std::vector<double>>> z_values
//)
//{
//    std::shared_ptr<std::vector<std::vector<BasePolynomial>>> base_polynomials;
//
//    for (int j = 0; j < y_values->size(); ++j) {
//        base_polynomials->emplace_back();
//        for (int i = 0; i < x_values->size(); ++i) {
//            base_polynomials->back().emplace_back(i, j, x_values, y_values);
//        }
//    }
//
//    return [x_values, y_values, z_values, base_polynomials](double x, double y) -> double {
//        double result = 0;
//        for (int i = 0; i < x_values.size(); ++i) {
//            for (int j = 0; j < y_values.size(); ++j) {
//                result += z_values[j][i] * base_polynomials[j][i](x, y);
//            }
//        }
//        return result;
//    };
//}

int main(int argc, char **argv)
{
    CLI::App app;

    app.add_option("-f,--file", config.filename_option, "Read linear equations system from file (extended matrix format)")->check(
            CLI::ExistingFile
    );

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

    auto x_values = std::make_shared<std::vector<double>>();
    auto y_values = std::make_shared<std::vector<double>>();
    auto z_values = std::make_shared<std::vector<std::vector<double>>>();

    std::string line;
    while (line.empty()) {
        std::getline(*in_stream, line);
    }
    std::istringstream iss(line);
    double current;
    while (iss >> current) {
        x_values->push_back(current);
    }

    while (!in_stream->eof())
    {
        std::getline(*in_stream, line);
        if (line.empty()) continue;
        iss = std::istringstream(line);

        iss >> current;
        y_values->push_back(current);

        z_values->emplace_back();
        while (iss >> current)
        {
            z_values->back().push_back(current);
        }
    }

    LagrangePolynomial polynomial(x_values, y_values, z_values);
    auto [minX, maxX] = std::minmax_element(x_values->begin(), x_values->end());
    auto [minY, maxY] = std::minmax_element(y_values->begin(), y_values->end());

    std::vector<double> result_x_values, result_y_values;
    std::vector<std::vector<double>> result_z_values;

    double step = 0.2;
    for (double j = *minY; j <= *maxY; j += step) {
        result_y_values.push_back(j);
    }
    for (double i = *minX; i <= *maxX; i += step){
        result_x_values.push_back(i);
    }

    for (int j = 0; j < result_y_values.size(); ++j) {
        result_z_values.emplace_back();
        for (int i = 0; i < result_x_values.size(); ++i){
            result_z_values.back().push_back(polynomial.evaluate(i, j));
        }
    }

    std::cout << ", ";
    for (double result_x_value : result_x_values) {
        std::cout << result_x_value << ", ";
    }
    std::cout << '\n';

    for (int j = 0; j < result_y_values.size(); ++j) {
        std::cout << result_y_values[j] << ", ";
        for (int i = 0; i < result_x_values.size(); ++i) {
            std::cout << result_z_values[j][i] << ", ";
        }
        std::cout << '\n';
    }

    return 0;
}