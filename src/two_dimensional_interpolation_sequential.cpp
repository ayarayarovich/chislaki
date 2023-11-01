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

double base_polynomial(
        const std::vector<double> &x_values,
        const std::vector<double> &y_values,
        int i,
        int j,
        double x,
        double y
)
{
    double divider = 1, result = 1;
    for (int q = 0; q < y_values.size(); ++q)
    {
        for (int p = 0; p < x_values.size(); ++p)
        {
            if (p != i && q != j)
            {
                result *= (x - x_values[p]) * (y - y_values[q]);
                divider *= (x_values[i] - x_values[p]) * (y_values[j] - y_values[q]);
            }
        }
    }
    return result / divider;
}

class LagrangePolynomial
{
  public:
    LagrangePolynomial(
            const std::vector<double> &xValues,
            const std::vector<double> &yValues,
            const std::vector<std::vector<double>> &zValues
            )
            : x_values(xValues),
              y_values(yValues),
              z_values(zValues)
    {}

    double evaluate(double x, double y)
    {
        double result = 0;
        for (int i = 0; i < x_values.size(); ++i)
        {
            for (int j = 0; j < y_values.size(); ++j)
            {
                result += z_values[j][i] * base_polynomial(x_values, y_values, i, j, x, y);
            }
        }
        return result;
    }

  private:
    const std::vector<double> &x_values;
    const std::vector<double> &y_values;
    const std::vector<std::vector<double>> &z_values;
};

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

    std::vector<double> x_values;
    std::vector<double> y_values;
    std::vector<std::vector<double>> z_values;

    std::string line;
    while (line.empty())
    {
        std::getline(*in_stream, line);
    }
    std::istringstream iss(line);
    double current;
    while (iss >> current)
    {
        x_values.push_back(current);
    }

    while (!in_stream->eof())
    {
        std::getline(*in_stream, line);
        if (line.empty()) continue;
        iss = std::istringstream(line);

        iss >> current;
        y_values.push_back(current);

        z_values.emplace_back();
        while (iss >> current)
        {
            z_values.back().push_back(current);
        }
    }

    LagrangePolynomial polynomial(x_values, y_values, z_values);
    auto [minX, maxX] = std::minmax_element(x_values.begin(), x_values.end());
    auto [minY, maxY] = std::minmax_element(y_values.begin(), y_values.end());

    std::vector<double> result_x_values, result_y_values;
    std::vector<std::vector<double>> result_z_values;

    double step = 0.05;
    for (double j = *minY; j <= *maxY + 0.00001; j += step)
    {
        result_y_values.push_back(j);
    }
    for (double i = *minX; i <= *maxX + 0.00001; i += step)
    {
        result_x_values.push_back(i);
    }

    for (double result_y_value : result_y_values)
    {
        result_z_values.emplace_back();
        for (double result_x_value : result_x_values)
        {
            result_z_values.back().push_back(polynomial.evaluate(result_x_value, result_y_value));
        }
    }


    std::ofstream ofs("result.txt");
    std::stringstream sstream;

    sstream << "; ";
    for (double result_x_value: result_x_values)
    {
        sstream << result_x_value << "; ";
    }
    sstream << '\n';

    for (int j = 0; j < result_y_values.size(); ++j)
    {
        sstream << result_y_values[j] << "; ";
        for (int i = 0; i < result_x_values.size(); ++i)
        {
            sstream << result_z_values[j][i] << "; ";
        }
        sstream << '\n';
    }

    std::string result = sstream.str();
    std::replace(result.begin(), result.end(), '.', ',');

    ofs << result;

    return 0;
}