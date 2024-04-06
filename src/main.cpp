#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <algorithm>
#include <filesystem>

/* TODO:
 * Разбить на файлы
 */

struct Method {
    std::string path;
    std::pair<int, int> order;
};

std::string findFileInDirectory(const std::string& filename) {
    for (std::filesystem::recursive_directory_iterator i("data"), end; i != end; ++i)
        if (!is_directory(i->path()) && i->path().filename() == filename)
            return i->path().string();
    return "";
}

std::vector<std::string> split(const std::string& source, char delimiter) {
    std::vector<std::string> result;
    std::stringstream ss(source);
    std::string word;
    while (!ss.eof()) {
        getline(ss, word, delimiter);
        result.push_back(word);
    }
    return result;
}

std::vector<std::vector<double>> readButcherTableauFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Unable to read tableau.");

    std::vector<std::vector<double>> tableau;
    std::string row_string;
    double element = 0;
    while (getline(file, row_string)) {
        std::stringstream ss(row_string);
        std::vector<double> row;
        while (ss >> element)
            row.push_back(element);
        tableau.push_back(row);
    }
    file.close();

    if (tableau[0].size() >= tableau.size())
        throw std::invalid_argument("Incorrect tableau.");

    return tableau;
}

template <class T>
void write2DVectorToFile(const std::string& filename, const std::vector<std::vector<T>>& sequence,
                      const std::string& head, const std::string& delimiter = ",") {
    std::ofstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Unable to open file for write.");

    file << head << '\n';
    for (size_t i = 0; const auto& row : sequence) {
        if (i++ != 0)
            file << '\n';

        for (size_t j = 0; const auto& element : row) {
            if (j++ != 0)
                file << delimiter;
            file << element;
        }
    }

    file.close();
}

std::unordered_map<std::string, Method> getAvailableMethods() {
    std::string methods_path = findFileInDirectory("methods.txt");
    if (methods_path.empty())
        throw std::runtime_error("File 'methods.txt' not found in data folder.");

    std::ifstream file(methods_path);
    if (!file.is_open())
        throw std::runtime_error("Unable to open file 'methods.txt'.");

    std::string row;
    std::string element;
    std::unordered_map<std::string, Method> methods;
    while (getline(file, row)) {
        auto data = split(row, ':');
        auto order_str = split(data[1], '/');
        std::pair<int, int> order {stoi(order_str[0]), order_str.size() > 1 ? stoi(order_str[1]) : 0};
        methods[data[0]] = {data[2], order};
    }
    file.close();
    return methods;
}

std::unordered_map<std::string, double> readStateFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Unable to open file for read.");

    std::unordered_map<std::string, double> state;
    std::string row;
    std::string key;
    while (getline(file, row)) {
        std::stringstream ss(row);
        ss >> key >> state[key];
    }
    file.close();
    return state;
}

// du/dt = u - (u^3)/3 - v + I
double voltage(double u, const std::unordered_map<std::string, double>& params) {
    return u - pow(u, 3) / 3 - params.at("v") + params.at("I");
}

// T*(dv/dt) = u + a - bv
double recovery(double v, const std::unordered_map<std::string, double>& params) {
    return (params.at("u") + params.at("a") - params.at("b") * v) / params.at("T");
}

std::vector<double> rungeKutta(double(*function)(double, const std::unordered_map<std::string, double>&),
                               double value, const std::unordered_map<std::string, double>& params,
                               const std::vector<std::vector<double>>& tableau) {
    std::vector<double> k(tableau.back().size());
    std::vector<double> weighted_sum(tableau.size() - tableau.back().size(), 0);

    for (int i = 0; i < tableau.back().size(); ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += tableau[i][j] * k[j];
        }
        k[i] = function(value + params.at("h") * sum, params);

        for (int j = 0; j < weighted_sum.size(); ++j) {
            weighted_sum[j] += tableau[tableau.back().size() + j][i] * k[i];
        }
    }
    std::vector<double> answer(weighted_sum.size());
    std::transform(weighted_sum.cbegin(), weighted_sum.cend(), answer.begin(),
                   [&](double weight) -> double { return value + params.at("h") * weight; });
    return answer;
}

std::vector<std::vector<double>> fitzHughNagumoModel(double duration, const std::string& method,
                                                     std::unordered_map<std::string, double> state) {
    if (duration <= 0)
        throw std::invalid_argument("Duration must be more than zero.");

    auto available_methods = getAvailableMethods();
    auto tableau = readButcherTableauFromFile(available_methods[method].path);

    const auto tolerance = [](std::pair<double, double> values, double absolute, double relative){
        return absolute + relative * std::max(std::abs(values.first), std::abs(values.second));
    };

    std::vector<std::vector<double>> points;
    double t = 0;
    while (t <= duration) {
        points.push_back({t, state.at("u"), state.at("v")});

        auto u_values = rungeKutta(voltage, state.at("u"), state, tableau);
        auto v_values = rungeKutta(recovery, state.at("v"), state, tableau);

        // Попробуем подобрать шаг
        if (u_values.size() > 1 && v_values.size() > 1) {
            // err = sqrt(1/2*((epsilon_u/tol_u)^2+(epsilon_v/tol_v)^2)) - standard deviation
            double err = sqrt((pow((u_values[1] - u_values[0]) /
                                   tolerance({u_values[0], u_values[1]}, state.at("Atol_u"),
                                             state.at("Rtol_u")), 2) + pow((v_values[1] - v_values[0]) /
                                                                           tolerance({v_values[0], v_values[1]},
                                                                                     state.at("Atol_v"),
                                                                                     state.at("Rtol_v")), 2)) / 2);
            // h_opt = h * (1/err)^(1/(min(p, p^)+1))
            double h_opt = state.at("h") * pow(1 / err, 1.0 / (std::min(available_methods[method].order.first,
                                                                 available_methods[method].order.second) + 1));
            state["h"] = h_opt;
            u_values = rungeKutta(voltage, state.at("u"), state, tableau);
            v_values = rungeKutta(recovery, state.at("v"), state, tableau);
        }

        state["u"] = u_values[0];
        state["v"] = v_values[0];
        t += state.at("h");
    }
    return points;
}

int main(int argc, char** argv) {
    auto state = readStateFromFile("data/input/a.txt");
    auto points = fitzHughNagumoModel(120, "dp87", state);
    write2DVectorToFile("data/result/dp.txt", points, "t,u,v");
    return 0;
}