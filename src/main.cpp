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
 * Подбор шага
 * Разбить на файлы
 * Реализовать FSAL так чтобы они make sense
 */

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

std::unordered_map<std::string, std::string> getAvailableMethods() {
    std::string methods_path = findFileInDirectory("methods.txt");
    if (methods_path.empty())
        throw std::runtime_error("File 'methods.txt' not found in data folder.");

    std::ifstream file(methods_path);
    if (!file.is_open())
        throw std::runtime_error("Unable to open file 'methods.txt'.");

    std::string row;
    std::string element;
    std::unordered_map<std::string, std::string> methods;
    while (getline(file, row)) {
        auto data = split(row, ':');
        methods[data[0]] = data[1];
    }
    file.close();
    return methods;
}

std::unordered_map<char, double> readStateFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Unable to open file for read.");

    std::unordered_map<char, double> state;
    std::string row;
    char key;
    while (getline(file, row)) {
        std::stringstream ss(row);
        ss >> key >> state[key];
    }
    file.close();
    return state;
}

// du/dt = u - (u^3)/3 - v + I
double voltage(double u, const std::unordered_map<char, double>& params) {
    return u - pow(u, 3) / 3 - params.at('v') + params.at('I');
}

// T*(dv/dt) = u + a - bv
double recovery(double v, const std::unordered_map<char, double>& params) {
    return (params.at('u') + params.at('a') - params.at('b') * v) / params.at('T');
}

std::vector<double> rungeKutta(double(*function)(double, const std::unordered_map<char, double>&),
                               double value, const std::unordered_map<char, double>& params,
                               const std::vector<std::vector<double>>& tableau) {
    std::vector<double> k(tableau.back().size());
    std::vector<double> weighted_sum(tableau.size() - tableau.back().size(), 0);

    for (int i = 0; i < tableau.back().size(); ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += tableau[i][j] * k[j];
        }
        k[i] = function(value + params.at('h') * sum, params);

        for (int j = 0; j < weighted_sum.size(); ++j) {
            weighted_sum[j] += tableau[tableau.back().size() + j][i] * k[i];
        }
    }
    std::vector<double> answer(weighted_sum.size());
    std::transform(weighted_sum.cbegin(), weighted_sum.cend(), answer.begin(), [&](double weight) -> double { return value + params.at('h') * weight; });
    return answer;
}

std::vector<std::vector<double>> fitzHughNagumoModel(double duration, const std::string& method,
                                                     std::unordered_map<char, double> state, double accuracy = 0.001) {
    if (duration <= 0)
        throw std::invalid_argument("Duration must be more than zero.");

    auto available_methods = getAvailableMethods();
    auto tableau = readButcherTableauFromFile(available_methods[method]);

    std::vector<std::vector<double>> points;
    double t = 0;
    while (t <= duration) {
        points.push_back({t, state.at('u'), state.at('v')});
        t += state.at('h');
        state['u'] = rungeKutta(voltage, state.at('u'), state, tableau)[0];
        state['v'] = rungeKutta(recovery, state.at('v'), state, tableau)[0];
    }
    return points;
}

int main(int argc, char** argv) {
    auto state = readStateFromFile("data/input/a.txt");
    auto points = fitzHughNagumoModel(160, "dp8", state);
    write2DVectorToFile("data/result/dp.txt", points, "t,u,v");
    return 0;
}