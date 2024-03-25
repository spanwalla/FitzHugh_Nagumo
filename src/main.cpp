#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <unordered_map>

// du/dt = u - (u^3)/3 - v + I
double voltage(double u, const std::unordered_map<char, double>& params) {
    return u - pow(u, 3) / 3 - params.at('v') + params.at('I');
}

// T*(dv/dt) = u + a - bv
double recovery(double v, const std::unordered_map<char, double>& params) {
    return (params.at('u') + params.at('a') - params.at('b') * v) / params.at('T');
}

double rungeKutta(double(*function)(double, const std::unordered_map<char, double>&), double value,
                  const std::unordered_map<char, double>& params, const std::vector<std::vector<double>>& tableau) {
    std::vector<double> k(tableau[0].size());
    double weighted_sum = 0;
    for (int i = 0; i < tableau[0].size(); ++i) {
        double sum = 0;
        // для использования неявных методов можно будет заменить j < i на j < tableau[0].size()
        for (int j = 0; j < i; ++j) {
            sum += tableau[i][j] * k[j];
        }
        k[i] = function(value + params.at('h') * sum, params);
        weighted_sum += tableau.back()[i] * k[i];
    }
    return value + params.at('h') * weighted_sum;
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
    return tableau;
}

int main() {
    auto rk4_table = readButcherTableauFromFile("../data/tableau/rk4.txt");

    std::unordered_map<char, double> state = {{'a', 0.7}, {'b', 0.8}, {'I', 0.5}, {'T', 12.5},
                                              {'u', 0}, {'v', 0}, {'h', 0.1}};
    std::cout << "t u\n";
    double t = 0;
    while (t <= 3) {
        std::cout << t << "," << 100 * state.at('u') << '\n';
        state['u'] = rungeKutta(voltage, state.at('u'), state, rk4_table);
        state['v'] = rungeKutta(recovery, state.at('v'), state, rk4_table);
        t += state.at('h');
    }
    return 0;
}