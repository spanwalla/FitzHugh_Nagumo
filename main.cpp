#include <iostream>
#include <cmath>
#include <unordered_map>

// du/dt = u - (u^3)/3 - v + I
double voltage(double u, const std::unordered_map<char, double>& values) {
    double v = values.at('v');
    return u - pow(u, 3) / 3 - v + values.at('I');
}

// T*(dv/dt) = u + a - bv
double recovery(double v, const std::unordered_map<char, double>& values) {
    return (values.at('u') + values.at('a') - values.at('b') * v) / values.at('T');
}

double rungeKutta(double(*function)(double, const std::unordered_map<char, double>&), double value, const std::unordered_map<char, double>& params) {
    double table[5][4] = {
            {0, 0, 0, 0},
            {0.5, 0, 0, 0},
            {0, 0.5, 0, 0},
            {0, 0, 1, 0},
            {0.1667, 0.3333, 0.3333, 0.1667}
    };

    double k[4];
    double weighted_sum = 0;
    for (int i = 0; i < 4; ++i) {
        double sum = 0;
        // для использования неявных методов можно будет заменить j < i на j < 4
        for (int j = 0; j < i; ++j) {
            sum += table[i][j] * k[j];
        }
        k[i] = function(value + params.at('h') * sum, params);
        weighted_sum += table[4][i] * k[i];
    }
    return value + params.at('h') * weighted_sum;
}

int main() {
    std::unordered_map<char, double> state = {{'a', 0.7}, {'b', 0.8}, {'I', 0.5}, {'T', 12.5},
                                              {'u', 0}, {'v', 0}, {'h', 0.1}};
    std::cout << voltage(state.at('u'), state) << ' ' << recovery(state.at('v'), state) << '\n';
    std::cout << "t u\n";
    double t = 0;
    while (t <= 50) {
        std::cout << t << "," << 100 * state.at('u') << '\n';
        state['u'] = rungeKutta(voltage, state.at('u'), state);
        state['v'] = rungeKutta(recovery, state.at('v'), state);
        t += state.at('h');
    }
    return 0;
}
