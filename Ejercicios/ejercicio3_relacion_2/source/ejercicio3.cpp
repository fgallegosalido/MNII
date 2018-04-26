#include <iostream>
#include <utility>
#include <vector>
#include <array>
#include <cmath>

#define UNICODE_SUPPORT
#include "polynomial.hpp"

using punto = std::pair<double, double>;

auto formula_interpolacion_3_puntos(std::array<punto, 3> &&p){
    double c0 = p[0].second / ((p[0].first-p[1].first)*(p[0].first-p[2].first));
    double c1 = p[1].second / ((p[1].first-p[0].first)*(p[1].first-p[2].first));
    double c2 = p[2].second / ((p[2].first-p[0].first)*(p[2].first-p[1].first));

    detail::polynomial_double ret{p[1].first*p[2].first, (-p[1].first-p[2].first), 1}; ret *= c0;
    ret += c1*detail::polynomial_double{p[0].first*p[2].first, (-p[0].first-p[2].first), 1};
    ret += c2*detail::polynomial_double{p[0].first*p[1].first, (-p[0].first-p[1].first), 1};

    return ret;
}

auto formula_derivacion_3_puntos(std::array<punto, 3> &&p){
    return formula_interpolacion_3_puntos(std::forward<std::array<punto, 3>>(p)).differentiate();
}

int main(){
    std::array<punto, 4> v{{{2.9, -4.827866},
                            {3.0, -4.240058},
                            {3.1, -3.496909},
                            {3.2, -2.596792}}};

    auto p1 = formula_interpolacion_3_puntos(std::array<punto, 3>{v[0], v[1], v[2]});
    auto p2 = formula_interpolacion_3_puntos(std::array<punto, 3>{v[0], v[1], v[3]});
    auto p3 = formula_interpolacion_3_puntos(std::array<punto, 3>{v[0], v[2], v[3]});
    auto p4 = formula_interpolacion_3_puntos(std::array<punto, 3>{v[1], v[2], v[3]});

    std::cout << "Las fórmulas interpolatorias son:\n"
                << "   " << p1 << " sin (3.2,-2.596792)" << std::endl
                << "   " << p2 << " sin (3.1,-3.496909)" << std::endl
                << "   " << p3 << " sin (3.0,-4.240058)" << std::endl
                << "   " << p4 << " sin (2.9,-4.827866)" << std::endl;

    std::cout.precision(12);
    std::cout << "\nErrores con respecto al punto que falta:\n"
                << "   (3.2,-2.596792): " << std::abs(v[3].second - p1(v[3].first)) << std::endl
                << "   (3.1,-3.496909): " << std::abs(v[2].second - p2(v[2].first)) << std::endl
                << "   (3.0,-4.240058): " << std::abs(v[1].second - p3(v[1].first)) << std::endl
                << "   (2.9,-4.827866): " << std::abs(v[0].second - p4(v[0].first)) << std::endl;

    p3.differentiate();

    std::cout << "\nEl mejor polinomio parece ser el tercero. Su derivada es:\n"
                <<  "   " << p3 << std::endl
                << "\nLa tabla queda:\n\n"
                << "    x  |   f(x)    |    f'(x)\n";

    for(const auto &i : v){
        std::cout << "  ---------------------------------\n";
        std::cout << "   " << i.first;

        if(i.first-std::floor(i.first) == 0){
            std::cout << ".0";
        }

        std::cout << " | " << i.second << " | " << p3(i.first) << std::endl;
    }
}
