#pragma once

#include "aux.hpp"

#include <iostream>  // std::cout, std::cin, std::endl
#include <vector> // std::vector
#include <deque> // std::deque
#include <algorithm> // std::transform
#include <initializer_list>   // std::initializer_list
#include <iterator>  // std::iterator_traits
#include <utility>   // std::move
#include <complex>   // std::complex
#include <cmath>  // std::pow()
#include <cstdint>   // std::int8_t, std::int16_t, std::int32_t, std::int64_t

// Library for Z-module rings and fields
#ifdef Z_MODULE_SUPPORT
    #include "z_module.hpp" // detail::ZModule
    #include "z_module_prime.hpp" // detail::ZModulePrime
#endif

// Boost libraries for some interesting number types
#ifdef RATIONAL_SUPPORT
    #include <boost/rational.hpp> // boost::rational
#endif

#if defined(MULTIPRECISION_SUPPORT)
    #include <boost/multiprecision/cpp_int.hpp>  // cpp_int, cpp_rational
    #include <boost/multiprecision/cpp_bin_float.hpp>  // cpp_bin_float

    typedef boost::multiprecision::cpp_int multiprecision_int;
    typedef boost::multiprecision::cpp_rational multiprecision_rational;
    typedef boost::multiprecision::cpp_bin_float<100> multiprecision_float;
#elif defined(FAST_MULTIPRECISION_SUPPORT)
    #include <boost/multiprecision/gmp.hpp>   // mpz_int, mpq_rational, mpf_float_100

    typedef boost::multiprecision::mpz_int multiprecision_int;
    typedef boost::multiprecision::mpq_rational multiprecision_rational;
    typedef boost::multiprecision::mpf_float_100 multiprecision_float;
#endif

namespace detail{

    // Alias to get the category of an iterator (random access, bidirectional,...)
    template <typename Iterator>
    using traits = typename std::iterator_traits<Iterator>::iterator_category;

    // Class Polynomial. CType is the type of the coefficients
    template <typename CType, typename Container = std::deque<CType>>
    class Polynomial{

        private:

            Container coeffs;  // Actual coefficients of the polynomial
            char var = 'x';   // Letter that identifies the variable

            // Helper function to adjust the degree, so the last coefficient is not 0
            void adjust_degree (){
                while (coeffs.back()==static_cast<value_type>(0) && coeffs.size()>1){
                    coeffs.pop_back();
                }
            }

            // Helper constructor for input iterators acceptance
            template <typename InputIterator>
            Polynomial (InputIterator first, InputIterator last, std::input_iterator_tag)
                : coeffs(first, last) {adjust_degree();}

        public:

            // typedefs for the value_type (coefficients) and the size_type
            typedef typename Container::value_type value_type;
            typedef typename Container::size_type size_type;

            // typedefs for iterator validity
            typedef typename Container::iterator iterator;
            typedef typename Container::const_iterator const_iterator;
            typedef typename Container::reverse_iterator reverse_iterator;
            typedef typename Container::const_reverse_iterator const_reverse_iterator;

            // Some constructors
            explicit Polynomial (const value_type& x = value_type(0))
                : coeffs(1, x) {}
            Polynomial (std::initializer_list<value_type> l)
                : coeffs(l) {adjust_degree();}

            // Range constructor using tag dispatching for input iterators
            template <typename InputIterator>
            Polynomial (InputIterator first, InputIterator last)
                : Polynomial(first, last, traits<InputIterator>()) {}

            template <typename Cont>
            Polynomial (const Cont& cont)
                : Polynomial(cont.begin(), cont.end()) {}

            value_type get_coefficient (size_type i) const{
                return (i<coeffs.size()) ? coeffs[i] : static_cast<value_type>(0);
            }

            void set_coefficient (size_type i, value_type elem){
                if (i>=coeffs.size()){
                    coeffs.resize(i+1, 0);
                }

                coeffs[i] = elem;
                adjust_degree();
            }

            value_type operator[] (size_type i) const{
                return get_coefficient(i);
            }

            char get_variable () const{
                return var;
            }

            void set_variable (const char& c) {
                var = c;
            }

            size_type degree () const{
                return coeffs.size()-1;
            }

            /* Evaluates the polynomial for the value x using the Horner's
             * polynomial evaluation scheme.
             *
             * RType is the type of the evaluation
             */
            template<typename RType>
            RType evaluate_at (const RType& x) const{
                RType res = static_cast<RType>(coeffs.back());
                for (size_type i=coeffs.size()-1; i>0; --i){
                    res = static_cast<RType>(coeffs[i-1]) + res*x;
                }

            return res;
            }

            /* Enable natural evaluation of a mathematical function, so you
             * can write p(x) instead of p.evaluate_at(x). RType is the type
             * of the evaluation
             */
            template<typename RType>
            RType operator() (const RType& x) const{
                return evaluate_at(x);
            }

            // Operator overloadings for polynomials arithmetic
            Polynomial& operator+= (const Polynomial& pol){
                size_type size = pol.coeffs.size();
                if (size > coeffs.size()) coeffs.resize(size, 0);

                for (size_type i=0; i<size; ++i){
                    coeffs[i] += pol.coeffs[i];
                }

                if (size == coeffs.size()) adjust_degree();

                return *this;
            }

            Polynomial& operator-= (const Polynomial& pol){
                size_type size = pol.coeffs.size();
                if (size > coeffs.size()) coeffs.resize(size, 0);

                for (size_type i=0; i<size; ++i){
                    coeffs[i] -= pol.coeffs[i];
                }

                if (size == coeffs.size()) adjust_degree();

                return *this;
            }

            Polynomial& operator*= (const Polynomial& pol){
                coeffs.resize(coeffs.size()+pol.coeffs.size()-1, 0);

                for (int i=coeffs.size()-pol.coeffs.size(); i>=0; --i){
                    for (size_type j=pol.coeffs.size()-1; j>0; --j){
                        coeffs[i+j] += coeffs[i]*pol.coeffs[j];
                    }
                    coeffs[i] = coeffs[i]*pol.coeffs[0];
                }

                return *this;
            }

            Polynomial& operator/= (const Polynomial& pol){
                if (coeffs.size() < pol.coeffs.size()){
                    coeffs.resize(1);
                    coeffs[0] = 0;
                    return *this;
                }

                Container coc(coeffs.size()-pol.coeffs.size()+1, 0);

                for (int i=coc.size()-1; i>=0; --i){
                    coc[i] = coeffs[pol.coeffs.size()+i-1]/pol.coeffs.back();
                    for (int j=pol.coeffs.size()-2; j>=0; --j){
                        coeffs[i+j] -= pol.coeffs[j]*coc[i];
                    }
                }

                coeffs = std::move(coc);
                return *this;
            }

            Polynomial& operator%= (const Polynomial& pol){
                if (coeffs.size() < pol.coeffs.size()) return *this;

                value_type coc;

                for (int i=coeffs.size()-pol.coeffs.size(); i>=0; --i){
                    coc = coeffs[pol.coeffs.size()+i-1]/pol.coeffs.back();
                    for (int j=pol.coeffs.size()-1; j>=0; --j){
                        coeffs[i+j] -= pol.coeffs[j]*coc;
                    }
                }

                adjust_degree();
                return *this;
            }

            template <typename U>
            Polynomial& operator+=(const U& other){
                coeffs[0] += static_cast<value_type>(other);
                return *this;
            }
            template <typename U>
            Polynomial& operator-=(const U& other){
                coeffs[0] -= static_cast<value_type>(other);
                return *this;
            }
            template <typename U>
            Polynomial& operator*=(const U& other){
                std::transform(coeffs.begin(), coeffs.end(), coeffs.begin(),
                                [&other](const value_type& val){
                                    return val*static_cast<value_type>(other);
                                });
                return *this;
            }
            template <typename U>
            Polynomial& operator/=(const U& other){
                std::transform(coeffs.begin(), coeffs.end(), coeffs.begin(),
                                [&other](const value_type& val){
                                    return val/static_cast<value_type>(other);
                                });
                return *this;
            }
            template <typename U>
            Polynomial& operator%=(const U& other){
                coeffs.resize(1, 0);
                return *this;
            }

            bool operator== (const Polynomial& pol) const{
                if (coeffs.size() != pol.coeffs.size()) return false;

                for (size_type i=0; i<coeffs.size(); ++i){
                    if (coeffs[i] != pol.coeffs[i]) return false;
                }

                return true;
            }

            bool operator!= (const Polynomial& pol) const{
                return !(*this == pol);
            }

            /* Modifies the coefficients so they match the polynomial to the
             * power of n (n must be unsigned type)
             */
            template <typename UInt>
            Polynomial pow(const Polynomial& pol, const std::make_unsigned_t<UInt>& n) {
                if (n==0){
                    coeffs[0] = 1;
                    coeffs.resize(1);
                    return *this;
                }

                for (unsigned i=1; i<n; ++i){
                    (*this) *= (*this);
                }

                return *this;
            }

            /* Modifies the coefficients so they match the derivative of
             * the polynomial defined by *this
             */
            Polynomial& differentiate (){
                if (coeffs.size() == 1){
                    coeffs[0] = 0;
                }
                else{
                    for (size_type i=1; i<coeffs.size(); ++i){
                        coeffs[i-1] = coeffs[i]*static_cast<value_type>(i);
                    }
                    coeffs.resize(coeffs.size()-1);
                }

                return *this;
            }

            /* Modifies the coefficients so they match the antiderivative, with
             * integration constant equals to c
             */
            Polynomial& integrate_const (const value_type& c = 0){
                if (coeffs.size()==1 && coeffs[0]==0){
                    coeffs[0] = c;
                    return *this;
                }

                coeffs.resize(coeffs.size()+1);

                for (size_type i=coeffs.size()-1; i>0; --i){
                    coeffs[i] = coeffs[i-1]/static_cast<value_type>(i);
                }
                coeffs[0] = c;

                return *this;
            }

            /* Modifies the coefficients so they match the antiderivative that
             * meets the condition of F(x)=y (where F is the antiderivative).
             * RType is the type of the evaluation
             */
            template<typename RType>
            Polynomial& integrate_point (const RType& x, const RType& y){
                coeffs[0] = static_cast<value_type>(y - (*this).integrate_const().evaluate_at(x));
                return *this;
            }

            /* Calculates the definite integral between two values of the
             * polynomial. RType is the type of the evaluation
             */
            template <typename RType>
            RType definite_integral (const RType& lower_bound, const RType& upper_bound){
                Polynomial p(*this);
                p.integrate_const();
                return p.evaluate_at(upper_bound) - p.evaluate_at(lower_bound);
            }

            // Swap two polynomials
            friend void swap (Polynomial& lhs, Polynomial& rhs){
                using std::swap;

                swap(lhs.coeffs, rhs.coeffs);
                swap(lhs.var, rhs.var);
            }

            /* TODO: Decide a good way to input a polynomial (or let the user
             * implement it itself).
             *
             * One idea could be to implement a parser, so it would work for std::cin
             * and a new constructor.
             */
            // friend std::istream& operator>>(std::istream& is, Polynomial& pol);

            /* Pretty print for the polynomial
             *
             * If UNICODE_SUPPORT is enabled, the polynomial will be printed
             * with superscript characters instead of the expresion "^n"
             */
            friend std::ostream& operator<< (std::ostream& os, const Polynomial& pol){
                if (pol.degree() == 0){
                    os << pol[0];
                }
                else{
                    if (pol[pol.degree()] == -1){
                        os << "-";
                    }
                    else if (pol[pol.degree()] != 1){
                        os << pol[pol.degree()];
                    }
                    os << pol.get_variable() << aux::exponent_string(pol.degree());

                    for (size_type i=pol.degree()-1; i>0; --i){
                        if (pol[i] != 0){
                            if (pol[i] > 0){
                                os << "+";

                                if (pol[i] != 1){
                                    os << pol[i];
                                }
                            }
                            else if (pol[i] < 0){
                                if (pol[i] == -1){
                                    os << "-";
                                }
                                else{
                                    os << pol[i];
                                }
                            }

                            os << pol.get_variable() << aux::exponent_string(i);
                        }
                    }

                    if (pol[0] != 0){
                        if (pol[0] > 0){
                            os << "+";
                        }
                        os << pol[0];
                    }
                }

                return os;
            }

            /* Iterator functions to iterate through a polynomial
             *
             * It's based on the vector's iterator, so this is just a wrapper
             * for polynomials
             */
            iterator begin(){
                return coeffs.begin();
            }
            const_iterator begin() const{
                return coeffs.begin();
            }

            iterator end(){
                return coeffs.end();
            }
            const_iterator end() const{
                return coeffs.end();
            }

            reverse_iterator rbegin(){
                return coeffs.rbegin();
            }
            const_reverse_iterator rbegin() const{
                return coeffs.rbegin();
            }

            reverse_iterator rend(){
                return coeffs.rend();
            }
            const_reverse_iterator rend() const{
                return coeffs.rend();
            }

            const_iterator cbegin(){
                return coeffs.cbegin();
            }

            const_iterator cend(){
                return coeffs.cend();
            }

            const_reverse_iterator crbegin(){
                return coeffs.crbegin();
            }

            const_reverse_iterator crend(){
                return coeffs.crend();
            }

    };

    template <typename CType>
    Polynomial<CType> operator+(Polynomial<CType> rhs){
        return rhs;
    }
    template <typename CType>
    Polynomial<CType> operator-(Polynomial<CType> rhs){
        return rhs *= -1;
    }

    template <typename CType>
    Polynomial<CType> operator+(Polynomial<CType> lhs, const Polynomial<CType>& rhs){
        return lhs += rhs;
    }
    template <typename CType, typename U>
    Polynomial<CType> operator+(Polynomial<CType> lhs, const U& rhs){
        return lhs += rhs;
    }
    template <typename CType, typename U>
    Polynomial<CType> operator+(const U& lhs, Polynomial<CType> rhs){
        return rhs += lhs;
    }

    template <typename CType>
    Polynomial<CType> operator-(Polynomial<CType> lhs, const Polynomial<CType>& rhs){
        return lhs -= rhs;
    }
    template <typename CType, typename U>
    Polynomial<CType> operator-(Polynomial<CType> lhs, const U& rhs){
        return lhs -= rhs;
    }
    template <typename CType, typename U>
    Polynomial<CType> operator-(const U& lhs, Polynomial<CType> rhs){
        return (rhs -= lhs) *= -1;
    }

    template <typename CType>
    Polynomial<CType> operator*(Polynomial<CType> lhs, const Polynomial<CType>& rhs){
        return lhs *= rhs;
    }
    template <typename CType, typename U>
    Polynomial<CType> operator*(Polynomial<CType> lhs, const U& rhs){
        return lhs *= rhs;
    }
    template <typename CType, typename U>
    Polynomial<CType> operator*(const U& lhs, Polynomial<CType> rhs){
        return rhs *= lhs;
    }

    template <typename CType>
    Polynomial<CType> operator/(Polynomial<CType> lhs, const Polynomial<CType>& rhs){
        return lhs /= rhs;
    }
    template <typename CType, typename U>
    Polynomial<CType> operator/(Polynomial<CType> lhs, const U& rhs){
        return lhs /= rhs;
    }
    template <typename CType, typename U>
    Polynomial<CType> operator/(const U& lhs, Polynomial<CType> rhs){
        return Polynomial<CType>(lhs) /= rhs;
    }

    template <typename CType>
    Polynomial<CType> operator%(Polynomial<CType> lhs, const Polynomial<CType>& rhs){
        return lhs %= rhs;
    }
    template <typename CType, typename U>
    Polynomial<CType> operator%(Polynomial<CType> lhs, const U& rhs){
        return lhs %= rhs;
    }
    template <typename CType, typename U>
    Polynomial<CType> operator%(const U& lhs, Polynomial<CType> rhs){
        return Polynomial<CType>(lhs) %= rhs;
    }

   /* Like Polynomial::pow, but returning an rvalue
    * (not modifying the original)
    */
   template <typename CType>
   const Polynomial<CType> pow(const Polynomial<CType>& pol, const unsigned n) {
      return Polynomial<CType>(pol).pow(n);
   }

   /* Like Polynomial::differentiate, but returning an rvalue
    * (not modifying the original)
    */
   template <typename CType>
   const Polynomial<CType> differentiate (const Polynomial<CType>& pol){
      return Polynomial<CType>(pol).differentiate();
   }

   /* Like Polynomial::integrate_const, but returning an rvalue
    * (not modifying the original)
    */
   template <typename CType>
   const Polynomial<CType> integrate_const (const Polynomial<CType>& pol, const CType& c = 0){
      return Polynomial<CType>(pol).integrate_const(c);
   }

   /* Like Polynomial::integrate_point, but returning an rvalue
    * (not modifying the original)
    */
   template <typename CType, typename RType>
   const Polynomial<CType> integrate_point (const Polynomial<CType>& pol, const RType& x, const RType& y){
      return Polynomial<CType>(pol).integrate_point(x, y);
   }

   // GCD of two polynomials using Euclidean's algorithm
   template <typename CType>
   const Polynomial<CType> gcd (Polynomial<CType> lhs, Polynomial<CType> rhs){
      if (lhs==Polynomial<CType>()) return rhs;
      if (rhs==Polynomial<CType>()) return lhs;

      while (rhs.degree()>0 || rhs[0]!=0){
         swap(lhs, rhs);
         rhs %= lhs;
      }

      return lhs;
   }

   // LCM of two polynomials using Euclidean's algorithm
   template <typename CType>
   const Polynomial<CType> lcm (const Polynomial<CType>& lhs, const Polynomial<CType>& rhs){
      return (lhs/gcd(lhs, rhs))*rhs;
   }

   // Typedefs for coefficients with integer types
   typedef Polynomial<std::int8_t> polynomial_int8;
   typedef Polynomial<std::int16_t> polynomial_int16;
   typedef Polynomial<std::int32_t> polynomial_int32;
   typedef Polynomial<std::int64_t> polynomial_int64;

   // Typedefs for coefficients in floating point (real numbers)
   typedef Polynomial<float> polynomial_float;
   typedef Polynomial<double> polynomial_double;
   typedef Polynomial<long double> polynomial_long_double;

   // Typedefs for extra types (boost::rational and std::complex)
#ifdef RATIONAL_SUPPORT
   typedef Polynomial<boost::rational<std::int8_t>> polynomial_rational_int8;
   typedef Polynomial<boost::rational<std::int16_t>> polynomial_rational_int16;
   typedef Polynomial<boost::rational<std::int32_t>> polynomial_rational_int32;
   typedef Polynomial<boost::rational<std::int64_t>> polynomial_rational_int64;
#endif

   typedef Polynomial<std::complex<float>> polynomial_complex_float;
   typedef Polynomial<std::complex<double>> polynomial_complex_double;
   typedef Polynomial<std::complex<long double>> polynomial_complex_long_double;

   /* Typedefs for multiprecision polynomial using boost multiprecision
    * library or fast multiprecision using gmp library
    */
#if defined(MULTIPRECISION_SUPPORT) || defined(FAST_MULTIPRECISION_SUPPORT)
   typedef Polynomial<multiprecision_int> polynomial_multiprecision_int;
   typedef Polynomial<multiprecision_rational> polynomial_multiprecision_rational;
   typedef Polynomial<multiprecision_float> polynomial_multiprecision_float;
   typedef Polynomial<std::complex<multiprecision_float>> polynomial_multiprecision_complex;
#endif

#ifdef Z_MODULE_SUPPORT
   template<auto N>
   using polynomial_modular = Polynomial<ZModule<N>>;
   template<auto N>
   using polynomial_modular_prime = Polynomial<ZModulePrime<N>>;
#endif

}  // namespace detail
