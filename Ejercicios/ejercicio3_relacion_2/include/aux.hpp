#pragma once

#include <string> // std::string, std::to_string
#include <limits> // std::numeric_limits

namespace detail{
    namespace aux{
        namespace unicode{
            /* Array that contains numbers from 0 to 9 (including empty string)
             *
             * If the user defines the macro UNICODE_SUPPORT, an array of
             * superscript/subscript numbers will be generated instead of normal ones
             */
        #ifdef UNICODE_SUPPORT
            constexpr const char* superscript_chars[] = {
            "\u2070", "\u00B9", "\u00B2", "\u00B3", "\u2074", "\u2075", "\u2076", "\u2077", "\u2078", "\u2079", ""};
            constexpr const char* subscript_chars[] = {
            "\u2080", "\u2081", "\u2082", "\u2083", "\u2084", "\u2085", "\u2086", "\u2087", "\u2088", "\u2089", ""};
            constexpr const char* integer_set_symbol = "\u2124";
        #else
            constexpr const char* superscript_chars[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ""};
            constexpr const char* subscript_chars[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ""};
            constexpr const char* integer_set_symbol = "Z";
        #endif

            /* Functions that returns a number from 0 to 9 in superscript/subscript
            * (if enabled by the user through the UNICODE_SUPPORT macro)
            */
            constexpr auto superscript (std::size_t n){
                return (n<=9) ? superscript_chars[n]
                        : superscript_chars[10];   // Any other number will return an empty string
            }
            constexpr auto subscript (std::size_t n){
                return (n<=9) ? subscript_chars[n]
                        : subscript_chars[10];   // Any other number will return an empty string
            }
        }  // namespace unicode

        /* Function that returns any unsigned number in superscript (exponents)
         *
         * If UNICODE_SUPPORT is not enabled, it will return the string "^n"
         */
        std::string exponent_string (std::size_t  n){
            std::string exponent = "";

        #ifdef UNICODE_SUPPORT
            if (n==1 || n==0) return exponent;

            while (n>0){
                exponent = unicode::superscript(n%10) + exponent;
                n /= 10;
            }
        #else
            exponent = "^" + std::to_string(n);
        #endif

            return exponent;
        }

        /* Function that returns any unsigned number in subscript (subindex)
         *
         * If UNICODE_SUPPORT is not enabled, it will return the string "_n"
         */
        std::string subindex_string (std::size_t  n){
            std::string subindex = "";

        #ifdef UNICODE_SUPPORT
            if (n==0) return std::string(unicode::subscript(0));

            while (n>0){
                subindex = unicode::subscript(n%10) + subindex;
                n /= 10;
            }
        #else
            subindex = "_" + std::to_string(n);
        #endif

            return subindex;
        }
    }
}  // namespace detail::aux

namespace detail{
    namespace aux{
        /* Primality test for compile-time evaluation
         *
         * Source code taken from Casey's answer at Stack Overflow:
         * https://stackoverflow.com/questions/18303632/compile-time-prime-checking
         */
        namespace primes{

            template<typename U>
            constexpr U mid(U low, U high) {
                return (low + high) / 2;
            }

            // precondition: low*low <= n, high*high > n.
            template<typename U>
            constexpr U ceilsqrt (U n, U low, U high){
                return low + 1 >= high
                    ? high
                    : (mid(low, high) * mid(low, high) == n)
                        ? mid(low, high)
                        : (mid(low, high) * mid(low, high) <  n)
                            ? ceilsqrt(n, mid(low, high), high)
                            : ceilsqrt(n, low, mid(low, high));
            }

            // returns ceiling(sqrt(n))
            template<typename U>
            constexpr U ceilsqrt (U n){
                return n < 3
                    ? n
                        : ceilsqrt(n, U(1), U(1) << (std::numeric_limits<U>::digits / 2));
            }


            // returns true if n is divisible by an odd integer in
            // [2 * low + 1, 2 * high + 1).
            template<typename U>
            constexpr bool find_factor (U n, U low, U high){
                return low + 1 >= high
                    ? (n % (2 * low + 1)) == 0
                        : (find_factor(n, low, mid(low, high))
                            || find_factor(n, mid(low, high), high));
            }

        }  // namespace primes

        template<typename U>
        constexpr bool is_prime (U n){
            return n > 1
                && (n == 2
                    || (n % 2 == 1
                        && (n == 3
                            || !primes::find_factor(n, U(1), (primes::ceilsqrt(n) + 1) / 2))));
        }
    }
}  // namespace detail::aux
