#include <bits/stdc++.h>

/* Use as a normal C++ dist.
 * Note that it returns a "sorted" zipfian distribution
 * You may want to hash its outputs
 */
class zipf_distribution {
    private:
        uint64_t n;
        long double theta, alpha, zetan, eta;

        long double zeta(uint64_t n, long double theta) {
            long double out=0;
            for(uint64_t i=1; i<=n; ++i) {
                long double powr = pow((long double)1/(long double)i, theta);
                out += powr;
                if(powr<1e-6) break;
            }
            return out;
        }

    public:
        /* n_ is the maximum number
         * theta_ is the zipf distribution parameter (decay)
         */
        zipf_distribution(uint64_t n_, long double theta_): n(n_), theta(theta_) {
            alpha = 1/(1-theta);
            zetan = zeta(n, theta);
            eta = (1-pow((long double)2/n, (long double)1-theta))/(1-zeta(2, theta) / zetan);
        }

        template< class Generator >
            uint64_t operator()( Generator& g ) {
                uniform_real_distribution<> dist(0.0, 1.0);
                long double u = dist(g);
                long double uz = u*zetan;
                if (uz < 1) return 1;
                if (uz < 1 + pow(0.5, theta)) return 2;
                return 1 + (uint64_t)(n * pow(eta*u - eta + 1, alpha));
            }
};
