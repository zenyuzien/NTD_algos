
// the algorithm to find the solution is not optimal because the program focuses in describing the process which takes extra resources including time
// un-comment the function you wish to use, modify the param in the main function only. Don't disturb the code outside main function

#include <iostream>
#include <stdio.h>
#include <unordered_map>
#include <vector>
#include <cmath>
using namespace std;

int GCD(int,int); // greatest common divisor for 2 integers
int MI_EEA(int,int); // Multiplicative inverse of first parameter under modulus of 2nd parameter
int mod_exp(int,int,int); // modulo exponentiation method for finding mod, first param - base, 2nd param - exponent, 3rd param - mod 
bool MILLER_ROBIN(int); // MIller robin test to check primality of given integer 
int totient(int); // euler's totient function for given integer ϕ()
void fillPrimes();
bool is_prime(int);
int PRIMES[1000];
bool euler_th(int , int );
int CRT(vector<int>& , vector<int>& );


int main()
{
    fillPrimes(); // init function 
    
//    GCD(1701,3768); // gcd of give n 2 numbers
   
//    MI_EEA(3,5); // to find the Multiplicative inverse of 3 under mod 5 
   
//   MILLER_ROBIN(137); // check if 137 is prime
   
//   mod_exp(11,19,23); // 11^19 mod 23
   
//   totient(589); // finds ϕ(589)
   
   //euler_th(2,10);
   
   //mod_exp(3,4,10);
   
    
    
    //X ≡ 2 MOD 3 
    //X ≡ 3 MOD 5
    //X ≡ 2 MOD 7
    /*
    vector<int> a{2,3,2}; 
    vector<int> m{3,5,7};
    CRT(a,m); // chinese reminder theorum to solve above equations
    */
    return 0;
}


bool euler_th(int a, int n) {
    cout << "Euler's theorum: For every positive integer 'a' and 'n' which are said to be in relatively prime, then a^ϕ(n) ≡ 1 mod n \n";
    cout << "Firstly, checking if " << a << ", " << n << " are relatively primes:\n";
    if (GCD(a, n) == 1) cout << "Since they are relatively prime, checking furthur condition\n";
    else {
        cout << "They aren't relatively prime, hence euler theorum fails\nProof:\n";
    }
    cout << "finding totient function of " << n << ": ϕ(" << n << ")\n";
    int tot_n = totient(n);
    cout << "Now that we have ϕ(" << n << ") = " << tot_n << ", we check a^ϕ(n) mod n \n";
    cout << a << "^" << tot_n << " mod " << n << endl;
    cout << "using modular exponentiation to find the above value\n";
    int val = mod_exp(a, tot_n, n);
    if (val == 1) {
        cout << "Hence prooved Euler theory for " << a << " and " << n << endl;
        return true;
    } else cout << "Euler theory failed for " << a << " and " << n << endl;
    return false;
}
int totient(int n) {
    if (n == 1) {
        cout << "there are no positive integers less than 1\n";
        return 0;
    } else if (is_prime(n)) {
        cout << n << " is a prime number itself, hence we use the formula \n";
        cout << "ϕ(n) = (n-1)\n";
        cout << "=> ϕ(" << n << ") = " << n - 1;
        return n - 1;
    } else if (n < 10) {
        cout << "since there are only " << n - 1 << " numbers less than " << n << ", we can do brute force by counting relatively prime numbers among them\n";
        int r = 0;
        for (int i = 1; i < n; i++) {
            cout << "checking if " << i << " and " << n << " are relatively prime: \n";
            if (GCD(i, n) == 1) {
                r++;
                cout << "=> " << i << " and " << n << " are relatively prime\n";
            } else cout << "=> " << i << " and " << n << " aren't relatively prime\n";
        }
        cout << r;
        return r;
    }
    int N = n;
    vector < int > v;
    unordered_map < int, int > mp;
    cout << "factorising " << n << endl;
    if (n % 2 == 0) v.push_back(2);
    while (n % 2 == 0) {
        mp[2]++;
        cout << n << "/2 = ";
        n /= 2;
        cout << n << endl;
    }
    int ptr = 1;
    while (PRIMES[ptr] < N) {
        if (n % PRIMES[ptr] == 0) v.push_back(PRIMES[ptr]);
        while (n % PRIMES[ptr] == 0) {
            mp[PRIMES[ptr]]++;
            cout << n << "/" << PRIMES[ptr] << " = ";
            n /= PRIMES[ptr];
            cout << n << endl;
        }
        ptr++;
    }
    cout << "reult: " << N << " = ";
    for (auto i: v) cout << i << "^" << mp[i] << " . ";
    cout << endl;
    if (v.size() == 2 && (mp[v[0]] == 1 && mp[v[1]] == 1)) {
        cout << N << " can be written as product of 2 primes: " << v[0] << " and " << v[1] << endl;
        cout << "hence we use the formula: ϕ(p.q) = (p-1)*(q-1) \n";
        cout << "ϕ(" << N << ") = " << v[0] - 1 << " * " << v[1] - 1 << " = " << (v[0] - 1) * (v[1] - 1) << endl;
        return (v[0] - 1) * (v[1] - 1);
    }
    cout << "Hence we use the formula ϕ(n) = n*(1- 1/p1)*(1- 1/p2)..." << endl;
    cout << "the distinct primes are: ";
    for (auto i: v) cout << i << " ";
    cout << endl << "Now applying the above formula \n";
    cout << "ϕ(" << N << ") = " << N << "* ";
    for (auto i: v) cout << "(1- (1/" << i << "))";
    cout << endl << "ϕ(" << N << ") = " << N << "* ";
    for (auto i: v) cout << "(( " << i << " -1 )/" << i << ")";
    cout << endl << "ϕ(" << N << ") = " << N << "* ";
    double nr = 1, dr = 1;
    for (auto i: v) cout << "( " << i - 1 << "/ " << i << ")";
    cout << endl << "ϕ(" << N << ") = " << N << "* ";
    for (auto i: v) {
        nr *= (i - 1);
        dr *= i;
    }
    nr = (N * nr) / dr;
    cout << nr << "/" << dr << " = " << nr << endl;
    cout << "Hence ϕ(" << N << ") = " << nr << endl;
    return (N * nr) / dr;
}
bool is_prime(int n) {
    if (n <= 1) return false;
    int i = 2;
    while (i * i <= n) {
        if (n % i == 0) return false;
        i++;
    }
    return true;
}
void fillPrimes() {
    int i = 0, p = 3;
    PRIMES[i++] = 2;
    while (i < 1000) {
        if (is_prime(p)) PRIMES[i++] = p;
        p += 2;
    }
}
int GCD(int A, int B) {
    int a = A, b = B;
    cout << "GCD(" << a << ", " << b << ")  \n";
    if (b > a) {
        int tmp = a;
        a = b;
        b = tmp;
        cout << "= GCD(" << a << ", " << b << ") {GCD is commutative function}\n";
    }
    while (b) {
        cout << "= GCD(" << b << ", " << a << " mod " << b << ")\n";
        a = a % b;
        int tmp = a;
        a = b;
        b = tmp;
        cout << "= GCD(" << a << ", " << b << ")\n";
    }
    cout << "the GCD of " << A << " ," << B << " is " << a << endl;
    return a;
}
int MI_EEA(int a, int m) {
    int Q = m / a, A = m, B = a, R = A % B, T1 = 0, T2 = 1, T = -Q;
    B = B % m;
    if (B == 0) {
        cout << "Multiplicative inverse of 0 doesn't exist\n";
        return -1;
    }
    printf("\n      Q      |       A      |       B      |       R      |      T1     |      T2    |       T     |\n\n");
    while (B) {
        printf(" %6d      |  %6d      |  %6d      |  %6d      |  %6d     |  %6d    |  %6d     |\n", Q, A, B, R, T1, T2, T);
        cout << "T = T1 - (T2*Q) \n";
        cout << "T = " << T1 << " - ( " << T1 << " * " << Q << " )\n";
        cout << "T = " << T << endl << endl;
        A = B;
        B = R;
        if (B) {
            Q = A / B;
            R = A % B;
        }
        T1 = T2;
        T2 = T;
        T = T1 - (T2 * Q);
    }
    printf("      X      |  %6d      |  %6d      |      X       |  %6d     |  %6d    |      X      |\n", A, B, T1, T2);
    cout << "We got B=0 which ends the algorithm\nHence, the MI of " << a << " under mod " << m << " is " << T1 << endl;
    return T1;
}
int mod_exp(int a, int P, int m) {
    int r = 1;
    cout << a << "^" << P << " mod " << m << endl << endl;
    int pr = 2, tmp;
    cout << a << " mod " << m << " ≡ " << a << " mod " << m << endl;
    if (a >= m) {
        a = a % m;
        if (a == 0) {
            cout << "0 under any mod";
            return -1;
        }
        cout << a << " mod " << m << " ≡ " << a;
    }
    cout << endl;
    if (abs(a - m) < a) {
        cout << a << " mod " << m << " ≡ ";
        a = (a - m);
        cout << a << " mod " << m << endl;
    }
    unordered_map < int, int > mp;
    int mp_ptr = 1;
    mp[mp_ptr++] = a;
    cout << a << "^" << pr << " mod " << m << " ≡ ";
    int prev = a * a;
    cout << prev << " mod " << m;
    if (prev >= m) {
        prev = prev % m;
        cout << " ≡ " << prev << " mod " << m;
    }
    if (abs(prev - m) < prev) {
        prev = (prev - m);
        cout << " ≡ " << prev << " mod " << m << endl;
    }
    mp[mp_ptr++] = prev;
    cout << endl;
    pr *= 2;
    while (pr <= P) {
        if (a == 1 || a == -1) {
            cout << a << "^" << pr << " mod " << m << " ≡ " ;
            cout << " 1 " << "mod " << m << endl;
            mp[mp_ptr++] = 1;
            pr *= 2;
            continue;
        }
        cout << a << "^" << pr << " mod " << m << " ≡ ";
        printf("(%d^%d)^2 ", a, pr / 2);
        cout << " mod " << m << endl;
        cout << "substituing ";
        printf("%d^%d mod %d", a, pr / 2, m);
        cout << " from above, which is ";
        cout << prev << " mod " << m << endl;
        cout << a << "^" << pr << " mod " << m << " ≡ ";
        cout << prev << "*" << prev << " mod " << m;
        prev = prev * prev;
        cout << " ≡ " << prev << " mod " << m;
        if (prev >= m) {
            prev = prev % m;
            cout << " ≡ " << prev << " mod " << m;
        }
        cout << endl;
        if (abs(prev - m) < prev) {
            prev = (prev - m);
            cout << " ≡ " << prev << " mod " << m << endl;
        }
        mp[mp_ptr++] = prev;
        pr *= 2;
    }
    pr = P;
    cout << a << "^" << P << " mod " << m << " ≡ ";
    int mi = 1, mj = 1;
    while (pr) {
        if (pr % 2 == 1) cout << a << "^" << mi << "  ";
        pr = pr / 2;
        mi *= 2;
    }
    cout << "mod " << m << endl << "substituting corresponding mod values\n";
    pr = P;
    mi = 1;
    vector < int > v;
    bool inc;
    while (pr) {
        inc = false;
        if (pr % 2 == 1) {
            cout << mp[mj];
            v.push_back(mp[mj]);
            inc = true;
        }
        pr = pr / 2;
        if (pr && inc) cout << " . ";
        else if (pr == 0 && v.size()) cout << " mod " << m << endl;
        mj++;
    }
    int vs = v.size();
    while (vs >= 2) {
        v[vs - 2] *= v[vs - 1];
        vs--;
        for (auto i = 0; i < vs; i++) cout << v[i] << " . ";
        cout << "mod " << m << endl;
        if (v[vs - 1] >= m) {
            v[vs - 1] = v[vs - 1] % m;
            for (auto i = 0; i < vs; i++) cout << v[i] << " . ";
            cout << "mod " << m << endl;
        }
        if (v[vs - 1] < 0) {
            while (v[vs - 1] < 0) v[vs - 1] += m;
            for (auto i = 0; i < vs; i++) cout << v[i] << " . ";
            cout << "mod " << m << endl;
        }
    }
    if (v[0] < 0) {
        v[0] += m;
        cout << v[0] << " mod " << m << endl;
    }
    return v[0];
}
bool MILLER_ROBIN(int P) {
    cout << "in this solution, to compute modulus, modular exponentiation is used everytime for full description, if calculator is permitted, it can be used\n";
    int p = P;
    if (p <= 2) {
        cout << "MILLER ROBIN not applicable, choose a positive integer greater than 2 \n";
        return false;
    }
    cout << "MILLER ROBIN TEST FOR " << p << endl;
    cout << "Step 1: find k,m in\nn-1 = 2^k . m\n";
    p--;
    int k = 1;
    cout << "trying k= " << k << ": " << p << "/" << "2^" << k << " = ";
    float r = (float) p / ((float) pow(2, k));
    cout << r << endl;
    while (ceilf(r) == floorf(r)) {
        cout << "trying k= " << (++k) << ": " << p << "/" << "2^" << k << " = ";
        r = (float) p / ((float) pow(2, k));
        cout << r << endl;
    }
    k--;
    int pow_ = pow(2, k);
    cout << "hence k= " << k << ". Substituting in above eqn.." <<
        endl;
    cout << p << " = " << "2^" << k << " . " << "m" << endl;
    cout << p << " = " << pow_ << " . " << "m" << endl;
    cout << "m = " << p << "/" << pow_ << endl;
    int m = p / pow_;
    cout << "hence m = " << m << " and k = " << k << endl;
    cout << "step 2: choosing a from (1, " << p << "): choosing a = 2 {1<2<" << p << "}\n";
    int a = 2;
    cout << "step 3: compute b0 = a^m (mod n)\n";
    cout << "b0 = 2^" << m << " mod " << p + 1 << endl;
    cout << "finding 2^" << m << " mod " << p + 1 << endl;
    int b0 = mod_exp(2, m, p + 1);
    if (b0 - (p + 1) == -1) b0 = -1;
    cout << "b0 = " << b0 << endl;
    int bptr = 0;
    while (1) {
        if (abs(b0) == 1) {
            if (b0 == 1) {
                cout << "Hence, since we got b" << bptr << " as 1 which means it is composite\n";
                return false;
            } else {
                cout << "Hence, since we got b" << bptr << " as -1 which means it is prolly prime\n";
                return true;
            }
        } else {
            cout << "since b" << bptr << " isn't +1 or -1, we calculate b" << ++bptr << endl;
            cout << "b" << bptr << " = b" << bptr - 1 << "^2 mod " << p + 1 << endl;
            b0 = mod_exp(b0, 2, p + 1);
            if (b0 - (p + 1) == -1) b0 = -1;
        }
    }
}
int CRT(vector < int > & a, vector < int > & m) {
    if (a.size() == 0 || m.size() == 0 || (a.size() != m.size())) {
        cout << "incorrect params\n";
        return -1;
    }
    cout << "CRT: all above equations have a unique solution if the moduli are relatively prime\n";
    cout << "Firstly, given equations:\n";
    for (auto i = 0; i < a.size(); i++) cout << "X ≡ " << a[i] << " mod " << m[i] << endl;
    cout << "(this step is optional) - checking if the moduli are relatively prime, first.. \n";
    vector < int > tmp(m.begin(), m.end());
    int i = a.size() - 1;
    while (i >= 1) {
        int tp = GCD(tmp[i], tmp[i - 1]);
        if (tp != 1) {
            cout << "Given moduli aren't relatively prime, there exists no unique solution\n";
            return -1;
        }
        tmp[(i--) - 1] = tp;
    }
    cout << "Now confirm that moduli are relatively prime, proceeding:\n";
    cout << "We have \n";
    for (auto j = 0; j < a.size(); j++) cout << "a" << j + 1 << " = " << a[j] << ", m" << j + 1 << " = " << m[j] << endl;
    cout << "To find: M,\n";
    for (auto j = 0; j < a.size(); j++) cout << "M" << j + 1 << ", M" << j + 1 << "'\n";
    cout << "X = (a1M1M1' + a2M2M2' +..) mod M  where\n";
    cout << "M = product of all moduli\n";
    cout << "Mi = M/ mi \n";
    cout << "Mi' = Multiplicative inverse of Mi under mod mi\n";
    cout << "Calculating M: product of ";
    int M = 1;
    for (auto j: m) {
        M *= j;
        cout << j << " ";
    }
    cout << "= " << M << endl;
    vector < int > Mi, Minv;
    Mi.resize(m.size());
    Minv.resize(m.size());
    cout << "Mi = M/ mi \n";
    for (auto j = 0; j < m.size(); j++) {
        cout << "M" << j + 1 << " = " << M << "/" << m[j] << "= " << M / m[j] << endl;
        Mi[j] = M / m[j];
    }
    cout << "Mi' = Multiplicative inverse of Mi under mod mi\n";
    for (auto j = 0; j < m.size(); j++) {
        cout << "M" << j + 1 << "' = MI(" << Mi[j] << "," << m[j] << "): ";
        Minv[j] = MI_EEA(Mi[j], m[j]);
    }
    cout << "Now substituting values in above eqn X = (a1M1M1' + a2M2M2' +..) mod M  \nX = ";
    for (auto j = 0; j < Mi.size(); j++) cout << "+ ( " << a[j] << " * " << Mi[j] << " * " << Minv[j] << " )";
    int X = 0;
    cout << endl << "X = ";
    for (auto j = 0; j < Mi.size(); j++) {
        M = a[j] * Mi[j] * Minv[j];
        X += M;
        cout << "+(" << M << ") ";
    }
    cout << "\nX = " << X;
    return X;
}

