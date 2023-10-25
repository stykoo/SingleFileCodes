#ifndef UTILS_H
#define UTILS_H

// Compute a^b with b a positive integer
// From https://stackoverflow.com/a/18734205
template<typename T, typename U>
T mypow(const T base, U degree) {
    T result = 1;
    T term = base;
    while (degree)
    {
        if (degree & 1)
            result *= term;
        term *= term;
        degree = degree >> 1;
    }
    return result;
}  

// Wrap x into [0,L).
template<typename T>
long periodicBC(const T x, const T L) {
	return (x + L) % L;
	//return (x < 0) ? (x % L + L) : (x % L);
}

// Wrap x into (-L/2,L/2].
// TODO improve
template<typename T>
long periodicBCsym(const T x, const T L) {
	long y = periodicBC(x, L);
	return (y > L / 2) ? (y - L) : y;
}

// Add 1 in periodic boundary conditions (assuming 0 <= x < L)
template<typename T>
long periodicAdd1(const T x, const T L) {
	return (x + 1) % L;
}

// Substract 1 in periodic boundary conditions (assuming 0 <= x < L)
template<typename T>
long periodicSubs1(const T x, const T L) {
	return (x - 1 + L) % L;
}

// Increment in periodic boundary conditions (assuming 0 <= x < L)
template<typename T>
void periodicInc(T &x, const T L) {
	if (x < L - 1) {
		x++;
	} else {
		x = 0;
	}
}

// Decrement in periodic boundary conditions (assuming 0 <= x < L)
template<typename T>
void periodicDec(T &x, const T L) {
	if (x > 0) {
		x--;
	} else {
		x = L - 1;
	}
}

#endif
