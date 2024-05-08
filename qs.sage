# Returns a heuristic smoothness bound B
def smoothness_bound(N):
    L = e^sqrt(ln(N) * ln(ln(N)))
    B = L^(1/sqrt(2))
    return next_prime(ceil(B))

# Returns true if n is a quadratic residue mod p (p odd)
def euler_criterion(n, p):
    return mod(n, p)^((p-1)/2)

# Returns a list of primes up to B which are quadratic residues mod N
def generate_factor_base(N, B):
    factor_base = list()
    factor_base.append(2) # Euler criterion only works with odd primes, all integers are quad residues of 2

    # This loop checks if N is a quadratic residue for each prime p
    for p in prime_range(3, B+1):
        if euler_criterion(N, p) == 1:
            factor_base.append(p)
    
    return factor_base

# Constructs a sieve using the chosen polynomial Y(X) = (X + ceil(sqrt(N))^2 - N
def construct_sieve(N):
    sieve = list()
    n_sqrt = ceil(sqrt(N))

    for i in range(0, 100): # Abitrary bound
        sieve.append((i + n_sqrt)^2 - N)

    return sieve

# Uses Tonelli-Shanks algorithm to solve relation of the form X^2 === N (mod p)
def tonelli(n, p):
    q = p - 1
    s = 0
    while mod(q, 2) == 0:
        q //= 2
        s += 1
    if s == 1:
        r = mod(n, p)^((p + 1) // 4)
        return r, p - r
    for z in range(2, p):
        if p - 1 == euler_criterion(z, p):
            break
    c = mod(z, p)^q
    r = mod(n, p)^((q + 1) // 2)
    t = mod(n, p)^q
    m = s
    t2 = 0
    while mod(t - 1, p) != 0:
        t2 = mod(t * t, p)
        for i in range(1, m):
            if mod(t2 - 1, p) == 0:
                break
            t2 = mod(t2 * t2, p)
        b = mod(c, p)^(1 << (m - i - 1))
        r = mod(r * b, p)
        c = mod(b * b, p)
        t = mod(t * c, p)
        m = i

    return (r, p-r)

# Sieves value list for b-smooth values
def find_smooth(N, B, fb):
    values = construct_sieve(N)
    original_values = values.copy()
    smooth_values = list()

    # First check if we have to sieve by 2
    # If so, find starting index and sieve every other value by 2
    if fb[0] == 2:
        start = mod(isqrt(N) - ceil(sqrt(N)), 2)
        for i in range(start, len(values), 2):
            values[i] //= 2

    # For every p in the factor base, we sieve every pth value of the sieve starting from the values returned by Tonelli-Shanks
    for f in fb[1:]:
        residues = tonelli(N, f)
        indices = [mod(residues[0] - ceil(sqrt(N)), f), mod(residues[1] - ceil(sqrt(N)), f)]
        
        for i in indices:
            for j in range(i, len(values), f):
                values[j] //= f

    # Create a list of the found smooth values
    for i in range(0, len(values)):
        if values[i] == 1:
            smooth_values.append(original_values[i])

    return smooth_values

# Given factor base and list of smooth values, builds a matrix of exponent vectors
def build_matrix(smooth_values, fb):
    m = matrix(Zmod(2), len(smooth_values), len(fb))

    for i in range(0, len(smooth_values)):
        for j in range(0, len(fb)):
            m[i, j] = mod(valuation(smooth_values[i], fb[j]), 2)

    return m

# Returns the left nullspace of the given matrix m
def solve_matrix(m):
    soln = m.kernel()
    return list(soln[1])

# Using the matrix solution and smooth values, find factors of N
def find_factors(N, smv, soln):
    a = 1
    # Multiply smooth values together to get one side of a congruence of squares
    for v in range(0, len(smv)):
        if soln[v] == 1:
            a *= smv[v]
    a = isqrt(a)

    # Re-calculate initial values from polynomial in the form of X + 124
    vals = list()
    for sm in smv:
        vals.append(isqrt(N + sm))
    
    b = 1
    # Multiply initial values together to find the other side of the congruence of squares
    for v in vals:
        b *= v

    # Calculate GCD with N of each resulting part of the difference of squares
    f1 = gcd(b + a, N)
    f2 = gcd(b - a, N)
    
    return [f1, f2]

def qs(N):
    print(f"Factoring {N}")

    # Choose a smoothness bound
    B = smoothness_bound(N)
    print(f"Smoothness bound: {B}")

    # Use sieving to locate prime_pi(B) + 1 B-smooth numbers
    sieve = construct_sieve(N)

    # Factor the smooth numbers and generate exponent vectors mod 2
    fb = generate_factor_base(N, B)
    print(f"Factor base: {fb}")

    smv = find_smooth(N, B, fb)
    print(f"Smooth values: {smv}")

    m = build_matrix(smv, fb)
    print(f"Matrix of exponent vectors:\n{m}")

    # Find a subset of the exponent vectors that adds to the zero vector
    soln = solve_matrix(m)
    print(f"Solution to matrix: {soln}")

    # Use difference of squares identity and Euclidean algorithm to find a factor
    factors = find_factors(N, smv, soln)
    print(f"{N} = {factors[0]} * {factors[1]}")