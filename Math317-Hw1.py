#Scratch code for HW1 of Math 317

import math

def subst_encrpy(plaintext, key):
    L = list()
    for char in plaintext:
        L.append(chr(key[ord(char) - 65] + 65))
    return str().join(L)

def subst_decrypt(ciphertext, key):
    L = list()
    for T in ciphertext:
        L.append(chr(key.index(ord(T) - 65) + 65))
    return str().join(L)

def is_prime(num):
    if (num > 1):
        for i in range(2,int(math.sqrt(num) + 1)):
            if num % i == 0:
                return False
    else:
        return False
    return True

#Counts primes below n where p ≅ b(mod m)
def count_congruent_primes_below(n, b, m):
    total = 0
    for i in range(0,n):
        if(is_prime(i) and (i % m == b)):
            total += 1
    print("Total number of primes below " + str(n) + " congruent to "+ str(b) + "(mod " + str(m) + "): " + str(total))

def main():
    '''
    key = [15, 21, 25, 17, 5, 16, 22, 13, 19, 23, 18, 2, 12, 20, 6, 3, 1, 14, 24, 7, 8, 10, 0, 9, 4, 11]
    M = 'ILOVESPAGHETTI'
    C = subst_encrpy(M, key)
    print(C)

    P = subst_decrypt(C, key)
    print(P)
    

    for i in range(0,50):
        if is_prime(i):
            print(i)
    '''
    count_congruent_primes_below(100, 1, 4)

    count_congruent_primes_below(100, 3, 4)

if __name__ == "__main__":
    main()