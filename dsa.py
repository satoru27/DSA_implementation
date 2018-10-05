import gmpy2 as gmp
import random
import hashlib
from gmpy2 import mpz

L = 2048
N = 256

#https://en.wikibooks.org/wiki/Algorithm_Implementation/Mathematics/Extended_Euclidean_algorithm


def egcd(a, b):
    if a == 0:
        return b, 0, 1
    else:
        g, y, x = egcd(b % a, a)
        return g, x - (b // a) * y, y


def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m


def no_bits(n):
    return len(bin(n))-2


def print_info(item_lst, name_lst):
    for item, name in zip(item_lst, name_lst):
        print(f'{name} = {item}')
        print(f'{name} type = {type(item)}')
        print(f'{name} bits = {no_bits(item)}\n')


def generate_q():
    q = gmp.mpz_urandomb(gmp.random_state(random.randint(0, 367263292)), N)
    while not gmp.is_prime(q):
        q = gmp.next_prime(q)
        if no_bits(q) != N:
            q = gmp.mpz_urandomb(gmp.random_state(random.randint(0, 367263292)), N)

    return q


def generate_p(q):
    # p = gmp.mpz_urandomb(gmp.random_state(random.randint(0, 367263292)), L)
    # p = gmp.next_prime(p)
    p = mpz(2**L)
    p = gmp.next_prime(p)
    i = 0

    while ((p - 1) % q) != 0:
        p = gmp.next_prime(p)
        # if no_bits(p) != L:
        #     p = gmp.mpz_urandomb(gmp.random_state(random.randint(0, 367263292)), L)
        #     p = gmp.next_prime(p)

    return p


def generate_g(q, p):
    g = 1
    n = mpz(p-1)
    m = mpz(p-1//q)
    while g == 1:
        h = mpz(random.randint(1, n))
        g = gmp.powmod(h, m, p)

    return g


def generate_keys(q, p, g):
    # generates private key x and public key y
    x = gmp.mpz_random(gmp.random_state(random.randint(0, 367263292)), q)
    while x == 0:
        x = gmp.mpz_random(random.randint(0, 367263292), q)

    y = gmp.powmod(g, x, p)

    return x, y


def sign(m, g, x, p, q):
    k = gmp.mpz_random(gmp.random_state(random.randint(0, 367263292)), q)
    r = gmp.powmod(g, k, p) % q

    if r == 0:
        while r != 0:
            k = gmp.mpz_random(gmp.random_state(random.randint(0, 367263292)), q)
            r = gmp.powmod(g, k, p) % q

    inv_k = modinv(k, q)

    # s = inv_k * (H(m) + xr) mod q

    h_m = hashlib.sha3_256()
    h_m.update(m.encode('utf-8'))
    h_m = mpz(h_m.hexdigest(), base=16)

    s = (inv_k * (h_m + x*r)) % q

    return r, s


def verify(r, s, q, m, y, g, p):
    w = modinv(s, q)
    # test
    test_var = w * s % q

    h = hashlib.sha3_256()
    h.update(m.encode('utf-8'))
    h = mpz(h.hexdigest(), base=16)

    u1 = (h*w) % q

    u2 = (r*w) % q

    # v = (g^u1 * y^u2 mod p) mod q

    v = ((gmp.powmod(g, u1, p) * gmp.powmod(y, u2, p)) % p) % q

    return v == r


def main():
    q = generate_q()
    p = generate_p(q)
    g = generate_g(q, p)

    print_info([q, p, g], ['q', 'p', 'g'])

    #
    # x, y = generate_keys(q, p, g)
    #
    # print_info([x, y], ['x', 'y'])
    #
    # # k, r, inv_k = sign("Hello World!",g , x, p, q)
    # #
    # # print_info([k, r, inv_k], ['k', 'r', 'inv_k'])
    # #
    # # print(f'RESULT k*inv_k mod q = {k*inv_k % q}')
    #
    # m = 'Oh God Why'
    #
    # r, s = sign(m, g, x, p, q)
    #
    # print_info([r, s], ['r', 's'])
    #
    # print('******Verification 1******')
    # print(verify(r, s, q, m, y, g, p))
    #
    # print('******Verification 2******')
    # m2 = 'Oh God Why Me'
    # print(verify(r, s, q, m2, y, g, p))
    #
    # print('******Verification 3******')
    # print(verify(r, s, q, m, y+1, g, p))


if __name__ == "__main__":
    main()
