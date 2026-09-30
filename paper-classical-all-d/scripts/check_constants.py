"""Check of the explicit constants in the lemmas "junction function" and "junction inequality" of the paper.

Everything is evaluated in interval arithmetic (mpmath.iv, 60 decimal digits); the zeta values are enclosed
rigorously: zeta(2) = pi^2/6 and zeta(4) = pi^4/90 exactly, and for s = 3, 5 the partial sum up to M plus the
integral tail bounds  1/((s-1)(M+1)^(s-1)) <= zeta(s) - sum_{m<=M} m^-s <= 1/((s-1) M^(s-1)).

Claims checked:
 (C1) |eps(1)| = 2 log 2 - 1 and |eps(5)| = |1/5 - (6 log 6 - 10 log 5 + 4 log 4)|;  |eps(1)| - |eps(5)| >= 0.38493.
 (C2) tail constant  T = sum_{i>=3} 4^(3-i)/(i(2i-1)) < 0.0773 < 4/45  (T enclosed by a partial sum + geometric tail).
 (C3) c_sc = (zeta(2) - zeta(3))/6 + (4/45)(zeta(4) - zeta(5)) <= 0.07785,
      c_fw = (zeta(2) - 1)/6 + (4/45)(zeta(4) - 1)             <= 0.11481.
 (C4) final junction bound  -(w-1)*0.38493 + 0.07785*w + 0.11481 = 0.49974 - 0.30708 w < 0 for w >= 2.
 (C5) numerical sanity checks: series vs closed form of eps(u); |eps(u)| <= 1/(6u^3) + (4/45) u^-5 for 2 <= u <= 10^4;
      |eps| strictly decreasing on 1..10^4; the limits sum_{m>=2} (m-1)|eps(m)| and sum_{m>=2} m|eps(m)|.
"""
from mpmath import mp, iv, mpf, log, nsum, inf

iv.dps = 60
mp.dps = 40


def zeta_enclosure(s, M=2000):
    part = iv.mpf(0)
    for m in range(1, M + 1):
        part += iv.mpf(1) / iv.mpf(m) ** s
    lo = part + iv.mpf(1) / ((s - 1) * iv.mpf(M + 1) ** (s - 1))
    hi = part + iv.mpf(1) / ((s - 1) * iv.mpf(M) ** (s - 1))
    return lo.hull(hi) if hasattr(lo, "hull") else iv.mpf([lo.a, hi.b])


def eps_closed_iv(u):
    u = iv.mpf(u)
    return 1 / u - ((u + 1) * iv.log(u + 1) - 2 * u * iv.log(u) + (u - 1) * iv.log(u - 1))


if __name__ == "__main__":
    # (C1)
    e1 = 2 * iv.log(iv.mpf(2)) - 1
    e5 = -eps_closed_iv(5)
    diff = e1 - e5
    print(f"(C1) |eps(1)| = 2 log 2 - 1 in {iv.nstr(e1, 12)};  |eps(5)| in {iv.nstr(e5, 12)};  "
          f"|eps(1)| - |eps(5)| in {iv.nstr(diff, 12)}  (>= 0.38493: {bool(diff >= mpf('0.38493'))})")
    # (C2)  partial sum i = 3..40 plus tail <= (1/(41*81)) * 4^(3-41) * 4/3
    T = iv.mpf(0)
    for i in range(3, 41):
        T += iv.mpf(4) ** (3 - i) / (i * (2 * i - 1))
    tail = iv.mpf(4) ** (3 - 41) / (41 * 81) * iv.mpf(4) / 3
    T = iv.mpf([T.a, (T + tail).b])      # enclosure of the full series (all terms positive)
    print(f"(C2) T = sum_(i>=3) 4^(3-i)/(i(2i-1)) in {iv.nstr(T, 12)};  T < 0.0773: "
          f"{bool(T < mpf('0.0773'))};  0.0773 < 4/45 = {mp.nstr(mpf(4)/45, 12)}: {bool(mpf('0.0773') < mpf(4)/45)}")
    # (C3)
    pi = iv.pi
    z2, z4 = pi ** 2 / 6, pi ** 4 / 90
    z3, z5 = zeta_enclosure(3), zeta_enclosure(5)
    csc_ = (z2 - z3) / 6 + iv.mpf(4) / 45 * (z4 - z5)
    cfw = (z2 - 1) / 6 + iv.mpf(4) / 45 * (z4 - 1)
    print(f"(C3) zeta(3) in {iv.nstr(z3, 15)}, zeta(5) in {iv.nstr(z5, 15)}")
    print(f"     c_sc in {iv.nstr(csc_, 12)}  (<= 0.07785: {bool(csc_ <= mpf('0.07785'))});  "
          f"c_fw in {iv.nstr(cfw, 12)}  (<= 0.11481: {bool(cfw <= mpf('0.11481'))})")
    # (C4)
    ok = True
    for w in range(2, 200):
        val = -(w - 1) * mpf('0.38493') + mpf('0.07785') * w + mpf('0.11481')
        ok &= val < 0
    print(f"(C4) -(w-1)*0.38493 + 0.07785 w + 0.11481 = 0.49974 - 0.30708 w;  at w = 2: "
          f"{mp.nstr(-(1) * mpf('0.38493') + mpf('0.07785') * 2 + mpf('0.11481'), 8)};  negative for w = 2..199: {ok} "
          f"(linear and decreasing in w, hence negative for all w >= 2)")
    # (C5)
    def eps_closed(u):
        u = mpf(u)
        if u == 1:
            return 1 - 2 * log(2)
        return 1 / u - ((u + 1) * log(u + 1) - 2 * u * log(u) + (u - 1) * log(u - 1))

    def eps_series(u):
        u = mpf(u)
        return -nsum(lambda i: u ** (-(2 * i - 1)) / (i * (2 * i - 1)), [2, inf])
    ser = max(abs(eps_closed(u) - eps_series(u)) for u in (1, 2, 3, 5, 10, 37, 100))
    bound_ok = all(abs(eps_closed(u)) <= 1 / (6 * mpf(u) ** 3) + mpf(4) / 45 / mpf(u) ** 5 for u in range(2, 10001))
    mono = all(abs(eps_closed(u + 1)) < abs(eps_closed(u)) for u in range(1, 10000))
    lim_sc = nsum(lambda m: (m - 1) * abs(eps_closed(int(m))), [2, inf])
    lim_fw = nsum(lambda m: m * abs(eps_closed(int(m))), [2, inf])
    print(f"(C5) series = closed form (max diff {mp.nstr(ser, 3)}); |eps(u)| <= 1/(6u^3) + (4/45)u^-5 for u = 2..10^4: "
          f"{bound_ok}; |eps| strictly decreasing on 1..10^4: {mono};")
    print(f"     sum_(m>=2) (m-1)|eps(m)| = {mp.nstr(lim_sc, 10)} (<= c_sc),  sum_(m>=2) m|eps(m)| = {mp.nstr(lim_fw, 10)} (<= c_fw)")
