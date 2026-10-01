"""Decompress the xz-compressed certificates cert_tracial_d{13..20}.pkl.xz and check their SHA-256 sums.

Usage (from this folder):  python decompress.py
Writes cert_tracial_d{d}.pkl next to the .xz files and verifies them against SHA256SUMS_d13-20.txt.
Only the Python standard library is needed.
"""
import hashlib
import lzma
import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
sums = {}
with open(os.path.join(here, 'SHA256SUMS_d13-20.txt')) as f:
    for line in f:
        h, name = line.split()
        sums[name] = h

ok = True
for name, h in sorted(sums.items()):
    xz = os.path.join(here, name + '.xz')
    out = os.path.join(here, name)
    if not os.path.exists(out):
        with open(xz, 'rb') as fi, open(out, 'wb') as fo:
            fo.write(lzma.decompress(fi.read()))
    got = hashlib.sha256(open(out, 'rb').read()).hexdigest()
    status = 'OK' if got == h else 'MISMATCH'
    ok &= got == h
    print(f'{name}: sha256 {got[:16]}... {status}')
sys.exit(0 if ok else 1)
