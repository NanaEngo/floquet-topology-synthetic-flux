"""Derive Fock-space weights of the prepared even cat state (alpha = 1).

The production initial state is the normalized even cat |alpha>+|-alpha> with
alpha = 1 placed on the first site. This script computes the exact Fock
weights of that state and the weight discarded by the production Fock
truncation N_max = 2, writing data/cat_fock_weights.json.

The script is self-checking and aborts on any mismatch (weights must sum to
one within 1e-12 and must reproduce the reference values at 1e-6).
"""

import json
import math
from pathlib import Path

ALPHA = 1.0
N_MAX_PRODUCTION = 2
N_CUTOFF = 40
OUT = Path(__file__).resolve().parents[1] / "data" / "cat_fock_weights.json"

# even-cat amplitudes: c_n = 2 exp(-a^2/2) a^n / sqrt(n!) / sqrt(2(1+exp(-2a^2))),
# nonzero for even n only.
norm2 = 2.0 * (1.0 + math.exp(-2.0 * ALPHA**2))
amp = {}
for n in range(0, N_CUTOFF, 2):
    amp[n] = 2.0 * math.exp(-0.5 * ALPHA**2) * ALPHA**n / math.sqrt(
        math.factorial(n)
    )
total = sum(a * a for a in amp.values())
assert abs(total - norm2) < 1e-12, (total, norm2)

weights = {str(n): (amp[n] ** 2) / total for n in amp}
w_n0 = weights["0"]
w_n2 = weights["2"]
w_n4 = weights["4"]
w_ge4 = sum(weights[str(n)] for n in amp if n >= 4)
w_truncated = sum(weights[str(n)] for n in amp if n > N_MAX_PRODUCTION)
w_kept = sum(weights[str(n)] for n in amp if n <= N_MAX_PRODUCTION)

assert abs(w_kept + w_truncated - 1.0) < 1e-12
assert abs(w_ge4 - w_truncated) < 1e-12
assert abs(w_n0 - 0.648055) < 1e-6, w_n0
assert abs(w_n2 - 0.324027) < 1e-6, w_n2
assert abs(w_n4 - 0.027002) < 1e-6, w_n4
assert abs(w_ge4 - 0.027919) < 1e-6, w_ge4

payload = {
    "status": "CAT_FOCK_WEIGHTS_DERIVED",
    "scope": "exact Fock decomposition of the prepared even cat; analytic, no sampling",
    "design": {
        "alpha": ALPHA,
        "state": "normalized even cat |alpha>+|-alpha> on site 0",
        "production_truncation": N_MAX_PRODUCTION,
        "cutoff_for_sum": N_CUTOFF,
    },
    "weights": weights,
    "aggregates": {
        "n0": w_n0,
        "n2": w_n2,
        "n4": w_n4,
        "n_ge_4_discarded_by_nmax2": w_truncated,
        "n_le_2_kept_by_nmax2": w_kept,
    },
    "software": {"python": "3.10+"},
}
OUT.write_text(json.dumps(payload, indent=1) + "\n")
print(f"wrote {OUT}")
print(json.dumps(payload["aggregates"], indent=1))
