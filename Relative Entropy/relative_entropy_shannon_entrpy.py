### This code calculates the relative entropy for Shannon entropy values

import numpy as np
from sklearn.neighbors import KernelDensity

def kde_kl_estimate(samples_q, samples_p, bandwidth=None, kernel='gaussian'):  # samples_q & p: arrays of shape (n_samples, d)
    if bandwidth is None:                             # Choose bandwidth by rule-of-thumb if not given (Silverman)
        n, d = samples_q.shape
        std = np.std(samples_q, axis=0).mean()
        bandwidth = 1.06 * std * n ** (-1 / (d + 4))  # Silverman's rule

    kde_q = KernelDensity(bandwidth=bandwidth, kernel=kernel).fit(samples_q)
    kde_p = KernelDensity(bandwidth=bandwidth, kernel=kernel).fit(samples_p)

    logq = kde_q.score_samples(samples_q)             # log density at q-samples
    logp_at_q = kde_p.score_samples(samples_q)        # log p evaluated at q-samples

    if np.any(np.isneginf(logp_at_q)):                # If any logp_at_q is -inf (zero density) raise warning
        print("Warning: p_hat ~ 0 at some q samples — KL estimate may be +inf or unstable.")

    kl_hat = np.mean(logq - logp_at_q)
    kl_hat_2 = np.log2(np.exp(kl_hat))
    return kl_hat_2


# MAIN
with open("shannon_entropy_natural_L20.txt", 'r') as f:
    entropies_natural = [(float(line.strip().split()[1])) for line in f.readlines()]  

with open("shannon_entropy_rand_L20.txt", 'r') as f:
    entropies_random = [(float(line.strip().split()[1])) for line in f.readlines()]

# Convert to numpy arrays with shape (n_samples, 1) for sklearn
natural_data = np.array(entropies_natural).reshape(-1, 1)
rand_data = np.array(entropies_random).reshape(-1, 1)

print(f'{kde_kl_estimate(natural_data, rand_data):.2f}')