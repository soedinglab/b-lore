import numpy as np


def unique_matrix(a):
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx = np.unique(b, return_index=True)
    idx_sorted = np.sort(idx)
    unique_a = a[idx_sorted]
    return unique_a, idx_sorted
