from scipy.stats import gmean
from sklearn.metrics import r2_score
import numpy as np


def get_alpha_gm(x, y):
    alpha = np.divide((y * (1 - x)) , (x * (1 - y)))
    alpha_gm = gmean(alpha)
    if alpha_gm < 1:
        alpha_gm = 1/alpha_gm
    y_alpha = alpha_gm * x / (1 + (alpha_gm - 1) * x)
    acc = r2_score(y, y_alpha)
    if acc < 0.80:
        alpha_gm = 0
    return alpha_gm

