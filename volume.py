import pandas as pd
from density import get_density

molecular_weights = pd.read_csv("molecularWeights.txt", sep='!', names=['Compounds', 'MW'])
molecular_weights = dict(zip(molecular_weights['Compounds'], molecular_weights['MW']))


def get_volume(s, T):
    return 0.001 * molecular_weights[s] / get_density(s, T)
