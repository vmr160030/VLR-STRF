import numpy as np

def build_vlrPrior(name: str, nk: np.ndarray, **kwargs):
    """
    Initialise structure for variational low-rank strf prior.

    Args:
        name (str): 'ASD','ALD','TRD','RR'
        nk (np.ndarray): [nk1 nk2] vector with RF dimensionality
        hprs (dict): initial hyperparameter vector, see kernel functions
%                        for detials on order 
    """