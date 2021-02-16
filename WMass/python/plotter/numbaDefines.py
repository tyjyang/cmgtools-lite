import ROOT
import numpy as np

@ROOT.Numba.Declare(["RVec<float>"], "float")
def maxPt(pts):
    return np.max(pts) if len(pts) else 0
