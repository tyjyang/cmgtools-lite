import ROOT
import numpy as np
import os
script_dir = os.path.dirname(os.path.realpath(__file__))

@ROOT.Numba.Declare(["RVec<float>"], "float")
def maxPt(pts):
    return np.max(pts) if len(pts) else 0

correctionsWp_file = np.load(f"{script_dir}/../postprocessing/data/gen/fiducialWp_RatioN3LL.npz")
binsWp = correctionsWp_file['bins']
correctionsWp = correctionsWp_file['hist']

@ROOT.Numba.Declare(["float"], "float")
def correctN3LL_Wp(ptW):
    corrbin = np.digitize(np.array([ptW]), binsWp)
    if corrbin[0] >= len(correctionsWp):
        return 1.
    return correctionsWp[corrbin][0]
