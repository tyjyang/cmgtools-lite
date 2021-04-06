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
    corrbin = np.digitize(np.array([ptW]), binsWp, right=True)
    if corrbin[0] >= len(correctionsWp):
        return 1.
    return correctionsWp[corrbin][0]

@ROOT.Numba.Declare(["RVec<int>","RVec<int>","RVec<int>","RVec<int>", "RVec<float>", ], "RVec<bool>")
def prefsrLeptons(status, statusFlags, pdgId, motherIdx, pts):
    leptons = (np.abs(pdgId) >= 11) & (np.abs(pdgId) <= 14) & (motherIdx >= 0)
    status746 = status == 746
    status23 = status == 23
    motherV = (pdgId[motherIdx] == 23) | (np.abs(pdgId[motherIdx]) == 24)
    fromHardProcess = ((statusFlags >> 8 ) & 1).astype(np.bool_)

    # Some leptons in MadGraph have no W/Z in the history, but have status 23
    others = leptons & (motherV | status23)
    if np.count_nonzero(others) > 2:
        others = others & fromHardProcess
    
    # If there are status = 746 leptons, they came from photos and are pre-FSR
    photos = leptons & status746
    photos = status746
    nphotos = np.count_nonzero(photos)
    nothers = np.count_nonzero(others)

    if nphotos == 2:
        return photos
    elif nphotos == 1 and nothers == 1:
        return photos | others

    return others

