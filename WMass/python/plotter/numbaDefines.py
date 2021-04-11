import ROOT
import numpy as np
import numba
import os
script_dir = os.path.dirname(os.path.realpath(__file__))

@ROOT.Numba.Declare(["RVec<float>"], "float")
def maxPt(pts):
    return np.max(pts) if len(pts) else 0

correctionsWp_file = np.load(f"{script_dir}/../postprocessing/data/gen/fiducialWp_RatioN3LL.npz")
binsWp = correctionsWp_file['bins']
correctionsWp = correctionsWp_file['hist']

correctionsZ_file = np.load(f"{script_dir}/../postprocessing/data/gen/fiducialDY_RatioN3LL.npz", allow_pickle=True)
binsZ = correctionsZ_file['bins']
binsZ_m = binsZ[0]
binsZ_y = binsZ[1]
binsZ_pt = binsZ[2]
correctionsZ = correctionsZ_file['scetlibCorr3D_Z']

#@numba.jit(nopython=True)
#def correctN3LL(mV, yV, ptV, bins, corr):
@ROOT.Numba.Declare(["double", "double", "double"], "float")
def correctN3LL_Z(mV, yV, ptV):
    # Numba doesn't seem to have digitize with the first arg a scalar implemented
    binm = np.digitize(np.array([mV]), binsZ_m)[0]
    biny = np.digitize(np.array([yV]), binsZ_y)[0]
    binpt = np.digitize(np.array([ptV]), binsZ_pt)[0]
    if binm == 0 or binm == len(binsZ_m) or biny == 0 or biny == len(binsZ_y) or binpt == 0 or binpt == len(binsZ_pt):
        return 1.
    return correctionsZ[binm-1,biny-1,binpt-1]

#@ROOT.Numba.Declare(["float", "float", "float"], "float")
#def correctN3LL_Wp(mW, yW, ptW):
#    return correctN3LL(mW, yW, ptW, binsWp, correctionsWp)

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

