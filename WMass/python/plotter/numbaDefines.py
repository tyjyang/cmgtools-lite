import ROOT

@ROOT.Numba.Declare(["RVec<float>"], "float")
def sumPt(pts):
    sumpt = 0
    for pt in pts:
        sumpt += pt
    return sumpt

