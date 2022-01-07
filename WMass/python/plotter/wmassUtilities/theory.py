pdfMap = {
    "nnpdf31" : {
        "name" : "NNPDF31",
        "entries" : 103,
        "onlyW" : True, # Because the alpha_s sets are missing
        "truncate" : False,
        "weight" : "LHEPdfWeight",
        "alphas" : ["LHEPdfWeightAltSet5[0]", "LHEPdfWeightAltSet6[0]"],
        "alphaRange" : "001",
    },
    "nnpdf30" : {
        "name" : "NNPDF30",
        "entries" : 101,
        "onlyW" : True,
        "truncate" : False,
        "weight" : "LHEPdfWeightAltSet13",
        "alphas" : ["LHEPdfWeightAltSet15[0]", "LHEPdfWeightAltSet16[0]"],
        "alphaRange" : "001",
    },
    "ct18" : {
        "name" : "CT18",
        "entries" : 59,
        "onlyW" : True,
        "truncate" : True,
        "weight" : "LHEPdfWeightAltSet18",
        "alphas" : ["LHEPdfWeightAltSet18[59]", "LHEPdfWeightAltSet18[60]"],
        "alphaRange" : "002",
    },
    "mmht" : {
        "name" : "MMHT",
        "entries" : 51,
        "onlyW" : True,
        "truncate" : False,
        "weight" : "LHEPdfWeightAltSet19",
        "alphas" : ["LHEPdfWeightAltSet20[1]", "LHEPdfWeightAltSet20[2]"],
        "alphaRange" : "002",
    },
}