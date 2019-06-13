import ROOT,os

if __name__ == "__main__":

    fileOut = ROOT.TFile.Open('zpt_weights.root','recreate')
    fileIn = ROOT.TFile.Open("zll_ressum_ratio_SMP-17-10_Fig6.root","read")
    c = fileIn.Get("c_3_ratio")
    primitives  = c.GetListOfPrimitives()
    
    hname = {'c_3_ratio_3':'zll_ressum_ratio_h_1',
             'c_3_ratio_2':'zmmPtSplit_h_1',
             'c_3_ratio_1':'zmmPtSplit_h_1'}

    hnameNew = {'c_3_ratio_3':'amcnlo',
                'c_3_ratio_2':'resbos',
                'c_3_ratio_1':'geneva'}

    ratios = {}
    for p in primitives:
        padname = p.GetName()
        print "pad = ", padname
        for ratio in p.GetListOfPrimitives():
            print "\t",ratio
            if ratio.GetName()==hname[padname]:
                 if padname not in ratios: 
                     ratios[padname] = ratio.Clone(hnameNew[padname])
                     ratios[padname].SetDirectory(None)
                 
    fileOut.cd()
    for k,r in ratios.iteritems():
        r.Write()
    fileOut.Close()


