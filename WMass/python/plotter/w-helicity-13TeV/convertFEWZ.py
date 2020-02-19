import ROOT, itertools, math
from array import array

def getHelicity(charge, rightleft, a0, a4):

    sign  = -1. if 'p' in charge    else 1.
    rsign = -1. if 'r' in rightleft else 1.

    retval = array('d', [])
    retnum = array('d', [])
    retpdf = array('d', [])
    rettot = array('d', [])

    ## fl = 0.25 * (2. - a0 + sign*a4)
    ## fr = 0.25 * (2. - a0 - sign*a4)

    for i,ii in enumerate(a0[1]):
        retval .append( 0.25 * (2. - ii + sign*rsign*a4[1][i]) )
        retnum .append( 0.)#math.sqrt( a0[4][i]**2 + a4[4][i]**2 ) )
        retpdf .append( 0.)#math.sqrt( a0[5][i]**2 + a4[5][i]**2 ) )
        rettot .append( 0.)#math.sqrt( a0[6][i]**2 + a4[6][i]**2 ) )
        #print charge, rightleft, 'rap', a0[0][i], 'value', retval[-1], '+-', rettot[-1]

    #print ch, rightleft, retval, rettot

    return a0[0], retval, a0[2], a0[3], retnum, retpdf, rettot

    

def returnArray(lines, ipdf):

    retx   = array('d', [])
    retval = array('d', [])
    retnum = array('d', [])
    retpdf = array('d', [])
    rettot = array('d', [])

    retxup = array('d', [])
    retxdn = array('d', [])

    tmp_iy = 0
    for line in lines:
        tmp_line = line.strip()

        if not tmp_line[0].isdigit():
            continue

        if not tmp_iy:
            tmp_iy = tmp_line[0]

        retx   .append( float(tmp_line.split()[0]) )
        retval .append( float(tmp_line.split()[1]) )
        retnum .append( float(tmp_line.split()[2]) )
        retpdf .append( float(tmp_line.split()[3]) )
        rettot .append( math.sqrt(retnum[-1]**2 + retpdf[-1]**2) )
        retxup.append(0.)
        retxdn.append(0.)

    

    return retx, retval, retxdn, retxup, retnum, retpdf, rettot
        

def addPdfIndex(lines,offset=1):

    newlines = []
    tmp_iy = -1
    for iline, line in enumerate(lines):
        if line[0] != tmp_iy:
            tmp_iy = line[0]
            ntmp = 0
        newlines.append([ntmp+offset]+[float(i) for i in line])
        ntmp += 1

    return newlines
    
def divideByDenom(lines, denlines):
    newlines = []

    for il, line in enumerate(lines):
        newlines.append([line[0], line[1], line[2]/denlines[il][2], line[3]/denlines[il][3]])

    return newlines


def getGraph(pdfset, centrallines, errorlines):

    xvals = array('d', [i[1] for i in centrallines])
    xerrs = array('d', [0.12 for i in centrallines])

    yvals = array('d', [i[2] for i in centrallines])

    ups, dns = [], []
    if pdfset == 'nnpdf31':
        for central in centrallines:
            cen = central[2]

            allerrs = [i[2] for i in errorlines if i[1] == central[1] ]
            pdferrs = allerrs[:-2]
            aserrs  = allerrs[-2:]
            pdferr  = math.sqrt(sum( [ (cen-err)**2 for err in pdferrs] ) )

            ups.append(math.sqrt(pdferr**2 + ( (cen-aserrs[0])*0.67)**2 + central[3]**2) )
            dns.append(math.sqrt(pdferr**2 + ( (cen-aserrs[1])*0.67)**2 + central[3]**2) )

    elif pdfset in ['ct18', 'hera']:
        for central in centrallines:
            cen = central[2]

            allerrs = [i[2] for i in errorlines if i[1] == central[1] ]

            pairs = zip(allerrs[0::2], allerrs[1::2])
            pairsup = allerrs[0::2]
            pairsdn = allerrs[1::2]
            upstmp = [max(0, x[0] - cen, x[1] - cen) for x in pairs]
            dnstmp = [max(0, cen - x[0], cen - x[1]) for x in pairs] 
            #upstmp = [cen - x for x in pairsup]
            #dnstmp = [x - cen for x in pairsdn] 

            ups.append(math.sqrt(sum( [(x/1.645)**2 for x in upstmp] ) + central[3]**2))
            dns.append(math.sqrt(sum( [(x/1.645)**2 for x in dnstmp] ) + central[3]**2))

    arr = []
    for ii, v in enumerate(xvals):
        #print v, yvals[ii], ups[ii], dns[ii]
        arr.append([v, yvals[ii], ups[ii], dns[ii]])

    return arr
    

def getForPDF(pdfset): ## pdfset can be 'ct18' or 'nnpdf31'
    ## basedir = '/afs/cern.ch/work/a/arapyan/public/fewz_smp18012/'
    basedir = '/eos/cms/store/group/phys_smp/arapyan/fewz_prediction/'
    
    outfile_name = '/afs/cern.ch/work/m/mdunser/public/cmssw/w-helicity-13TeV/CMSSW_8_0_25/src/CMGTools/WMass/data/theory/theory_cross_sections_FEWZ.root'
    outfile = ROOT.TFile(outfile_name, 'RECREATE')

    vals = {}
    
    for ch, pset in itertools.product('pm', [pdfset]):
    
    
        tmp_infile = open(basedir+'/my_w{ch}_18012_{pset}/NNLO.pdf.results.dat'.format(ch=ch, pset=pset), 'r')
    
        print 'at file', tmp_infile
    
    
        lines = [line for line in tmp_infile.readlines() if len(line.strip()) ]
        
        ind_rap = [il for il,l in enumerate(lines) if 'Z/W rapidity' in l][0]
        ind_a0  = [il for il,l in enumerate(lines) if 'A_0 numer'    in l][0]
        ind_a1  = [il for il,l in enumerate(lines) if 'A_1 numer'    in l][0]
        ind_a4  = [il for il,l in enumerate(lines) if 'A_4 numer'    in l][0]
        ind_de  = [il for il,l in enumerate(lines) if 'A_i denom'    in l][0]
        
        lines_rap = addPdfIndex([i.split() for i in lines[ind_rap+1 :ind_a0]])
        # lines_a0  = addPdfIndex([i.split() for i in lines[ind_a0+1  :ind_a1]])
        # lines_a4  = addPdfIndex([i.split() for i in lines[ind_a4+1  :ind_de]])
        # lines_de  = addPdfIndex([i.split() for i in lines[ind_de+1  :      ]])
    
        # lines_a0  = divideByDenom(lines_a0, lines_de)
        # lines_a4  = divideByDenom(lines_a4, lines_de)
    
        ### print 'this is lines rap', lines_rap
    
        tmp_infile.close()

        tmp_infile_cen = open(basedir+'/my_w{ch}_18012_{pset}/NNLO.results.dat'.format(ch=ch, pset=pset), 'r')
    
        print 'at file for central values ', tmp_infile_cen
    
        lines = [line for line in tmp_infile_cen.readlines() if len(line.strip()) ]
        
        ind_rap = [il for il,l in enumerate(lines) if 'Z/W rapidity' in l][0]
        ind_a0  = [il for il,l in enumerate(lines) if 'A_0 numer'    in l][0]
        
        lines_rap_cen = addPdfIndex([i.split() for i in lines[ind_rap+1 :ind_a0]], offset=0)
        #print lines_rap_cen

        ##npdf = max(i[0] for i in lines_rap)
    
        ## for ipdf in range(npdf+1):
    
        ##    graph_vals = [lines_rap_cen] + [i for i in lines_rap if i[0] is ipdf]

        new_lines_rap = []

        graph = getGraph(pdfset, lines_rap_cen, lines_rap)
        
        for central in lines_rap_cen:
            new_lines_rap.append(central)
            tmp_rap = central[1]


            new_lines_rap += [i for i in lines_rap if i[1] == tmp_rap]
            

        #vals[ch] = new_lines_rap
        vals[ch] = graph


    return vals
        
    #    tmp_iy = -1
    #    ntmp = 0
    #    for iline, line in enumerate(lines_a0):
    #        if line[0] != tmp_iy:
    #            tmp_iy = line[0]
    #            ntmp = 0
    #        lines_a0[iline].insert(0, ntmp)
    #        ntmp += 1
        
        #lines_rap = returnArray(lines[ind_rap:ind_a0], ipdf)
        #lines_a0  = returnArray(lines[ind_a0 :ind_a1], ipdf)
        #lines_a4  = returnArray(lines[ind_a4 :ind_de], ipdf)
        #lines_de  = returnArray(lines[ind_de :      ], ipdf)
        #
        ### print lines_a0
        ### print lines_a4
        #
        #lines_l = getHelicity(ch, 'l', lines_a0, lines_a4)
        #lines_r = getHelicity(ch, 'r', lines_a0, lines_a4)
        #
        #basename = 'graph_W{ch}_{pset}_'.format(ch='plus' if ch == 'p' else 'minus', pset=pset)
        #
        #tmp_graph_tot = ROOT.TGraphAsymmErrors(len(lines_rap[0]), lines_rap[0], lines_rap[1], lines_rap[2], lines_rap[3], lines_rap[6], lines_rap[6])
        #tmp_graph_a0  = ROOT.TGraphAsymmErrors(len(lines_a0 [0]), lines_a0 [0], lines_a0 [1], lines_a0 [2], lines_a0 [3], lines_a0 [6], lines_a0 [6])
        #tmp_graph_a4  = ROOT.TGraphAsymmErrors(len(lines_a4 [0]), lines_a4 [0], lines_a4 [1], lines_a4 [2], lines_a4 [3], lines_a4 [6], lines_a4 [6])
        #tmp_graph_fl  = ROOT.TGraphAsymmErrors(len(lines_l [0]), lines_l [0], lines_l [1], lines_l [2], lines_l [3], lines_l [6], lines_l [6])
        #tmp_graph_fr  = ROOT.TGraphAsymmErrors(len(lines_r [0]), lines_r [0], lines_r [1], lines_r [2], lines_r [3], lines_r [6], lines_r [6])
        #
        #tmp_graph_tot.SetName(basename+'total'); tmp_graph_tot.SetTitle(basename+'total')
        #tmp_graph_a0 .SetName(basename+'a0'   ); tmp_graph_a0 .SetTitle(basename+'a0'   )
        #tmp_graph_a4 .SetName(basename+'a4'   ); tmp_graph_a4 .SetTitle(basename+'a4'   )
        #tmp_graph_fl .SetName(basename+'fl'   ); tmp_graph_fl .SetTitle(basename+'fl'   )
        #tmp_graph_fr .SetName(basename+'fr'   ); tmp_graph_fr .SetTitle(basename+'fr'   )
        #
        #tmp_graph_tot . Write()
        #tmp_graph_a0  . Write()
        #tmp_graph_a4  . Write()
        #tmp_graph_fl  . Write()
        #tmp_graph_fr  . Write()
    
    #outfile.Close()
    
        
