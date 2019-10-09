import ROOT, itertools, math
from array import array


def returnArray(lines):

    retx   = array('d', [])
    retval = array('d', [])
    retnum = array('d', [])
    retpdf = array('d', [])
    rettot = array('d', [])

    retxup = array('d', [])
    retxdn = array('d', [])

    for line in lines:
        tmp_line = line.strip()
        if not tmp_line[0].isdigit(): continue

        retx   .append( float(tmp_line.split()[0]) )
        retval .append( float(tmp_line.split()[1]) )
        retnum .append( float(tmp_line.split()[2]) )
        retpdf .append( float(tmp_line.split()[3]) )
        rettot .append( math.sqrt(retnum[-1]**2 + retpdf[-1]**2) )
        retxup.append(0.)
        retxdn.append(0.)

    

    return retx, retval, retxdn, retxup, retnum, retpdf, rettot
        


basedir = '/afs/cern.ch/work/a/arapyan/public/fewz_smp18012/'

outfile_name = '/afs/cern.ch/work/m/mdunser/public/cmssw/w-helicity-13TeV/CMSSW_8_0_25/src/CMGTools/WMass/data/theory/theory_cross_sections_FEWZ.root'
outfile = ROOT.TFile(outfile_name, 'RECREATE')

for ch, o in itertools.product('pm', ['nnlo', 'nlo']):


    tmp_infile = open(basedir+'/w{ch}_{o}.dat'.format(ch=ch, o=o), 'r')


    lines = [line for line in tmp_infile.readlines() if len(line.strip()) ]

    ind_rap = [il for il,l in enumerate(lines) if 'Z/W rapidity' in l][0]
    ind_a0  = [il for il,l in enumerate(lines) if 'A0 vs'        in l][0]
    ind_a1  = [il for il,l in enumerate(lines) if 'A1 vs'        in l][0]
    ind_a4  = [il for il,l in enumerate(lines) if 'A4 vs'        in l][0]


    lines_rap = returnArray(lines[ind_rap:ind_a0])
    lines_a0  = returnArray(lines[ind_a0 :ind_a1])
    lines_a4  = returnArray(lines[ind_a4 :      ])

    basename = 'graph_W{ch}_{o}_'.format(ch='plus' if ch == 'p' else 'minus', o=o)

    tmp_graph_tot = ROOT.TGraphAsymmErrors(len(lines_rap[0]), lines_rap[0], lines_rap[1], lines_rap[2], lines_rap[3], lines_rap[6], lines_rap[6])
    tmp_graph_a0  = ROOT.TGraphAsymmErrors(len(lines_a0 [0]), lines_a0 [0], lines_a0 [1], lines_a0 [2], lines_a0 [3], lines_a0 [6], lines_a0 [6])
    tmp_graph_a4  = ROOT.TGraphAsymmErrors(len(lines_a4 [0]), lines_a4 [0], lines_a4 [1], lines_a4 [2], lines_a4 [3], lines_a4 [6], lines_a4 [6])
    tmp_graph_tot.SetName(basename+'total'); tmp_graph_tot.SetTitle(basename+'total')
    tmp_graph_a0 .SetName(basename+'a0'   ); tmp_graph_a0 .SetTitle(basename+'a0'   )
    tmp_graph_a4 .SetName(basename+'a4'   ); tmp_graph_a4 .SetTitle(basename+'a4'   )
    tmp_graph_tot . Write()
    tmp_graph_a0  . Write()
    tmp_graph_a4  . Write()

outfile.Close()

    
