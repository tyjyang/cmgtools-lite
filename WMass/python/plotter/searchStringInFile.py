import ROOT, os, sys, re, array, math
import time

# originally designed to read job_XXX.sh files looking for keywords (like a process's name)
#Z_mu_minus

indir = "cards/diffXsec_mu_2019_05_09_recoEta0p1_recoPt1_genEta0p2from1p3_last2p1to2p4_genPt2/jobs/" 
fileNameMatch = "job_.*"
keyword = " Z_mu_plus "
matchingFilesId = []

files = [f for f in os.listdir(indir) if re.match(fileNameMatch,f) and "bkg" not in f]
for f in files:
    #print f
    fnum = int((f.split("_")[1]).split(".sh")[0])
    ifound = 1
    with open(indir+f) as ff:
        lines = ff.readlines()
        for line in lines:
            if keyword in line:
                print "{n}: found '{k}' in file {ff}".format(n=ifound, k=keyword, ff=ff.name)
                if ifound == 1: matchingFilesId.append(fnum)
                ifound += 1

print "-"*30
if len(matchingFilesId):
    print "Found {k} in {n} files with these ids".format(k=keyword,n=len(matchingFilesId))
    print  matchingFilesId
else:
    print "No match found!"
print "-"*30
