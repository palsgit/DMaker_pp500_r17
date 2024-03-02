import ROOT
import os, sys

print ("Merging %s" % sys.argv[1])

print ("Max tree size",ROOT.TTree.GetMaxTreeSize())
ROOT.TTree.SetMaxTreeSize(2000000000000) # 2 Tb
print ("Updated tree size",ROOT.TTree.GetMaxTreeSize())

rm = ROOT.TFileMerger(False)
rm.SetFastMethod(True)


path = './'
file_output = '%s.root' % sys.argv[1]
file_list = []
for path, dirs, files in os.walk(path):
  for filename in files:
    rfile = ROOT.TFile.Open(filename)
    if not rfile: continue
    del rfile
    ##if not rm.AddFile(filename): continue
    if ('picoD0AnaMaker.root') in filename: file_list.append(path+filename)

print ("Input file list:",file_list)
print ("Output file:",file_output)


##rmfinal = ROOT.TFileMerger(False)
for F in file_list:
    
    print ("Adding ->",F)
    rm.AddFile(F)

rm.OutputFile(file_output)
rm.Merge()
