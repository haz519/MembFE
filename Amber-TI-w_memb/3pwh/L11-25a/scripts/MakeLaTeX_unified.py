#!/usr/bin/env python

if __name__ == "__main__":

    import matplotlib
    matplotlib.use("Agg")

    from tianalysis import DataLoc as dloc
    from tianalysis import UsePkaUnits
    from tianalysis import make_latex_document
    from collections import defaultdict as ddict
    import os
    import glob
    
    #UsePkaUnits()
    
    tequil=200
    tmax=1.e+10
    odir="latex"
    D = ddict( lambda: ddict( lambda: ddict( list ) ) )
    for i in ["LIG"]:
        unif=[]
        for itrial in range(9):
            trial="t%02i"%(itrial+1)
            d = "ligand/concerted/%s/results/"%(trial)
            if os.path.exists(d):
                if len(glob.glob(d+"dvdl_*.dat")) > 0:
                    unif.append( dloc(trial,d,"",tequil,tmax) )
                    
        if len(unif) > 0:
            D[i]["ref"]["unif"]  = unif

        unif=[]
        for itrial in range(9):
            trial="t%02i"%(itrial+1)
            d = "complex/concerted/%s/results/"%(trial)
            if os.path.exists(d):
                if len(glob.glob(d+"dvdl_*.dat")) > 0:
                    unif.append( dloc(trial,d,"",tequil,tmax) )
                    
        if len(unif) > 0:
            D[i]["bio"]["unif"]  = unif
    make_latex_document( odir, D, methods=["TI","TI3","BAR","MBAR"] )
    

