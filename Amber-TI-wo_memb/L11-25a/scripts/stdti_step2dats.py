#!/usr/bin/env python
import sys,os

def extract_traditional_ti( fname, write=False, manual_lambda=False ):
    import os
    from collections import defaultdict as ddict

    fh = open(fname,"r")
    if not fh:
        raise Exception("Could not open %s\n"%(fname))



    numexchg=0
    nstlim=None
    ntpr=None
    dt=None
    irest=0
    for line in fh:
        cmdstr,sepstr,comstr = line.partition("!")
        if "ntpr" in cmdstr:
            cols = cmdstr.replace("=","").replace(",","").strip().split()
            for icol in range(len(cols)-1):
                if cols[icol] == "ntpr":
                    ntpr = int( cols[icol+1] )
                    break
        if "dt" in cmdstr:
            cols = cmdstr.replace("=","").replace(",","").strip().split()
            for icol in range(len(cols)-1):
                if cols[icol] == "dt":
                    dt = float( cols[icol+1] )
                    break
        if "numexchg" in cmdstr:
            cols = cmdstr.replace("=","").replace(",","").strip().split()
            for icol in range(len(cols)-1):
                if cols[icol] == "numexchg":
                    numexchg = int( cols[icol+1] )
                    break
        if "nstlim" in cmdstr:
            cols = cmdstr.replace("="," ").replace(",","").strip().split()
            for icol in range(len(cols)-1):
                if cols[icol] == "nstlim":
                    nstlim = int( cols[icol+1] )
                    break
        if "irest" in cmdstr:
            cols = cmdstr.replace("="," ").replace(",","").strip().split()
            for icol in range(len(cols)-1):
                if cols[icol] == "irest":
                    irest = int( cols[icol+1] )
                    break

    if ntpr is None:
        raise Exception("Could not determine ntpr from %s"%(fname))

    if dt is None:
        raise Exception("Could not determine dt from %s"%(fname))

    if nstlim is None:
        raise Exception("Could not determine nstlim from %s"%(fname))

    if numexchg < 1:
        numexchg = 1

    dt = dt
    nstep_per_sim = nstlim * numexchg
    nframe_per_sim = nstep_per_sim / ntpr

    if nstep_per_sim % ntpr != 0:
        print ("num md steps per simulation is not a multiple of ntpr. Unclear how the simulation time works")

    t_per_frame = dt * ntpr
    t_per_sim = t_per_frame * nframe_per_sim


    fh = open(fname,"r")

    
    efeps = []
    dvdls = []
    efep = []
    reading_region_1 = False

    lam = None
    nlam = 0
    
    for line in fh:
        if "A V E R A G E S" in line:
            break
        if "MBAR Energy analysis:" in line:
            efep = []
        if "clambda" in line:
            if lam is None:
                cols = line.replace("="," ").replace(","," ").split()
                for i in range(len(cols)):
                    if cols[i] == "clambda":
                        lam = float(cols[i+1])
                        break
        elif "Energy at " in line:
            #print (line)
            val = line.strip().split()[-1]
            if "****" in val:
                val = 10000
                #if len(efep) > 0:
                #   if efep[-1] < 0:
                #       val = -val
            else:
                val = float(val)
            efep.append( val )
        elif "TI region  2" in line:
            reading_region_1 = True
            dvdl = 0
        elif "| TI region  1" in line:
            #print (line)
            reading_region_1 = False
            #print (dvdl)
            dvdls.append( dvdl )
            if len( efep ) > 0:
                efeps.append( efep )
                nlam = len(efep)
        elif "TI region " in line:
            reading_region_1 = False

        if "DV/DL  =" in line and reading_region_1:
            #print (line)
            cols = line.strip().split()
            #print (cols)
            dvdl = float( cols[-1] )
            #dvdls.append( float( cols[-1] ) )
            #if len( efep ) > 0:
            #    efeps.append( efep )
            #    nlam = len(efep)
    if write:
        if manual_lambda:
            lams = manual_lambda
        else:
            lams = [ float(i) / ( nlam-1. ) for i in range(nlam) ]
        print(("MANUAL LAMBDA: ", manual_lambda))
        print(("LAMS: ", lams))
        for l in lams:
            if abs(l-lam) < 0.001:
                lam = l
                break
        head, tail = os.path.split(fname)
        dvdl_fname = os.path.join( head, "dvdl_%.8f.dat"%( lam ) )

        if irest == 0:
           dvdls=dvdls[1:]

        fh = open(dvdl_fname,"w")
        for i in range(len(dvdls)):
            fh.write("%.4f %18.6f\n"%((i+1)*t_per_frame,dvdls[i]))
        fh.close()
        for ilam,plam in enumerate(lams):
            efep_fname = os.path.join( head, "efep_%.8f_%.8f.dat"%( lam, plam ) )
            fh = open(efep_fname,"w")
            for i in range(len(efeps)):
                fh.write("%.4f %18.6f\n"%((i+1)*t_per_frame,efeps[i][ilam]))
            fh.close()

    return dvdls,efeps


# be careful to set lambdas!
#manual_lambda = False
manual_lambda = [0.0000, 0.0479, 0.1151, 0.2063, 0.3161, 0.4374, 0.5626, 0.6839, 0.7937, 0.8850, 0.9521, 1.0000]
for arg in sys.argv[1:]:
    if os.path.isfile( arg ):
        if ".out" in arg:
            extract_traditional_ti( arg, write=True, manual_lambda=manual_lambda )
        else:
            print(("File does not end in .out: %s"%(arg)))
    else:
        print(("File not found: %s"%(arg)))

