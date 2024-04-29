#!/usr/bin/env python

import numpy as np
import pymbar
from collections import defaultdict as ddict



def GetAvgAndStdErr(vals,errs):
    import numpy as np
    avg=0
    err=0
    if len(vals) > 0:
        avg = np.mean( vals )
        if len(vals) > 1:
            err = np.var( vals, ddof=1 )
        for e in errs:
            err += e*e
        err = np.sqrt( err / len(vals) )
    return avg,err
    
def GetSumAndStdErr(vals,errs):
    import numpy as np
    return np.sum(vals),np.sqrt( np.sum( [ e*e for e in errs ] ) )

def splitall(path):
    import os
    allparts = []
    while True:
        parts = os.path.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        elif parts[1] == path: # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        else:
            path = parts[0]
            allparts.insert(0, parts[1])
    if allparts[0][0] == "/":
        allparts[0] = allparts[0][1:]
    return allparts


def get_minmax(ymin,ymax,y,dy=None):
    if dy is None:
        dy = [ 0 ]*len(y)
    for i in range(len(y)):
        u = y[i]+1.5*dy[i]
        p = (u-ymin)/(ymax-ymin)
        if p > 0.9:
            ymax = (u-ymin)/0.9 + ymin
        u = y[i]-1.5*dy[i]
        p = (u-ymin)/(ymax-ymin)
        if p < 0.12:
            ymin = (0.12*ymax-u)/(0.12-1.)
    return ymin,ymax


def get_energy_minmax_from_trial(ymin,ymax,trial,method):
    n = len(trial.times)
    
    y   = [ trial.forward[i].ene[method][0] for i in range(n) ]
    dy  = [ 1.96 * trial.forward[i].ene[method][1] for i in range(n) ]
    ymin,ymax = get_minmax(ymin,ymax,y,dy)
    
    y   = [ trial.reverse[i].ene[method][0] for i in range(n) ]
    dy  = [ 1.96 * trial.reverse[i].ene[method][1] for i in range(n) ]
    ymin,ymax = get_minmax(ymin,ymax,y,dy)

    y   = [ trial.segment[i].ene[method][0] for i in range(n) ]
    ymin,ymax = get_minmax(ymin,ymax,y,dy)
    return ymin,ymax

def get_dvdl_minmax_from_trial(ymin,ymax,trial,ilams):
    n = len(trial.times)
    for ilam in ilams:
        n   = len(trial.times)
        y   = [ trial.forward[i].dvdl[ilam][0] for i in range(n) ]
        dy  = [ 1.96 * trial.forward[i].dvdl[ilam][1] for i in range(n) ]
        ymax = max( ymax, max( [ y[i]+1.5*dy[i] for i in range(n) ] ) )
        ymin = min( ymin, min( [ y[i]-1.5*dy[i] for i in range(n) ] ) )
        y   = [ trial.reverse[i].dvdl[ilam][0] for i in range(n) ]
        dy  = [ 1.96 * trial.reverse[i].dvdl[ilam][1] for i in range(n) ]      
        ymax = max( ymax, max( [ y[i]+1.5*dy[i] for i in range(n) ] ) )
        ymin = min( ymin, min( [ y[i]-1.5*dy[i] for i in range(n) ] ) )
        y   = [ trial.segment[i].dvdl[ilam][0] for i in range(n) ]
        ymax = max( ymax, max( [ y[i] for i in range(n) ] ) )
        ymin = min( ymin, min( [ y[i] for i in range(n) ] ) )
        
    for ilam in ilams:
        y   = [ trial.forward[i].dvdl[ilam][0] for i in range(n) ]
        dy  = [ 1.96 * trial.forward[i].dvdl[ilam][1] for i in range(n) ]
        ymin,ymax = get_minmax(ymin,ymax,y,dy)
    
        y   = [ trial.reverse[i].dvdl[ilam][0] for i in range(n) ]
        dy  = [ 1.96 * trial.reverse[i].dvdl[ilam][1] for i in range(n) ]
        ymin,ymax = get_minmax(ymin,ymax,y,dy)

        y   = [ trial.segment[i].dvdl[ilam][0] for i in range(n) ]
        ymin,ymax = get_minmax(ymin,ymax,y,dy)
    return ymin,ymax


def plot_energy_timeseries(ax,trial,method,ymin=1.e+10,ymax=-1.e+10):
    n   = len(trial.times)
    x   = [ time/1000. for time in trial.times ]

    if method not in trial.forward[0].ene:
        for alts in [ "BAR","IEXP","DEXP" ]:
            if alts in trial.forward[0].ene:
                method = alts
                break


            
    y   = [ trial.forward[i].ene[method][0] for i in range(n) ]
    #print (y)
    #print (x)
    dy  = [ 1.96 * trial.forward[i].ene[method][1] for i in range(n) ]
    ymax = max( ymax, max( [ y[i]+1.5*dy[i] for i in range(n) ] ) )
    ymin = min( ymin, min( [ y[i]-1.5*dy[i] for i in range(n) ] ) )
    ax.errorbar(x,y,yerr=dy,capsize=1,c='k',linestyle='solid',linewidth=1,marker='o',ms=2)

    y   = [ trial.reverse[i].ene[method][0] for i in range(n) ]
    dy  = [ 1.96 * trial.reverse[i].ene[method][1] for i in range(n) ]
    ymax = max( ymax, max( [ y[i]+1.5*dy[i] for i in range(n) ] ) )
    ymin = min( ymin, min( [ y[i]-1.5*dy[i] for i in range(n) ] ) )
    ax.errorbar(x,y,yerr=dy,capsize=1,c='r',linestyle='solid',linewidth=1,marker='o',ms=2)

    y   = [ trial.segment[i].ene[method][0] for i in range(n) ]
    dy  = [ 1.96 * trial.segment[i].ene[method][1] for i in range(n) ]
    ymax = max( ymax, max( [ y[i] for i in range(n) ] ) )
    ymin = min( ymin, min( [ y[i] for i in range(n) ] ) )
    ax.errorbar(x,y,capsize=0,c='g',linestyle='dotted',linewidth=1,marker='x',ms=2)
    
    ymin,ymax = get_energy_minmax_from_trial(ymin,ymax,trial,method)
    ax.set_xlim(0,x[-1]+x[0])
    delta=0
    if ymax == ymin:
        delta = 0.1
    ax.set_ylim(ymin-delta,ymax+delta)
    if UsesPkaUnits():
        ax.set_ylabel(r"$\Delta$pK${}_{\mathregular{a}}$")
    else:
        ax.set_ylabel(r"$\Delta$G (kcal/mol)")
    ax.set_xlabel("Time (ns)")
    ax.grid(linestyle='dotted',linewidth=0.5)
    
    return ymin,ymax
            

class naturalcubicspline:

    def __init__(self, x):

        # define some space
        L = len(x)
        H = np.zeros([L,L],float)
        M = np.zeros([L,L],float)
        BW = np.zeros([L,L],float)
        AW = np.zeros([L,L],float)
        DW = np.zeros([L,L],float)

        x = np.array(x)
        h = x[1:L]-x[0:L-1]
        ih = 1.0/h
        
        # define the H and M matrix, from p. 371 "applied numerical methods with matlab, Chapra"
        H[0,0] = 1
        H[L-1,L-1] = 1
        for i in range(1,L-1):
            H[i,i] = 2*(h[i-1]+h[i])
            H[i,i-1] = h[i-1]
            H[i,i+1] = h[i]
            
            M[i,i] = -3*(ih[i-1]+ih[i])
            M[i,i-1] = 3*(ih[i-1])
            M[i,i+1] = 3*(ih[i])
            
        CW = np.dot(np.linalg.inv(H),M)
        # this is the matrix translating c to weights in f.
        # each row corresponds to the weights for each c.

        # from CW, define the other coefficient matrices
        for i in range(0,L-1):
            BW[i,:]    = -(h[i]/3)*(2*CW[i,:]+CW[i+1,:])
            BW[i,i]   += -ih[i]
            BW[i,i+1] += ih[i]
            DW[i,:]    = (ih[i]/3)*(CW[i+1,:]-CW[i,:])
            AW[i,i]    = 1

        # Make copies of the arrays we'll be using in the future.
        self.x  = x.copy()
        self.AW = AW.copy()
        self.BW = BW.copy()
        self.CW = CW.copy()
        self.DW = DW.copy()

    def interpolate(self,y,xnew):
        if len(self.x) != len(y):
            parser.error("naturalcubicspline interpolate size mismatch")
        a = np.dot(self.AW,y)
        b = np.dot(self.BW,y)
        c = np.dot(self.CW,y)
        d = np.dot(self.DW,y)
        
        N = len(xnew)
        ynew = np.zeros([N],float)
        for i in range(N):
            # Find the index of 'xnew[i]' it would have in 'self.x'.
            j = np.searchsorted(self.x, xnew[i]) - 1
            lamw = xnew[i] - self.x[j]
            ynew[i] = d[j]*lamw**3 + c[j]*lamw**2 + b[j]*lamw + a[j]
        # Preserve the terminal points.
        ynew[0] = y[0]
        ynew[-1] = y[-1]
        return ynew

    def integrate(self,y,dy):
        N = len(self.x)
        if N != len(y):
            parser.error("naturalcubicspline integrate size mismatch")
        a = np.dot(self.AW,y)
        b = np.dot(self.BW,y)
        c = np.dot(self.CW,y)
        d = np.dot(self.DW,y)

        avec = np.zeros( [N], float )
        bvec = np.zeros( [N], float )
        cvec = np.zeros( [N], float )
        dvec = np.zeros( [N], float )
        for j in range(N-1):
            delx = self.x[j+1]-self.x[j]
            c4=(delx**4)/4.
            c3=(delx**3)/3.
            c2=(delx**2)/2.
            c1=delx
            avec[j] = c1
            bvec[j] = c2
            cvec[j] = c3
            dvec[j] = c4
        avec = np.dot(avec,self.AW)
        bvec = np.dot(bvec,self.BW)
        cvec = np.dot(cvec,self.CW)
        dvec = np.dot(dvec,self.DW)
        wvec = avec+bvec+cvec+dvec
        inty   = np.dot( wvec, y )
        intdy2 = np.dot( wvec**2, dy**2 )
        
        GvsL = np.zeros( [N], float )
        DvsL= np.zeros( [N], float )
        for j in range(N-1):
            delx = self.x[j+1]-self.x[j]
            c4=(delx**4)/4.
            c3=(delx**3)/3.
            c2=(delx**2)/2.
            c1=delx
            avec = np.zeros( [N], float )
            bvec = np.zeros( [N], float )
            cvec = np.zeros( [N], float )
            dvec = np.zeros( [N], float )
            avec[j] = c1
            bvec[j] = c2
            cvec[j] = c3
            dvec[j] = c4
            avec = np.dot(avec,self.AW)
            bvec = np.dot(bvec,self.BW)
            cvec = np.dot(cvec,self.CW)
            dvec = np.dot(dvec,self.DW)
            wvec = avec+bvec+cvec+dvec
            GvsL[j+1] = GvsL[j] + np.dot( wvec, y )
            DvsL[j+1] = DvsL[j] + np.dot( wvec**2, dy**2 )
        for j in range(N-1):
            DvsL[j+1] = np.sqrt( DvsL[j+1] )
        return inty,np.sqrt(intdy2),GvsL,DvsL


def CalcTI_spline( lams, dvdl, dvdl_stderrs ):
    return naturalcubicspline( lams ).integrate( dvdl, dvdl_stderrs )


def CalcTI_trapez( lams, dvdl, dvdl_stderrs ):
    f=0
    df2=0.
    n=len(lams)
    wts = np.zeros( [n], float )
    GvsL = np.zeros( [n], float )
    DvsL = np.zeros( [n], float )
    for k in range(n-1):
        c         = 0.5*(lams[k+1]-lams[k])
        wts[k+1] += c
        wts[k]   += c
        GvsL[k+1] = GvsL[k] + c * ( dvdl[k] + dvdl[k+1] )
        DvsL[k+1] = DvsL[k] + c*c*( dvdl_stderrs[k]**2 + dvdl_stderrs[k+1]**2 )
    for k in range(n-1):
        DvsL[k+1] = np.sqrt( DvsL[k+1] ) 
    for k in range(n):
        c    = wts[k]
        f   += c   * dvdl[k]
        df2 += c*c * dvdl_stderrs[k]**2

    return f,np.sqrt(df2),GvsL,DvsL






def plotOverlapMatrix(O,K,fname):
    import numpy
    from matplotlib import pyplot as pl
    """Plots the probability of observing a sample from state i (row) in state j (column).
    For convenience, the neigboring state cells are fringed in bold."""
    max_prob = O.max()
    fig = pl.figure(figsize=(K/2.,K/2.))
    fig.add_subplot(111, frameon=False, xticks=[], yticks=[])
    
    for i in range(K):
        if i!=0:
            pl.axvline(x=i, ls='-', lw=0.5, color='k', alpha=0.25)
            pl.axhline(y=i, ls='-', lw=0.5, color='k', alpha=0.25)
        for j in range(K):
            if O[j,i] < 0.005:
                ii = ''
            elif O[j,i] > 0.995:
                ii = '1.00'
            else:
                ii = ("%.2f" % O[j,i])[1:]
            alf = O[j,i]/max_prob
            #pl.fill_between([i,i+1], [K-j,K-j], [K-(j+1),K-(j+1)], color='k', alpha=alf)
            #pl.annotate(ii, xy=(i,j), xytext=(i+0.5,K-(j+0.5)), size=8, textcoords='data', va='center', ha='center', color=('k' if alf < 0.5 else 'w'))
            pl.annotate(ii, xy=(i,j), xytext=(i+0.5,K-(j+0.5)), size=8, textcoords='data', va='center', ha='center', color=('k' if alf < 0.5 else 'k'))
            

            
    ks = list(range(K))
    for i in range(K):
        pl.annotate(ks[i], xy=(i+0.5, 1), xytext=(i+0.5, K+0.5), size=10, textcoords=('data', 'data'), va='center', ha='center', color='k')
        pl.annotate(ks[i], xy=(-0.5, K-(j+0.5)), xytext=(-0.5, K-(i+0.5)), size=10, textcoords=('data', 'data'), va='center', ha='center', color='k')
    pl.annotate('$\lambda$', xy=(-0.5, K-(j+0.5)), xytext=(-0.5, K+0.5), size=10, textcoords=('data', 'data'), va='center', ha='center', color='k')
    pl.plot([0,K], [0,0], 'k-', lw=4.0, solid_capstyle='butt')
    pl.plot([K,K], [0,K], 'k-', lw=4.0, solid_capstyle='butt')
    pl.plot([0,0], [0,K], 'k-', lw=2.0, solid_capstyle='butt')
    pl.plot([0,K], [K,K], 'k-', lw=2.0, solid_capstyle='butt')
    
    cx = sorted(2*list(range(K+1)))
    cy = sorted(2*list(range(K+1)), reverse=True)
    pl.plot(cx[2:-1], cy[1:-2], 'k-', lw=2.0)
    pl.plot(numpy.array(cx[2:-3])+1, cy[1:-4], 'k-', lw=2.0)
    pl.plot(cx[1:-2], numpy.array(cy[:-3])-1, 'k-', lw=2.0)
    pl.plot(cx[1:-4], numpy.array(cy[:-5])-2, 'k-', lw=2.0)
    
    pl.xlim(-1, K)
    pl.ylim(0, K+1)
    pl.savefig(fname, bbox_inches='tight', pad_inches=0.0)
    pl.close(fig)
    return









class HamiltoniansAtThisFrame(object):
    
    def __init__(self,num_states):
        self.num_states=num_states
        self.ene  = [ None ] * self.num_states
        self.dvdl = None


        
class FramesInThisTraj(object):
    
    def __init__(self,num_states):
        self.num_states=num_states
        self.timesteps = ddict( lambda: HamiltoniansAtThisFrame(self.num_states) )


    def uncorrelated_frames(self,tequil=0,tmax=1.e+10,state=None):
        #
        # returns a new FramesInThisTraj object containing the
        # UNCORRELATED data within tequil <= time <= tmax
        # The autocorrelation function uses all time >= tequil
        #
        times = [ k for k in self.timesteps ]
        times.sort()
        if len(times) == 0:
            return FramesInThisTraj(self.num_states)
        if tequil >= times[-1]:
            tequil = times[-1]/2.
        if tmax < tequil:
            tmax = tequil*2
        idxs = self.get_uncorrelated_indexes(tequil=tequil,state=state)
        #print ("ui",len(idxs))
        a = FramesInThisTraj(self.num_states)
        for i,time in enumerate(sorted(self.timesteps)):
            if time >= tequil and i in idxs and time <= tmax:
                a.timesteps[time] = self.timesteps[time]
        if len(a.timesteps) == 0:
            for i,time in enumerate(sorted(self.timesteps)):
                if time >= tequil and i in idxs:
                    a.timesteps[time] = self.timesteps[time]
                    break
        return a

    
    def extract_frames(self,tequil=0,tmax=1.e+10):
        #
        # returns a new FramesInThisTraj object containing all
        # data within tequil <= time <= tmax
        #
        times = [ k for k in self.timesteps ]
        times.sort()
        if len(times) == 0:
            return FramesInThisTraj(self.num_states)
        if tequil >= times[-1]:
            tequil = times[-1]/2.
        if tmax < tequil:
            tmax = tequil*2
            
        a = FramesInThisTraj(self.num_states)
        for i,time in enumerate(sorted(self.timesteps)):
            if time >= tequil and time <= tmax:
                a.timesteps[time] = self.timesteps[time]
        if len(a.timesteps) == 0:
            for i,time in enumerate(sorted(self.timesteps)):
                if time >= tequil:
                    a.timesteps[time] = self.timesteps[time]
                    break
        return a

    
    def get_uncorrelated_indexes(self,tequil=0,state=None):
        #
        # returns a list of indexs corresponding to the
        # uncorrelated times with time >= tequil
        #        
        times = [ k for k in self.timesteps ]
        times.sort()
        if tequil >= times[-1]:
            tequil = times[-1]/2.
            
        for itime,t in enumerate(sorted(self.timesteps)):
            istart = itime
            if t >= tequil:
                break                
        g,dat = self.get_statistical_inefficiency(tequil,state=state,with_data=True)
        #print (g,len(dat),len(self.timesteps),tequil,times[-1])
        idxs = np.array(pymbar.timeseries.subsampleCorrelatedData(dat,g=g))
        #print (len(idxs))
        return istart + idxs
    

    def has_dvdl(self):
        has=True
        if None in [ self.timesteps[t].dvdl for t in self.timesteps ]:
            has=False
        return has

    
    def get_statistical_inefficiency(self,tequil=0,tmax=1.e+10,state=None,with_data=False):
        #
        # if state=None, then the dvdl data is used to compute the inefficiency
        # if state=istate, then the forward (or reverse) dE's are used
        #     
        times = [ k for k in self.timesteps ]
        times.sort()
        if tequil >= times[-1]:
            tequil = times[-1]/2.
        if tmax < tequil:
            tmax = tequil*2
            
        if state is None:
            if None in [  self.timesteps[t].dvdl \
                          for t in sorted(self.timesteps) if t >= tequil and t <= tmax ]:
                raise Exception("Can't use dvdl for autocorrelation times because not all dvdl data is available; rerun using the --no-ti option")
                
        if state is None:
            dat = [ self.timesteps[t].dvdl \
                    for t in sorted(self.timesteps) if t >= tequil and t <= tmax ]
            #print (len(dat))
        else:
            if state < self.num_states - 1:
                dat = [ self.timesteps[t].ene[state+1]-self.timesteps[t].ene[state] \
                        for t in sorted(self.timesteps) if t >= tequil and t <= tmax ]
            else:
                dat = [ self.timesteps[t].ene[state]-self.timesteps[t].ene[state-1] \
                        for t in sorted(self.timesteps) if t >= tequil and t <= tmax ]
        dat = np.array(dat)
        if with_data:
            try:
                return pymbar.timeseries.statisticalInefficiency(dat),dat
            except Exception as e:
                print (e)
                return 1,dat
        else:
            try:
                return pymbar.timeseries.statisticalInefficiency(dat)
            except Exception as e:
                print (e)
                return 1


    
    def get_correlation_time(self,tequil=0,tmax=1.e+10,state=None):
        #
        # reports the statistical inefficiency as an autocorrelation time in ps
        #
        if len(self.timesteps) > 1:
            s = self.get_statistical_inefficiency(tequil=tequil,tmax=tmax,state=state)
            times = [ t for t in self.timesteps ]
            times.sort()
            dt = times[1]-times[0]
            tau = 0.5*(s-1.)
            return dt * tau
        else:
            return 0

    

def UsePkaUnits(pka=True):
    AlchemicalTransform.global_pka_units = pka

def UsesPkaUnits():
    return AlchemicalTransform.global_pka_units

    
class AlchemicalTransform(object):

    global_pka_units = False
    global_pka_per_kcal = 1. / 1.364
    
    def __init__(self,lams,directory=None,label=None):
        self.lams = list(set(lams))
        self.lams.sort()
        self.nlam = len(self.lams)
        self.clear()
        self.label = "unlabeled"
        
        self.mbar=True
        self.bar=True
        self.ti=True
        self.ti3=True
        self.iexp=True
        self.dexp=True

        if directory is not None:
            self.read_from_files(directory,label=label)
        
        
    def clear(self):
        self.state = []
        for lam in self.lams:
            self.state.append( FramesInThisTraj(self.nlam) )

            
    def get_methods(self):
        methods=[]
        if self.ti:
            methods.append("TI")
        if self.ti3:
            methods.append("TI3")
        if self.bar:
            methods.append("BAR")
        if self.mbar:
            methods.append("MBAR")
        if self.iexp:
            methods.append("IEXP")
        if self.dexp:
            methods.append("DEXP")
        return methods
    
            
    def read_from_files(self,directory,label=None):
        import os
        self.clear()
        if len(label) == 0:
            label = None
        self.label=label
        if label is None:
            self.label = "unlabeled"
        for itraj,trajlam in enumerate(self.lams):
            if label is None:
                fname = os.path.join(directory,"dvdl_%.8f.dat"%(trajlam))
            else:
                fname = os.path.join(directory,"dvdl_%s_%.8f.dat"%(label,trajlam))
            if os.path.isfile(fname):
                fh=open(fname,"r")
                for line in fh:
                    cols=line.strip().split()
                    if len(cols) == 2:
                        time = float(cols[0])
                        dvdl = float(cols[1])
                        self.state[itraj].timesteps[time].dvdl = dvdl
            #print (fname,len(self.state[itraj].timesteps))

            for iparm,parmlam in enumerate(self.lams):
                if label is None:
                    fname = os.path.join(directory,"efep_%.8f_%.8f.dat"%(trajlam,parmlam))
                else:
                    fname = os.path.join(directory,"efep_%s_%.8f_%.8f.dat"%(label,trajlam,parmlam))
                if os.path.isfile(fname):
                    fh=open(fname,"r")
                    for line in fh:
                        cols=line.strip().split()
                        if len(cols) == 2:
                            time = float(cols[0])
                            ener = float(cols[1])
                            #print ("read %10.3f %20.10e"%(time,ener))
                            self.state[itraj].timesteps[time].ene[iparm] = ener
                    #print (fname,len(self.state[itraj].timesteps))

            print(("lam traj state %.8f has %i data pts"%(trajlam,len(self.state[itraj].timesteps))))

    def has_dvdl(self):
        has=True
        for s in self.state:
            if not s.has_dvdl():
                has=False
                break
        return has
    
        
    def get_beta(self,temperature=300.):
        kB_kcal_per_mol = 1.98720661358036e-03
        beta = 1. / ( kB_kcal_per_mol * temperature )

        #
        # this is what alchemical analysis uses, and the results seem correct
        # although, I can't explain why this isn't a bug
        # -- it looks like it is using kJ/mol instead of kcal/mol, but there must
        # be something funny internal to pymbar
        #
        #kB = 1.3806488*6.02214129/1000.0 # Boltzmann's constant (kJ/mol/K).
        #beta = 1./(kB*temperature)
        return beta


    def get_maxtime(self):
        maxtime = -1
        for s in self.state:
            #print (len( [ t for t in s.timesteps ] ))
            try:
                #print  (max( t for t in s.timesteps ))
                maxtime = max(maxtime, max( t for t in s.timesteps ))
            except:
                maxtime = max(maxtime,0)
        return maxtime


    def get_autocorrelation_times(self,tequil,tmax=1.e+10,dvdl_autocor=True):
            
        tau=[]
        for k in range(self.nlam):
            if dvdl_autocor:
                tau.append( self.state[k].get_correlation_time(tequil=tequil,tmax=tmax,state=None) )
            else:
                tau.append( self.state[k].get_correlation_time(tequil=tequil,tmax=tmax,state=k) )
        return tau
    
      
    def _get_mbar_arrays(self,time0,time1,temperature=300.,dvdl_autocor=True,uncorrelate=False):
        #print (time0,time1)
        if uncorrelate:
            states = []
            for k,s in enumerate(self.state):
                if dvdl_autocor:
                    states.append( s.uncorrelated_frames(tequil=time0,tmax=time1) )
                else:
                    states.append( s.uncorrelated_frames(tequil=time0,tmax=time1,state=k) )
        else:
            #print ("correlated")
            states = []
            for k,s in enumerate(self.state):
                states.append( s.extract_frames(tequil=time0,tmax=time1) )
        #print ([ len(s.timesteps) for s in states ])
        #
        # MAKE THE U_KLN ARRAY
        #
        beta = self.get_beta(temperature)
        N_k   = np.array( [ len(s.timesteps) for s in states ] )
        u_kln = np.zeros( [self.nlam,self.nlam,N_k.max()], np.float64 )
        dvdl  = np.zeros( [self.nlam,N_k.max()], np.float64 )
        nstates = len(states)

        mbar_ok=True
        for itraj,s in enumerate(states):
            times = sorted(s.timesteps)
            if len(times) == 0:
                mbar_ok=False
            else:
                t0 = times[0]
                ene = s.timesteps[t0].ene
                for e in ene:
                    if e is None:
                        mbar_ok=False
                        break
                if not mbar_ok:
                    break

        bar_ok=True
        iexp_ok=True
        dexp_ok=True
        if not mbar_ok:
            for itraj,s in enumerate(states):
                times = sorted(s.timesteps)
                if len(times) == 0:
                    bar_ok=False
                    if len(states) == 2:
                        if itraj == 0:
                            dexp_ok = False
                        else:
                            iexp_ok = False
                    else:
                        dexp_ok = False
                        iexp_ok = False
                            
                else:
                    t0 = times[0]
                    ene = s.timesteps[t0].ene
                    if itraj > 0:
                        if ene[itraj-1] is None:
                            bar_ok=False
                            iexp_ok=False
                            #break
                    if itraj + 1 < len(ene):
                        if ene[itraj+1] is None:
                            bar_ok=False
                            dexp_ok=False
                            #break
                #print ("GREP",itraj,len(times),bar_ok,iexp_ok,dexp_ok)


        if mbar_ok:
            for itraj,s in enumerate(states):
                tidx=0
                for time in sorted(s.timesteps):
                    step = s.timesteps[time]
                    for istate in range(len(step.ene)):
                        #print (time,istate,itraj,len(step.ene),step.ene[istate],step.ene[itraj])
                        u_kln[itraj,istate,tidx] = beta * (step.ene[istate]-step.ene[itraj])
                    dvdl[itraj,tidx] = step.dvdl
                    tidx += 1
        elif bar_ok:
            self.mbar=False
            for itraj,s in enumerate(states):
                tidx=0
                for time in sorted(s.timesteps):
                    step = s.timesteps[time]
                    for istate in range(max(itraj-1,0),min(len(step.ene),itraj+2)):
                        u_kln[itraj,istate,tidx] = beta * (step.ene[istate]-step.ene[itraj])
                    dvdl[itraj,tidx] = step.dvdl
                    tidx += 1
        else:
            self.mbar=False
            self.bar=False

            if self.iexp:
                self.iexp=iexp_ok
            if self.dexp:
                self.dexp=dexp_ok
                
            
            if self.iexp:
                for itraj,s in enumerate(states):
                    tidx=0
                    for time in sorted(s.timesteps):
                        step = s.timesteps[time]
                        for istate in range(max(itraj-1,0),min(len(step.ene),itraj+1)):
                            u_kln[itraj,istate,tidx] = beta * (step.ene[istate]-step.ene[itraj])
                        dvdl[itraj,tidx] = step.dvdl
                        tidx += 1
            if self.dexp:
                for itraj,s in enumerate(states):
                    tidx=0
                    #print (itraj,len(s.timesteps))
                    for time in sorted(s.timesteps):
                        step = s.timesteps[time]
                        for istate in range(max(itraj,0),min(len(step.ene),itraj+2)):
                            u_kln[itraj,istate,tidx] = beta * (step.ene[istate]-step.ene[itraj])
                        dvdl[itraj,tidx] = step.dvdl
                        tidx += 1
            else:
                for itraj,s in enumerate(states):
                    tidx=0
                    for time in sorted(s.timesteps):
                        dvdl[itraj,tidx] = s.timesteps[time].dvdl
                        tidx += 1
                    
        #
        # MAKE THE U_KN ARRAY
        #
        # N_k   = np.array( [ len(s.timesteps) for s in states ] )
        # N_max = N_k.sum()
        # u_kn  = np.zeros( [self.nlam,N_max], np.float64 )
        # idx=0
        # for itraj,s in enumerate(states):
        #     for time in sorted(s.timesteps):
        #         step = s.timesteps[time]
        #         for istate in range(len(step.ene)):
        #             u_kn[istate,idx] = beta * (step.ene[istate]-step.ene[itraj])
        #         idx += 1
        
        return dvdl,u_kln,N_k
    
            
    def get_energies(self,tequil=0,tmax=1.e+10,temperature=300.,dvdl_autocor=True,overlap=""):
        return self._get_energies_in_range(tequil,tmax,temperature,dvdl_autocor,overlap=overlap)
        
    
    def get_forward_timeseries(self,tequil=0,dtime=10000,temperature=300.,dvdl_autocor=True,verbose=False):
        tseries = ddict( float )
        maxtime = self.get_maxtime()
        if tequil >= maxtime:
            return tseries
        
        if True:
            ntimes = 0  
            t0=tequil
            while t0 < maxtime:
                t1    = t0 + dtime
                time  = min( t1, maxtime )
                t0    = t1
                ntimes += 1

            itime = 0
            t0=tequil
            while t0 < maxtime:
                t1    = t0 + dtime
                time  = min( t1, maxtime )-tequil
                itime += 1
                if verbose:
                    print(("%16s time series %8.3f ns -> %8.3f ns (%3.0f%% complete)"%\
                        ( self.label,
                          tequil/1000.,
                          t1/1000.,
                          100.*float(itime)/float(ntimes) )))
                tseries[time],dvdl,taus = self._get_energies_in_range(tequil,t1,temperature,dvdl_autocor)
                t0 = t1

        return tseries



            
    def get_reverse_timeseries(self,tequil=0,dtime=10000,temperature=300.,dvdl_autocor=True,verbose=False):
        tseries = ddict( float )
        maxtime = self.get_maxtime()
        if tequil >= maxtime:
            return tseries
        if True:
            ntimes = 0  
            t0=maxtime
            while t0 > tequil:
                t1    = t0 - dtime
                time  = max( t1, tequil )
                t0    = t1
                ntimes += 1

            itime = 0
            t0=maxtime
            while t0 > tequil:
                t1    = t0 - dtime
                time  = maxtime - max( t1, tequil )
                itime += 1
                if verbose:
                    print(("%16s time series %8.3f ns -> %8.3f ns (%3.0f%% complete)"%\
                        ( self.label,
                          t1/1000.,
                          maxtime/1000.,
                          100.*float(itime)/float(ntimes) )))
                tseries[time],dvdl,taus = self._get_energies_in_range(t1,maxtime,temperature,dvdl_autocor)
                t0 = t1
        return tseries

    
    def get_segmented_timeseries(self,tequil=0,dtime=10000,temperature=300.,dvdl_autocor=True,verbose=False):
        tseries = ddict( float )
        maxtime = self.get_maxtime()
        if tequil >= maxtime:
            return tseries
        
        if True:
            ntimes = 0  
            t0=tequil
            while t0 < maxtime:
                t1    = t0 + dtime
                time  = min( t1, maxtime )
                t0    = t1
                ntimes += 1

            itime = 0
            t0=tequil
            while t0 < maxtime:
                t1    = t0 + dtime
                time  = min( t1, maxtime )-tequil
                itime += 1
                if verbose:
                    print(("%16s time series %8.3f ns -> %8.3f ns (%3.0f%% complete)"%\
                        ( self.label,
                          t0/1000.,
                          t1/1000.,
                          100.*float(itime)/float(ntimes) )))
                tseries[time],dvdl,taus = self._get_energies_in_range(t0,t1,temperature,dvdl_autocor)
                t0 = t1

        return tseries



    def _get_energies_in_range( self, tequil, tmax, temperature, dvdl_autocor, overlap="" ):
        beta = self.get_beta(temperature)
        energies = ddict(tuple)
        
        cdvdl,u,cn = self._get_mbar_arrays( tequil, tmax, \
                                            temperature=temperature,\
                                            dvdl_autocor=dvdl_autocor,\
                                            uncorrelate=False )
        cor = self._get_mbar_energies( u, cn, beta )
        
        udvdl,u,un = self._get_mbar_arrays( tequil, tmax, \
                                            temperature=temperature,\
                                            dvdl_autocor=dvdl_autocor,\
                                            uncorrelate=True )
        uncor = self._get_mbar_energies( u, un, beta, overlap=overlap )

        for method in cor:
            energies[method] = ( cor[method][0], uncor[method][1], cor[method][2], uncor[method][3] )

        dvdl=None
        if self.ti or self.ti3:
            sigs = np.zeros( self.nlam )
            avgs = np.zeros( self.nlam )
            for k in range(self.nlam):
                if un[k] > 1:
                    # shouldn't sigs be udvdl ... for the uncorrelated data?
                    sigs[k] = np.sqrt( np.var( udvdl[k,0:un[k]], ddof=1 ) / (un[k]) )
                elif cn[k] > 1:
                    sigs[k] = np.sqrt( np.var( cdvdl[k,0:cn[k]], ddof=1 ) / (cn[k]) )
                if cn[k] > 0:
                    avgs[k] = np.mean( cdvdl[k,0:cn[k]] )
            #sigs = np.array( [ np.sqrt( np.var( cdvdl[k,0:cn[k]], ddof=1 ) / (cn[k]) ) for k in range(self.nlam) ] )
            #avgs = np.array( [ np.mean( cdvdl[k,0:cn[k]] ) for k in range(self.nlam) ] )

            dvdl = [ (avgs[k],sigs[k]) for k in range(self.nlam) ]
            if len(dvdl) == 2:
                # if 2 states and 1 has no data, then integrate a rectangle as a 1-sided approx
                if cn[0] == 0:
                    dvdl[0]=dvdl[1]
                    avgs[0]=avgs[1]
                    sigs[0]=sigs[1]
                if cn[1] == 0:
                    dvdl[1]=dvdl[0]
                    avgs[1]=avgs[0]
                    sigs[1]=sigs[0]

        if self.ti:
            energies["TI"] = CalcTI_trapez( self.lams, avgs, sigs )

        if self.ti3:
            energies["TI3"] = CalcTI_spline( self.lams, avgs, sigs )

        taus = self.get_autocorrelation_times(tequil,tmax=tmax,dvdl_autocor=dvdl_autocor)


        if self.global_pka_units:
            c = self.global_pka_per_kcal
            for method in energies:
                x = energies[method]
                f = c*x[0]
                df = c*x[1]
                fs = []
                for y in x[2]:
                    if y is not None:
                        fs.append( y*c )
                    else:
                        fs.append( None )
                dfs = []
                for y in x[3]:
                    if y is not None:
                        dfs.append( y*c )
                    else:
                        dfs.append( None )
                energies[method] = (f,df,fs,dfs)
#            if dvdl is not None:
#                for i in range(len(dvdl)):
#                    dvdl[i] = ( dvdl[i][0]*c, dvdl[i][1]*c )
                        
        return energies,dvdl,taus
        

        
    def _get_mbar_energies( self, u_kln, N_k, beta, overlap="" ):
        import math
        energies = ddict( tuple )

        if self.iexp:
            # IEXP
            bar_fs  = [0.]*self.nlam
            bar_dfs = [0.]*self.nlam
            for k in range(self.nlam-1):
                #w_F = u_kln[k,k+1,0:N_k[k]]   - u_kln[k,k,0:N_k[k]]
                w_R = u_kln[k+1,k,0:N_k[k+1]] - u_kln[k+1,k+1,0:N_k[k+1]]
                f,df = pymbar.EXP(w_R)
                f /= beta
                df /= beta
                bar_fs[k+1] = bar_fs[k] - f
                bar_dfs[k+1] = np.sqrt( bar_dfs[k]**2 + df**2 )
            bar_f=bar_fs[-1]
            bar_df=bar_dfs[-1]
            if math.isnan(bar_df):
                bar_df=0.
            energies["IEXP"]  = ( bar_f, bar_df, bar_fs, bar_dfs )
            
        if self.dexp:
            # DEXP
            bar_fs  = [0.]*self.nlam
            bar_dfs = [0.]*self.nlam
            for k in range(self.nlam-1):
                w_F = u_kln[k,k+1,0:N_k[k]]   - u_kln[k,k,0:N_k[k]]
                #w_R = u_kln[k+1,k,0:N_k[k+1]] - u_kln[k+1,k+1,0:N_k[k+1]]
                f,df = pymbar.EXP(w_F)
                f /= beta
                df /= beta
                bar_fs[k+1] = bar_fs[k] + f
                bar_dfs[k+1] = np.sqrt( bar_dfs[k]**2 + df**2 )
            bar_f=bar_fs[-1]
            bar_df=bar_dfs[-1]
            if math.isnan(bar_df):
                bar_df=0.
            energies["DEXP"]  = ( bar_f, bar_df, bar_fs, bar_dfs )

            
            
        if self.mbar:

            ok = True
            if u_kln.shape[2] < 2:
                ok = False
            else:
                if ( (u_kln[:,:,0]-u_kln[:,:,-1])**2 ).sum() < 1.e-15:
                    ok = False
            if ok:
                try:
                    mbar = pymbar.MBAR( u_kln, N_k, initialize='BAR', relative_tolerance=1.e-10 )
                    (Deltaf_ij, dDeltaf_ij) = mbar.getFreeEnergyDifferences( uncertainty_method='svd-ew', return_theta = False )
                except Exception as e:
                    mbar = pymbar.MBAR( u_kln, N_k, initialize='zeros', relative_tolerance=1.e-10 )
                    (Deltaf_ij, dDeltaf_ij) = mbar.getFreeEnergyDifferences( uncertainty_method='svd-ew', return_theta = False )

                #print ("df",Deltaf_ij[0,20])

                Deltaf_ij  /= beta
                dDeltaf_ij /= beta
                mbar_fs     = Deltaf_ij[0,:]
                mbar_dfs    = [0]*self.nlam
                for k in range(self.nlam-1):
                    mbar_dfs[k+1] = np.sqrt( dDeltaf_ij[k,k+1]**2 + mbar_dfs[k]**2 )
                mbar_f = mbar_fs[-1]
                mbar_df= mbar_dfs[-1]
                energies["MBAR"] = ( mbar_f, mbar_df, mbar_fs, mbar_dfs )
                if len(overlap) > 0:
                    try:
                        results = mbar.computeOverlap()
                        O = results['matrix']
                        plotOverlapMatrix(O,self.nlam,overlap)
                    except:
                        print ("Could not calculate overlap matrix")

            else:
                energies["MBAR"] = ( 0, 0, [0]*self.nlam, [0]*self.nlam )

                    
        
        if self.bar:
            bar_fs  = [0.]*self.nlam
            bar_dfs = [0.]*self.nlam
            for k in range(self.nlam-1):
                f=0
                df=0
                if N_k[k] > 1:
                    w_F = u_kln[k,k+1,0:N_k[k]]   - u_kln[k,k,0:N_k[k]]
                    w_R = u_kln[k+1,k,0:N_k[k+1]] - u_kln[k+1,k+1,0:N_k[k+1]]
                    ok=True
                    if abs(w_F[0]-w_F[-1]) < 1.e-15 and abs(w_R[0]-w_R[-1]) < 1.e-15:
                        ok=False
                    if ok:
                        f,df = pymbar.BAR(w_F,w_R)
                f /= beta
                df /= beta
                bar_fs[k+1] = bar_fs[k] + f
                bar_dfs[k+1] = np.sqrt( bar_dfs[k]**2 + df**2 )
            bar_f=bar_fs[-1]
            bar_df=bar_dfs[-1]
            if math.isnan(bar_df):
                bar_df=0.
            energies["BAR"]  = ( bar_f, bar_df, bar_fs, bar_dfs )
        
        return energies


    def write_header(self,fh,methods=None):
        if methods is None:
            methods = self.get_methods()
        fh.write("%16s"%("label"))
        for method in methods:
            fh.write("%12s%12s"%(method,""))
        fh.write("%12s"%("max(tau)"))
        fh.write("\n")
        
    def write_row(self,fh,energies,taus,methods=None):
        if methods is None:
            methods = self.get_methods()
            
        fh.write("%16s"%(self.label))
        for method in methods:
            if method in energies:
                f=energies[method][0]
                df=energies[method][1]
                fh.write("%12.5f +- %8.5f"%(f,df))
            else:
                fh.write("%12s +- %8s"%("na","na"))
        fh.write("%12.0f"%(max(taus)))
        fh.write("\n")

    def generic_run(self,tequil,dtime=None,outputdir="./"):
        import os
        import sys
        
        methods = self.get_methods()
        energies,dvdl,taus = self.get_energies(tequil=tequil)
        #taus = self.get_autocorrelation_times(tequil=tequil)
                 
        fh = sys.stdout
        self.write_header(fh,methods)
        self.write_row(fh,energies,taus)

        if self.ti or self.ti3:
            fname = os.path.join( outputdir, "dvdl_%s.dat"%(self.label) )
            fh = open(fname,"w")
            for k in range(self.nlam):
                lam = self.lams[k]
                val,err = dvdl[k]
                fh.write("%.8f %19.10e %19.10e\n"%(lam,val,1.96*err))
            fh.close()

        if dtime is not None:
            if dtime <= 0:
                pass
            else:
                ts = self.get_forward_timeseries(tequil=tequil,dtime=dtime,verbose=True)
                times = [t for t in sorted(ts)]
                for method in methods:
                    fname = os.path.join( outputdir, "fwdGvsT_%s_%s.dat"%(self.label,method) )
                    fh = open(fname,"w")
                    for time in times:
                        val,err = ts[time][method]
                        fh.write("%10.3f %15.8f %15.8f\n"%(time/1000.,val,1.96*err))
                    fh.write("\n")
                fh.close()
                    



class DataLoc(object):
    def __init__(self,name,inpdir,filelabel="",tequil=0,tmax=1.e+10):
        self.name = name
        self.inpdir = inpdir
        self.filelabel = filelabel
        self.tequil = tequil
        self.tmax = tmax

class DataAnalysis(object):
    def __init__(self,ene_dvdl_taus=None):
        if ene_dvdl_taus is None:
            self.ene = ddict( tuple )
            self.dvdl = []
            self.taus = []
            self.maxtau = None
        else:
            self.ene=ene_dvdl_taus[0]
            self.dvdl=ene_dvdl_taus[1]
            self.taus=ene_dvdl_taus[2]
            self.maxtau = max(self.taus)


class WuKofkeAnalysis(object):
    def __init__(self):
        self.piab = 0
        self.piba = 0
        self.Kab = 0
        self.Kba = 0
        self.PercentOverlap = 0
        self.HistAB = []
        self.HistBA = []
        self.GauAB = []
        self.GauBA = []
        self.GauXs = []
            

class TransformTrial(object):
    def __init__(self,transform_name,dataloc,methods=["TI","TI3","BAR","MBAR","IEXP","DEXP"]):
        self.name = transform_name + "." + dataloc.name
        self.dataloc = dataloc
        self.methods = methods
        self.read()
        
    def read(self):
        import re
        import os.path
        import glob
        from collections import defaultdict as ddict
        
        dvdls=[]
        efeps=[]
        for f in glob.glob( os.path.join(self.dataloc.inpdir,"dvdl_%s*.dat"%(self.dataloc.filelabel)) ):
            dvdls.append(f)
        for f in glob.glob( os.path.join(self.dataloc.inpdir,"efep_%s*.dat"%(self.dataloc.filelabel)) ):
            efeps.append(f)

        self.lams = []
        for f in dvdls:
            dname,fname = os.path.split(f)
            m = re.match( r"dvdl(.*)_([0-9\.]+?)\.dat", fname )
            if m:
                label = m.group(1)
                if len(label) > 0:
                    if label[0] == "_":
                        label = label[1:]
                self.lams.append( m.group(2) )
        for f in efeps:
            dname,fname = os.path.split(f)
            m = re.match( r"efep(.*)_([0-9\.]+?)_([0-9\.]+?)\.dat", fname )
            if m:
                label = m.group(1)
                if len(label) > 0:
                    if label[0] == "_":
                        label = label[1:]
                self.lams.append( m.group(2) )
                self.lams.append( m.group(3) )

        self.lams = [ float(lam) for lam in set(self.lams) ]
        self.lams = list(set(self.lams))
        self.lams.sort()
                
        print ("DATA FILES LOCATED:\n")
        print(("name:     ",self.dataloc.name))
        print(("directory:",self.dataloc.inpdir))
        print(("label:    ",self.dataloc.filelabel))
        print(("lambdas:  ",self.lams))
        print ("")

        self.job = AlchemicalTransform(self.lams,
                                       directory=self.dataloc.inpdir,
                                       label=self.dataloc.filelabel)

        self.dvdl_autocor=True
        if "TI" not in self.methods and "TI3" not in self.methods:
            self.job.ti  = False
            self.job.ti3 = False
            self.dvdl_autocor=False
        if "TI"   not in self.methods:
            self.job.ti   = False
        if "TI3"  not in self.methods:
            self.job.ti3  = False
        if "BAR"  not in self.methods:
            self.job.bar  = False
        if "MBAR" not in self.methods:
            self.job.mbar = False
        if "IEXP" not in self.methods:
            self.job.iexp = False
        if "DEXP" not in self.methods:
            self.job.dexp = False

        if not self.job.has_dvdl():
            print(("WARNING: directory='%s' label='%s'"%(self.dataloc.inpdir,self.dataloc.filelabel)))
            print ("Some dvdl data is missing; skipping TI and TI3 analysis and dvdl autocorrelations")
            self.job.ti  = False
            self.job.ti3 = False
            self.dvdl_autocor=False

    def get_overlap_plot_name(self,plotdir=""):
        name=""
        if len(plotdir) > 0 and "MBAR" in self.methods:
            if len(self.dataloc.filelabel) > 0:
                name = ".".join( splitall( self.dataloc.inpdir ) ) + ".%s.S.eps"%(self.dataloc.filelabel)
            else:
                name = ".".join( splitall( self.dataloc.inpdir ) ) + ".S.eps"
            #name = os.path.join(plotdir,name)
        return name
            
    def analyze(self,tequil=0,tmax=1.e+10,plotdir=""):
        import os
        tmax  = min( self.job.get_maxtime(), tmax )
        print(("tequil,tmax=",tequil,tmax,self.job.get_maxtime()))
        name = self.get_overlap_plot_name(plotdir)
        if len(name) > 0:
            name = os.path.join(plotdir,name)
        self.full = self.get_energies( tequil=tequil, tmax=tmax, plotdir=name )
        #print ("dvdl=",len(self.full.dvdl))
        
        #self.taus = self.get_taus(tequil=tequil)
        #self.maxtau = max(self.taus)
        delta = (tmax-tequil)/12.
        #print (tmax,tequil)
        self.times   = []
        self.forward = []
        self.reverse = []
        self.segment = []

        dvdlname=""
        vals = [ len( s.timesteps ) for s in self.job.state ]
        vals = [ v for v in vals if v > 0 ]

        #print ("delta=",delta)
        
        if min( vals ) < 12:
            for i in range(12):
                self.times.append( (i+1)*delta )
                self.forward.append( self.full )
                self.reverse.append( self.full )
                self.segment.append( self.full )
        else:
            for i in range(12):
                self.times.append( (i+1)*delta )
                self.forward.append( self.get_energies(tequil=tequil+0*delta,tmax=tequil+(i+1)*delta) )
                #print ("forward",i,len(self.forward[-1].dvdl))
                self.reverse.append( self.get_energies(tequil=tequil+(11-i)*delta,tmax=tequil+12*delta) )
                #print ("reverse",i,len(self.reverse[-1].dvdl))
                self.segment.append( self.get_energies(tequil=tequil+i*delta,tmax=tequil+(i+1)*delta) )
                #print ("segment",i,len(self.segment[-1].dvdl))

            
    def get_methods(self):
        return self.job.get_methods()
            
    def get_energies(self,tequil=0,tmax=1.e+10,plotdir=""):
        return DataAnalysis(self.job.get_energies(tequil=tequil,tmax=tmax,dvdl_autocor=self.dvdl_autocor,overlap=plotdir))
    
    def get_taus(self,tequil=0,tmax=1.e+10):
        return self.job.get_autocorrelation_times(tequil=tequil,tmax=tmax,dvdl_autocor=self.dvdl_autocor)

    def get_wukofke_metrics(self,tequil=0,tmax=1.e+10):
        from scipy import special
        #import scipy.special.erf
        import numpy as np
        from scipy import integrate
        
        if not self.job.mbar and not self.job.bar:
            return None
        
        nlam = self.job.nlam
        beta = self.job.get_beta(300.)
        states=[]
        for k,s in enumerate(self.job.state):
            states.append( s.uncorrelated_frames(tequil=tequil,tmax=tmax,state=k) )
            #states.append( s.extract_frames(tequil=tequil,tmax=tmax) )
 
        if "BAR" in self.methods:
            m="BAR"
        elif "MBAR" in self.methods:
            m="MBAR"

        results = []
        for ilam in range(nlam-1):
            results.append( WuKofkeAnalysis() )
            
        for ilam in range(nlam-1):
            #dAab = self.full.ene[m][2][ilam+1] - self.full.ene[m][2][ilam]
            #dAba = self.full.ene[m][2][ilam]   - self.full.ene[m][2][ilam+1]
            dUabs = []
            dUbas = []
            paa=[]
            pab=[]
            pbb=[]
            pba=[]
            na = len( states[ilam].timesteps )
            nb = len( states[ilam+1].timesteps )
            for time in sorted(states[ilam].timesteps):
                step = states[ilam].timesteps[time]
                dUabs.append(  (step.ene[ilam+1]-step.ene[ilam]) )
                paa.append(  step.ene[ilam] )
                pab.append( step.ene[ilam+1] )
            for time in sorted(states[ilam+1].timesteps):
                step = states[ilam+1].timesteps[time]
                dUbas.append(  (step.ene[ilam]-step.ene[ilam+1]) )
                pbb.append(  step.ene[ilam+1] )
                pba.append(  step.ene[ilam] )


            dUab = np.mean( dUabs )
            dUba = np.mean( dUbas )
            Fab = -( np.log( np.average([np.exp(-beta*(U-dUab)) for U in dUabs]) ) -beta*dUab ) / beta
            Fba = -( np.log( np.average([np.exp(-beta*(U-dUba)) for U in dUbas]) ) -beta*dUba ) / beta

            #print ("%12.3f %12.3f %12.3f %12.3f"%(dUab,dUba,Fab,Fba))

            Fa2b = 0.5 * (Fab-Fba)
            Fb2a = -Fa2b

            sa = beta*(dUab - Fa2b)
            sb = beta*(dUba - Fb2a)
            Wab = float(special.lambertw( (na-1)**2 / (2.*np.pi) ))
            Wba = float(special.lambertw( (na-1)**2 / (2.*np.pi) ))
            #print ("sa,sb,Wab",sa,sb,Wab)
            piab = np.sqrt( (sa/sb)*Wab ) - np.sqrt(2.*sa)
            piba = np.sqrt( (sb/sa)*Wba ) - np.sqrt(2.*sb)
            results[ilam].piab = piab
            results[ilam].piba = piba

            za = 0.5 / np.var( paa, ddof=1 )
            ca = np.mean( paa )
            zb = 0.5 / np.var( pba, ddof=1 )
            cb = np.mean( pba )
            nb = np.sqrt( zb/np.pi )
            na = np.sqrt( za/np.pi )
            f = lambda x: (1.+special.erf(x*np.sqrt(za))) * nb * np.exp(-zb*(x-cb+ca)**2)
            k1 = integrate.quad( f, -np.inf, np.inf )[0]
            f = lambda x: (1.+special.erf(x*np.sqrt(zb))) * na * np.exp(-za*(x-ca+cb)**2)
            k2 = integrate.quad( f, -np.inf, np.inf )[0]

            results[ilam].Kba = 0.5*(k1+k2)
            
            za = 0.5 / np.var( pbb, ddof=1 )
            ca = np.mean( pbb )
            zb = 0.5 / np.var( pab, ddof=1 )
            cb = np.mean( pab )
            nb = np.sqrt( zb/np.pi )
            na = np.sqrt( za/np.pi )
            f = lambda x: (1.+special.erf(x*np.sqrt(za))) * nb * np.exp(-zb*(x-cb+ca)**2)
            k1 = integrate.quad( f, -np.inf, np.inf )[0]
            f = lambda x: (1.+special.erf(x*np.sqrt(zb))) * na * np.exp(-za*(x-ca+cb)**2)
            k2 = integrate.quad( f, -np.inf, np.inf )[0]
            results[ilam].Kab = 0.5*(k1+k2)
           
            dUbas = [ -x for x in dUbas ]
            cb = np.mean( dUbas )
            ca = np.mean( dUabs )
            cc = (ca+cb)/2.
            dUabs = [ x-cc for x in dUabs ]
            dUbas = [ x-cc for x in dUbas ]
            cb = np.mean( dUbas )
            ca = np.mean( dUabs )
            za = 0.5 / np.var( dUabs, ddof=1 )
            zb = 0.5 / np.var( dUbas, ddof=1 )
            na = np.sqrt(za/np.pi)
            nb = np.sqrt(zb/np.pi)
            zz = za*zb/(za+zb)
            den = max( np.sqrt( 0.5*za/np.pi ), np.sqrt( 0.5*zb/np.pi ) )
            results[ilam].PercentOverlap = 100 * np.sqrt(zz/np.pi) * np.exp(-zz*(ca-cb)**2) / den
            

            results[ilam].HistAB = np.histogram( dUabs, bins=max(1,min(15,len(dUabs)//10)), density=True )
            results[ilam].HistBA = np.histogram( dUbas, bins=max(1,min(15,len(dUbas)//10)), density=True )

            xlo = min( results[ilam].HistAB[1][0], results[ilam].HistBA[1][0] )
            xhi = max( results[ilam].HistAB[1][-1], results[ilam].HistBA[1][-1] )
            results[ilam].GauXs = np.linspace( xlo,xhi,100 )
            results[ilam].GauAB = np.array([ na*np.exp(-za*(x-ca)**2) for x in results[ilam].GauXs ])
            results[ilam].GauBA = np.array([ nb*np.exp(-zb*(x-cb)**2) for x in results[ilam].GauXs ])
            klo=0
            khi=0
            leftside=True
            for k in range(len(results[ilam].GauXs)):
                g = max( results[ilam].GauAB[k],results[ilam].GauBA[k] )
                if leftside:
                    if g < 0.1:
                        klo=k
                    else:
                        klo=k
                        leftside=False
                else:
                    if g < 0.1:
                        khi=k
                        break
            results[ilam].xlo = results[ilam].GauXs[klo]
            results[ilam].xhi = results[ilam].GauXs[khi]


            
        return results
    
            
    
    def write_plot(self,directory,latex=None,ncol=5):
        import os
        from matplotlib import pyplot as plt
        trial = self
        nlam = trial.job.nlam
        nrow = nlam//ncol
        if nlam % ncol != 0:
            nrow += 1
        fig = plt.figure()
        fig.set_size_inches( 7, (1) * (7./3), forward=True )
        method = "TI3"
        if method not in trial.forward[0].ene:
            for alts in [ "BAR","IEXP","DEXP" ]:
                if alts in trial.forward[0].ene:
                    method = alts
                    break
        if self.full.dvdl is None and method == "TI3":
            method = "BAR"

        dvdlname=""
            
        if True:
            ax  = plt.subplot(1,3,1)
            ymin=1.e+10
            ymax=-1.e+10
    
            n   = len(trial.times)
            x   = [ time/1000. for time in trial.times ]
            
            ymin,ymax = plot_energy_timeseries(ax,trial,method)
            props = dict(boxstyle='round', facecolor='white', linewidth=0, alpha=0.5)
            ax.text(0.5,0.975,method,transform=ax.transAxes,
                    verticalalignment='top',horizontalalignment='center')
            ax.text(0.5,0.01,r"max $\tau$=%.0f ps"%(max(trial.full.taus)),transform=ax.transAxes,
                    verticalalignment='bottom',horizontalalignment='center')

            ax = plt.subplot(1,3,2)
            x = []
            y = []

            if trial.full.dvdl is not None:
                x  = trial.lams
                #if trial.forward is not None:
                #    y  = [ trial.forward[5].dvdl[i][0] for i in range(len(x)) ]
                #    dy = [ trial.forward[5].dvdl[i][1] for i in range(len(x)) ]
                #else:
                y  = [ trial.full.dvdl[i][0] for i in range(len(x)) ]
                dy = [ trial.full.dvdl[i][1] for i in range(len(x)) ]
                ax.errorbar(x,y,yerr=dy,capsize=1,c='k',linestyle='solid',linewidth=1,marker='o',ms=2)

                # ----------------
                # Write dat file
                if len(trial.dataloc.filelabel) > 0:
                    dvdlname = ".".join( splitall( trial.dataloc.inpdir ) ) + ".%s.DVDL.dat"%(trial.dataloc.filelabel)
                else:
                    dvdlname = ".".join( splitall( trial.dataloc.inpdir ) ) + ".DVDL.dat"
                fname = os.path.join( directory, dvdlname )
                datafh = open(fname,"w")
                for ii in range(len(x)):
                    datafh.write("%18.8e %18.8e %18.8e\n"%(x[ii],y[ii],dy[ii]))
                datafh.close()
                # ----------------
                
            ax.set_xlim(-0.1,1.1)
            #if UsesPkaUnits():
            #    ax.set_ylabel(r"$\langle\partial U/\partial\lambda\rangle_{\lambda}$ ($\Delta$pK$_{\text{a}}$ units)")
            #else:
            ax.set_ylabel(r"$\langle\partial U/\partial\lambda\rangle_{\lambda}$ (kcal/mol)")
            ax.set_xlabel(r"$\lambda$")
            ax.grid(linestyle='dotted',linewidth=0.5)



            
            ax  = plt.subplot(1,3,3)
            x  = trial.lams
            y  = trial.full.taus
            ax.set_xlim(-0.1,1.1)
            ax.set_ylabel(r"$\tau$ (ps)")
            ax.set_xlabel(r"$\lambda$")
            ax.grid(linestyle='dotted',linewidth=0.5)
            ax.plot(x,y,c='k',linestyle='solid',linewidth=1,marker='o',ms=2)

            if len(trial.dataloc.filelabel) > 0:
                gname = ".".join( splitall( trial.dataloc.inpdir ) ) + ".%s.GvsT.eps"%(trial.dataloc.filelabel)
            else:
                gname = ".".join( splitall( trial.dataloc.inpdir ) ) + ".GvsT.eps"
            fname = os.path.join( directory, gname )
            plt.tight_layout()
            plt.savefig(fname,dpi=300)
            fig.clf()


            
            ## GvsL

            fig = plt.figure()
            fig.set_size_inches( 7, (1) * (7./1), forward=True )
            ax  = plt.subplot(1,1,1)
            ax.set_xlim(-0.1,1.1)
            if UsesPkaUnits():
                ax.set_ylabel(r"$\Delta$pK${}_{\mathregular{a}}$")
            else:
                ax.set_ylabel(r"Free energy (kcal/mol)")
            ax.set_xlabel(r"$\lambda$")
            methods = trial.get_methods()
            nm = len(methods)
            dl = (trial.lams[1]-trial.lams[0])
            dx = dl/(nm+2)
            if nm % 2 == 1:
                xoff = [ dx*(i-(nm+1)/2) for i in range(nm) ]
            else:
                xoff = [ dx*(i-(nm+1)/2+0.5) for i in range(nm) ]
            for im,method in enumerate(methods):
                xs = np.array(trial.lams) + xoff[im]
                ys = np.array(trial.full.ene[method][2])

                for i in range(nlam-2):
                    ys[ nlam-1-i ] = ys[ nlam-1-i ] - ys[ nlam-2-i ]
                
                ds = np.array(trial.full.ene[method][3]) * 1.96
                #print ("")
                #print (method,xs)
                #print (method,ys)
                #print (method,ds)
                #print ("")
                ax.bar(xs,ys,width=dx,yerr=ds,label=method)
            plt.legend(loc='best')

            if len(trial.dataloc.filelabel) > 0:
                glname = ".".join( splitall( trial.dataloc.inpdir ) ) + ".%s.GvsL.eps"%(trial.dataloc.filelabel)
            else:
                glname = ".".join( splitall( trial.dataloc.inpdir ) ) + ".GvsL.eps"
            fname = os.path.join( directory, glname )
            plt.tight_layout()
            plt.savefig(fname,dpi=300)
            fig.clf()

            ## overlap
            sname = trial.get_overlap_plot_name(directory)

            
            
            ## LvsT
            
            fig = plt.figure()
            fig.set_size_inches( 7, (nrow) * (7./ncol), forward=True )
            has_dvdl = trial.forward[0].dvdl is not None
            if has_dvdl:
                rowlams = []
                for row in range(nrow):
                    rowlams.append( [] )
                for ilam in range(nlam):
                    row = ilam // ncol
                    rowlams[row].append( ilam )
                for row in range(nrow):
                    ymin=1.e+10
                    ymax=-1.e+10
                    if trial.forward[0].dvdl is not None:
                        ymin,ymax = get_dvdl_minmax_from_trial(ymin,ymax,trial,rowlams[row])
                        ymax += 0.05 * (ymax-ymin)
                    else:
                        ymin=0
                        ymax=1
                        
                    oax = None
                    for ilam in rowlams[row]:
                        if oax is None:
                            ax  = plt.subplot(nrow,ncol,ilam+1)
                            oax = ax
                        else:
                            ax  = plt.subplot(nrow,ncol,ilam+1,sharex=oax,sharey=oax)
                        ax.tick_params(axis='both',which='major',labelsize=9)
                        n   = len(trial.times)
                        x   = [ time/1000. for time in trial.times ]
                        y   = [ trial.forward[i].dvdl[ilam][0] for i in range(n) ]
                        dy  = [ 1.96 * trial.forward[i].dvdl[ilam][1] for i in range(n) ]
                        ax.errorbar(x,y,yerr=dy,capsize=1,c='k',linestyle='solid',linewidth=1,marker='o',ms=2)
                        y   = [ trial.reverse[i].dvdl[ilam][0] for i in range(n) ]
                        dy  = [ 1.96 * trial.reverse[i].dvdl[ilam][1] for i in range(n) ]      
                        ax.errorbar(x,y,yerr=dy,capsize=1,c='r',linestyle='solid',linewidth=1,marker='o',ms=2)
                        y   = [ trial.segment[i].dvdl[ilam][0] for i in range(n) ]
                        dy  = [ 1.96 * trial.segment[i].dvdl[ilam][1] for i in range(n) ]      
                        ax.errorbar(x,y,capsize=0,c='g',linestyle='dotted',linewidth=1,marker='x',ms=2)

                        props = dict(boxstyle='round', facecolor='white', linewidth=0, alpha=0.5)
                        ax.text(0.5,0.975,"$\lambda$=%.8f"%(trial.lams[ilam]),transform=ax.transAxes,verticalalignment='top',horizontalalignment='center',fontsize=9)
                        ax.text(0.5,0.01,r"$\tau$=%.0f ps"%(trial.full.taus[ilam]),transform=ax.transAxes,verticalalignment='bottom',horizontalalignment='center',fontsize=9)
                    
                        ax.set_xlim(0,x[-1]+x[0])
                        ax.set_ylim(ymin,ymax)
                        if ilam%ncol != 0:
                            plt.setp(ax.get_yticklabels(),visible=False)
                        if ilam+ncol < nlam:
                            plt.setp(ax.get_xticklabels(),visible=False)
                        if ilam%ncol == 0:
                            #if UsesPkaUnits():
                            #    ax.set_ylabel("$\partial U/\partial\lambda$ ($\Delta$pK$_{\text{a}}$ units)")
                            #else:
                            ax.set_ylabel("$\partial U/\partial\lambda$ (kcal/mol)",fontsize=9)
                        if ilam/ncol+1 == nrow:
                            ax.set_xlabel("Time (ns)",fontsize=9)
                        ax.grid(linestyle='dotted',linewidth=0.5)
                
                plt.tight_layout()
                fig.subplots_adjust(wspace=0, hspace=0)
                if len(trial.dataloc.filelabel) > 0:
                    dvdlname = ".".join( splitall( trial.dataloc.inpdir ) ) + ".%s.DVDLvsT.eps"%(trial.dataloc.filelabel)
                else:
                    dvdlname = ".".join( splitall( trial.dataloc.inpdir ) ) + ".DVDLvsT.eps"
                fname = os.path.join( directory, dvdlname )
                plt.savefig(fname,dpi=300)
                fig.clf()

                ##################################


            if "BAR" in self.get_methods() or "MBAR" in self.get_methods():

                wk = self.get_wukofke_metrics()
                
                nrow = (nlam-1) // ncol
                if nlam % ncol != 0:
                    nrow += 1
                if nrow < 1:
                    nrow = 1
                fig = plt.figure()
                fig.set_size_inches( 7, (nrow) * (7./ncol), forward=True )

                rowlams = []
                for row in range(nrow):
                    rowlams.append( [] )
                for ilam in range(nlam-1):
                    row = ilam // ncol
                    rowlams[row].append( ilam )
                ymin=1.e+10
                ymax=-1.e+10
                xmin=1.e+10
                xmax=-1.e+10
                for ilam in range(nlam-1):
                    xmin = min(xmin,wk[ilam].xlo)
                    xmax = max(xmax,wk[ilam].xhi)
                    ymin = 0
                    ymax = max(max(ymax,max(wk[ilam].HistAB[0])),max(wk[ilam].HistBA[0]))
                    ymax = max(max(ymax,max(wk[ilam].GauAB)),max(wk[ilam].GauBA))

                for row in range(nrow):
                    oax = None
                    for ilam in rowlams[row]:
                        if oax is None:
                            ax  = plt.subplot(nrow,ncol,ilam+1)
                            ax.tick_params(axis='both',which='major',labelsize=9)
                            oax = ax
                        else:
                            ax  = plt.subplot(nrow,ncol,ilam+1,sharex=oax,sharey=oax)

                        x = [ 0.5*(wk[ilam].HistAB[1][ii+1]+wk[ilam].HistAB[1][ii]) for ii in range(len(wk[ilam].HistAB[0])) ]
                        y = wk[ilam].HistAB[0]
                        ax.scatter(x,y,c='k',marker='o',s=2)
                        x = wk[ilam].GauXs
                        y = wk[ilam].GauAB
                        ax.plot(x,y,c='k',linestyle='solid',linewidth=1)
                        
                        x = [ 0.5*(wk[ilam].HistBA[1][ii+1]+wk[ilam].HistBA[1][ii]) for ii in range(len(wk[ilam].HistBA[0])) ]
                        y = wk[ilam].HistBA[0]
                        ax.scatter(x,y,c='r',marker='o',s=2)
                        x = wk[ilam].GauXs
                        y = wk[ilam].GauBA
                        ax.plot(x,y,c='r',linestyle='solid',linewidth=1)


                        props = dict(boxstyle='round', facecolor='white', linewidth=0, alpha=0.5)
                        ax.text(0.5,0.975,r"%i $\leftrightarrow$ %i"%(ilam,ilam+1),transform=ax.transAxes,verticalalignment='top',horizontalalignment='center',fontsize=7)
                        ax.text(0.22,0.7,r"$\Pi_{ab}$=%.2f"%(wk[ilam].piab),transform=ax.transAxes,verticalalignment='bottom',horizontalalignment='center',fontsize=7)
                        ax.text(0.22,0.55,r"$\Pi_{ba}$=%.2f"%(wk[ilam].piba),transform=ax.transAxes,verticalalignment='bottom',horizontalalignment='center',fontsize=7)
                        ax.text(0.22,0.4,r"$K_{ab}$=%.2f"%(wk[ilam].Kab),transform=ax.transAxes,verticalalignment='bottom',horizontalalignment='center',fontsize=7)
                        ax.text(0.22,0.25,r"$K_{ba}$=%.2f"%(wk[ilam].Kba),transform=ax.transAxes,verticalalignment='bottom',horizontalalignment='center',fontsize=7)
                        ax.text(0.22,0.1,r"$S_{ba}$=%.1f"%(wk[ilam].PercentOverlap),transform=ax.transAxes,verticalalignment='bottom',horizontalalignment='center',fontsize=7)

                        
                        ax.set_xlim(xmin-0.7*(xmax-xmin),xmax+0.05*(xmax-xmin))
                        ax.set_ylim(ymin,ymax+0.1*(ymax-ymin))
                        if ilam%ncol != 0:
                            plt.setp(ax.get_yticklabels(),visible=False)
                        if ilam+ncol < nlam:
                            plt.setp(ax.get_xticklabels(),visible=False)
                        if ilam%ncol == 0:
                            ax.set_ylabel("Prob.",fontsize=9)
                        if ilam/ncol+1 == nrow:
                            ax.set_xlabel("$\Delta U$ (kcal/mol)",fontsize=9)
                        ax.grid(linestyle='dotted',linewidth=0.5)
                
                plt.tight_layout()
                fig.subplots_adjust(wspace=0, hspace=0)
                if len(trial.dataloc.filelabel) > 0:
                    histname = ".".join( splitall( trial.dataloc.inpdir ) ) + ".%s.hist.eps"%(trial.dataloc.filelabel)
                else:
                    histname = ".".join( splitall( trial.dataloc.inpdir ) ) + ".hist.eps"
                fname = os.path.join( directory, histname )
                plt.savefig(fname,dpi=300)
                fig.clf()



                ###################################
                
                if latex is not None:
                    latex.write(r"""
\clearpage
\pagebreak
\begin{figure*}
\includegraphics[clip,width=6in]{%s}\vspace{-0.3cm}
\includegraphics[clip,width=6in]{%s}\vspace{-0.3cm}
\caption{%s in %s}
\end{figure*}
"""%(gname,glname,
     trial.name.replace("_","\_"),
     trial.dataloc.inpdir.replace("_","\_")))
                    if os.path.exists(os.path.join(directory,sname)) and len(sname) > 0:
                        latex.write(r"""

\begin{figure*}
\includegraphics[clip,width=6in]{%s}\vspace{-0.3cm}
\caption{MBAR overlap matrix for %s in %s}
\end{figure*}

"""%(sname,
     trial.name.replace("_","\_"),
     trial.dataloc.inpdir.replace("_","\_")))
                    if len(dvdlname) > 0:
                        latex.write(r"""
\begin{figure*}
\includegraphics[clip,width=6in]{%s}\vspace{-0.3cm}
\caption{$\partial U/\partial\lambda$ time series for %s in %s}
\end{figure*}

"""%(dvdlname,
     trial.name.replace("_","\_"),
     trial.dataloc.inpdir.replace("_","\_")))

                    if "BAR" in self.get_methods() or "MBAR" in self.get_methods():
                        latex.write(r"""
\begin{figure*}
\includegraphics[clip,width=6in]{%s}\vspace{-0.3cm}
                        \caption{Wu and Kofke metrics for %s in %s. The black line is the energy distribution of Ub-Ua from ensemble of a, and the red line is Ua-Ub from the ensemble distribution of b. If the $\Pi$ metrics are negative, then more sampling is likely needed.}
\end{figure*}

"""%(histname,
     trial.name.replace("_","\_"),
     trial.dataloc.inpdir.replace("_","\_")))

                    


class TransformStage(object):
    def __init__(self,name,datalocs,methods=["TI","TI3","BAR","MBAR","IEXP","DEXP"]):
        self.name = name
        self.trials = []
        for dataloc in datalocs:
            self.trials.append( TransformTrial( name, dataloc, methods=methods ) )
            
    def get_methods(self):
        methods = []
        for i,trial in enumerate(self.trials):
            m = trial.get_methods()
            if i == 0:
                methods = m
            else:
                methods = [ x for x in m if x in methods ]
        return methods


    def analyze(self,tequil=0,tmax=1.e+10,plotdir=""):
        if len(self.trials) == 0:
            return
        
        for trial in self.trials:
            trial.analyze(tequil=tequil,tmax=tmax,plotdir=plotdir)

        methods = self.get_methods()
        self.full = DataAnalysis()
        self.forward=[]
        self.reverse=[]
        self.segment=[]

        
        for method in methods:
            self.full.ene[method] = GetAvgAndStdErr(
                [ trial.full.ene[method][0] for trial in self.trials ],
                [ trial.full.ene[method][1] for trial in self.trials ] )
        self.full.dvdl = None
        if trial.full.dvdl is not None:
            self.full.dvdl =  GetAvgAndStdErr(
                [ trial.full.dvdl[k][0] for trial in self.trials for k in range(len(trial.full.dvdl)) ],
                [ trial.full.dvdl[k][1] for trial in self.trials for k in range(len(trial.full.dvdl)) ] )
        self.full.maxtau = max( [ trial.full.maxtau for trial in self.trials ] )

        self.times = [0]*12
        for i in range(12):
            self.times[i] = sum( [ trial.times[i] for trial in self.trials ] ) / len(self.trials)
        
        
        for i in range(12):
            
            self.forward.append( DataAnalysis() )
            for method in methods:
                self.forward[i].ene[method] = GetAvgAndStdErr(
                    [ trial.forward[i].ene[method][0] for trial in self.trials ],
                    [ trial.forward[i].ene[method][1] for trial in self.trials ])
            self.forward[i].dvdl = None
            if trial.forward[i].dvdl is not None:
                self.forward[i].dvdl =  GetAvgAndStdErr(
                    [ trial.forward[i].dvdl[k][0] for trial in self.trials for k in range(len(trial.forward[i].dvdl)) ],
                    [ trial.forward[i].dvdl[k][1] for trial in self.trials for k in range(len(trial.forward[i].dvdl)) ])
            self.forward[i].maxtau = max( [ trial.forward[i].maxtau for trial in self.trials ] )

            self.reverse.append( DataAnalysis() )
            for method in methods:
                self.reverse[i].ene[method] = GetAvgAndStdErr(
                    [ trial.reverse[i].ene[method][0] for trial in self.trials ],
                    [ trial.reverse[i].ene[method][1] for trial in self.trials ])
            self.reverse[i].dvdl = None
            if trial.reverse[i].dvdl is not None:
                self.reverse[i].dvdl =  GetAvgAndStdErr(
                    [ trial.reverse[i].dvdl[k][0] for trial in self.trials for k in range(len(trial.reverse[i].dvdl)) ],
                    [ trial.reverse[i].dvdl[k][1] for trial in self.trials for k in range(len(trial.reverse[i].dvdl)) ])
            self.reverse[i].maxtau = max( [ trial.reverse[i].maxtau for trial in self.trials ] )

            self.segment.append( DataAnalysis() )
            for method in methods:
                self.segment[i].ene[method] = GetAvgAndStdErr(
                    [ trial.segment[i].ene[method][0] for trial in self.trials ],
                    [ trial.segment[i].ene[method][1] for trial in self.trials ])
            self.segment[i].dvdl = None
            if trial.segment[i].dvdl:
                self.segment[i].dvdl =  GetAvgAndStdErr(
                    [ trial.segment[i].dvdl[k][0] for trial in self.trials for k in range(len(trial.segment[i].dvdl)) ],
                    [ trial.segment[i].dvdl[k][1] for trial in self.trials for k in range(len(trial.segment[i].dvdl)) ])
            self.segment[i].maxtau = max( [ trial.segment[i].maxtau for trial in self.trials ] )

        
    def write_plot(self,directory,latex=None):
        if len(self.trials) > 1:
            pass
        for trial in self.trials:
            trial.write_plot(directory,latex=latex)
        
    # def get_energies(self,tequil=0,tmax=1.e+10):
    #     from collections import defaultdict as ddict
    #     ene
    #     return self.job.get_energies(tequil=tequil,tmax=tmax,dvdl_autocor=self.dvdl_autocor)



class Transform(object):
    def __init__(self,name,dict_of_datalocs_foreach_stage,methods=["TI","TI3","BAR","MBAR","IEXP","DEXP"]):
        from collections import defaultdict as ddict
        self.name = name
        self.stages = ddict(int)
        for stage in dict_of_datalocs_foreach_stage:
            datalocs = dict_of_datalocs_foreach_stage[stage]
            if datalocs is not None:
                self.stages[stage] = TransformStage( "%s.%s"%(self.name,stage), datalocs, methods=methods )
            
    def get_methods(self):
        methods = []
        for i,stage in enumerate(self.stages):
            m = self.stages[stage].get_methods()
            if i == 0:
                methods = m
            else:
                methods = [ x for x in m if x in methods ]
        return methods

    def analyze(self,tequil=0,tmax=1.e+10,plotdir=""):
        if len(self.stages) == 0:
            return
        
        for stage in self.stages:
            self.stages[stage].analyze(tequil=tequil,tmax=tmax,plotdir=plotdir)

        methods = self.get_methods()
        self.full = DataAnalysis()
        self.forward=[]
        self.reverse=[]
        self.segment=[]
        
        for method in methods:
            self.full.ene[method] = GetSumAndStdErr(
                [ self.stages[stage].full.ene[method][0] for stage in self.stages ],
                [ self.stages[stage].full.ene[method][1] for stage in self.stages ] )
        self.full.dvdl =  None
        self.full.maxtau = max( [ self.stages[stage].full.maxtau for stage in self.stages ] )

        self.times = [0]*12
        for i in range(12):
            self.times[i] = sum( [ self.stages[stage].times[i] for stage in self.stages ] ) / len(self.stages)
        
        for i in range(12):
            
            self.forward.append( DataAnalysis() )
            for method in methods:
                self.forward[i].ene[method] = GetSumAndStdErr(
                    [ self.stages[stage].forward[i].ene[method][0] for stage in self.stages ],
                    [ self.stages[stage].forward[i].ene[method][1] for stage in self.stages ])
            self.forward[i].dvdl =  None
            self.forward[i].maxtau = max( [ self.stages[stage].forward[i].maxtau for stage in self.stages ] )

            self.reverse.append( DataAnalysis() )
            for method in methods:
                self.reverse[i].ene[method] = GetSumAndStdErr(
                    [ self.stages[stage].reverse[i].ene[method][0] for stage in self.stages ],
                    [ self.stages[stage].reverse[i].ene[method][1] for stage in self.stages ])
            self.reverse[i].dvdl =  None
            self.reverse[i].maxtau = max( [ self.stages[stage].reverse[i].maxtau for stage in self.stages ] )

            self.segment.append( DataAnalysis() )
            for method in methods:
                self.segment[i].ene[method] = GetSumAndStdErr(
                    [ self.stages[stage].segment[i].ene[method][0] for stage in self.stages ],
                    [ self.stages[stage].segment[i].ene[method][1] for stage in self.stages ])
            self.segment[i].dvdl =  None
            self.segment[i].maxtau = max( [ self.stages[stage].segment[i].maxtau for stage in self.stages ] )

    def write_plot(self,directory,latex=None):
        import os
        from matplotlib import pyplot as plt
        nstages = len(self.stages)
        if nstages > 1:
            methods = self.get_methods()
            if "TI3" in methods:
                method="TI3"
            else:
                method="BAR"
            
            fig = plt.figure()
            w = (1) * (7./3)
            h = (1) * (7./3)
            fig.set_size_inches( w, h, forward=True )
            ax  = plt.subplot(1,1,1)
            plot_energy_timeseries(ax,self,method)
            ax.text(0.5,0.975,method,transform=ax.transAxes,
                    verticalalignment='top',horizontalalignment='center')
            ax.set_title(self.name)
            gname = "%s.net_GvsT.eps"%(self.name)
            fname = os.path.join( directory, gname )
            plt.tight_layout()
            plt.savefig(fname,dpi=300)
            fig.clf()

            fig = plt.figure()
            w = (nstages) * (7./3)
            h = (1) * (7./3)
            fig.set_size_inches( w, h, forward=True )
            for i,stage in enumerate(self.stages):
                ax  = plt.subplot(1,nstages,i+1)
                plot_energy_timeseries(ax,self.stages[stage],method)
                ax.text(0.5,0.975,method,transform=ax.transAxes,
                        verticalalignment='top',horizontalalignment='center')
                ax.set_title(self.stages[stage].name)
            sname = "%s.stages_GvsT.eps"%(self.name)
            fname = os.path.join( directory, sname )
            plt.tight_layout()
            plt.savefig(fname,dpi=300)
            fig.clf()

            if latex is not None:
                latex.write(r"""
\clearpage
\pagebreak
\begin{figure*}
\centering
\includegraphics[clip,width=2in]{%s}\vspace{-0.3cm}
\includegraphics[clip,width=6in]{%s}
\caption{%s}
\end{figure*}

"""%(gname,sname,self.name.replace("_","\_")))
            
        for stage in self.stages:
            self.stages[stage].write_plot(directory,latex=latex)
       
 
class RelativeTransform(object):
    def __init__(self,name,biotrans,reftrans):
        self.name = name
        self.stages = ddict(int)
        self.stages["bio"] = biotrans
        if reftrans is not None:
            self.stages["ref"] = reftrans
        self.analyze()
            
    def get_methods(self):
        methods = []
        for i,stage in enumerate(self.stages):
            m = self.stages[stage].get_methods()
            if i == 0:
                methods = m
            else:
                methods = [ x for x in m if x in methods ]
        return methods

    def analyze(self):
        methods = self.get_methods()
        self.full = DataAnalysis()
        self.forward=[]
        self.reverse=[]
        self.segment=[]
        signs = ddict(int)
        signs["bio"] =  1.
        signs["ref"] = -1.


        
        for method in methods:
            self.full.ene[method] = GetSumAndStdErr(
                [ signs[stage] * self.stages[stage].full.ene[method][0] for stage in self.stages ],
                [ self.stages[stage].full.ene[method][1] for stage in self.stages ] )
        self.full.dvdl =  None
        self.full.maxtau = max( [ self.stages[stage].full.maxtau for stage in self.stages ] )

        self.times = [ t for t in self.stages["bio"].times ]
        
        for i in range(12):
            
            self.forward.append( DataAnalysis() )
            for method in methods:
                self.forward[i].ene[method] = GetSumAndStdErr(
                    [ signs[stage] * self.stages[stage].forward[i].ene[method][0] for stage in self.stages ],
                    [ self.stages[stage].forward[i].ene[method][1] for stage in self.stages ])
            self.forward[i].dvdl =  None
            self.forward[i].maxtau = max( [ self.stages[stage].forward[i].maxtau for stage in self.stages ] )

            self.reverse.append( DataAnalysis() )
            for method in methods:
                self.reverse[i].ene[method] = GetSumAndStdErr(
                    [ signs[stage] * self.stages[stage].reverse[i].ene[method][0] for stage in self.stages ],
                    [ self.stages[stage].reverse[i].ene[method][1] for stage in self.stages ])
            self.reverse[i].dvdl =  None
            self.reverse[i].maxtau = max( [ self.stages[stage].reverse[i].maxtau for stage in self.stages ] )

            self.segment.append( DataAnalysis() )
            for method in methods:
                self.segment[i].ene[method] = GetSumAndStdErr(
                    [ signs[stage] * self.stages[stage].segment[i].ene[method][0] for stage in self.stages ],
                    [ self.stages[stage].segment[i].ene[method][1] for stage in self.stages ])
            self.segment[i].dvdl =  None
            self.segment[i].maxtau = max( [ self.stages[stage].segment[i].maxtau for stage in self.stages ] )


    def write_plot(self,directory,latex=None):
        import os
        from matplotlib import pyplot as plt
        nstages = len(self.stages)
        if nstages > 1:
            methods = self.get_methods()
            if "TI3" in methods:
                method="TI3"
            else:
                method="BAR"
            
            fig = plt.figure()
            w = (3) * (7./3)
            nstages = len(self.stages["bio"].stages)
            h = (nstages+1) * (7./3)
            fig.set_size_inches( w, h, forward=True )
            ax  = plt.subplot(nstages+1,3,1)
            plot_energy_timeseries(ax,self,method)
            ax.text(0.5,0.975,method,transform=ax.transAxes,
                    verticalalignment='top',horizontalalignment='center')
            ax.set_title(self.name)


            for imodel,model in enumerate(["bio","ref"]):
                icol = imodel+2
                ax  = plt.subplot(nstages+1,3,icol)
                t = self.stages[model]
                plot_energy_timeseries(ax,t,method)
                ax.text(0.5,0.975,method,transform=ax.transAxes,
                        verticalalignment='top',horizontalalignment='center')
                ax.set_title(t.name)
                for istage,stage in enumerate(sorted(t.stages)):
                    u = t.stages[stage]
                    ax  = plt.subplot(nstages+1,3,icol+(istage+1)*3)
                    plot_energy_timeseries(ax,u,method)
                    ax.text(0.5,0.975,method,transform=ax.transAxes,
                            verticalalignment='top',horizontalalignment='center')
                    ax.set_title(u.name)                
            gname = "%s.GvsT.eps"%(self.name)
            fname = os.path.join( directory, gname )
            plt.tight_layout()
            plt.savefig(fname,dpi=300)
            fig.clf()

            if latex is not None:
                latex.write(r"""
\clearpage
\pagebreak
\begin{figure*}
\centering
\includegraphics[clip,width=6in]{%s}
\caption{%s}
\end{figure*}

"""%(gname,self.name.replace("_","\_")))

        for model in self.stages:
            for stage in self.stages[model].stages:
                self.stages[model].stages[stage].write_plot(directory,latex=latex)











def collect_ene_vals(data_analysis,methods):
    vals = []
    valid_method = None
    for method in methods:
        if method in data_analysis.ene:
            vals.append( "%7.2f $\\pm$ %7.2f"%(data_analysis.ene[method][0],data_analysis.ene[method][1]) )
            valid_method = method
        else:
            vals.append( "%21s"%("na") )
    vals.append( "%7.0f"%( data_analysis.maxtau ) )
    return vals

            
def write_latex_header(latex):
    fh = open(latex,"w")
    fh.write(r"""
\makeatletter
\newcommand{\dontusepackage}[2][]{%
  \@namedef{ver@#2.sty}{9999/12/31}%
  \@namedef{opt@#2.sty}{#1}}
\makeatother

\dontusepackage{mciteplus}


% achemso loads mathptmx, but mathptmx messes \jmath, so save it here
\let\stdcoprod=\coprod
\let\stdamalg=\amalg
\let\stdjmath=\jmath

% \documentclass[journal=jctcce,manuscript=article,layout=twocolumn]{achemso}
\documentclass[journal=jctcce,manuscript=article,hyperref=false]{achemso}

% reset \jmath
\let\coprod=\stdcoprod
\let\amalg=\stdamalg
\let\jmath=\stdjmath



% Show all authors (don\'t use et al)
\makeatletter
\renewcommand*\acs@etal@firstonly{\acs@etal@truncatetrue}
\renewcommand*\acs@maxauthors{0}
\makeatother

\usepackage{amsmath}
\usepackage{graphicx}
\usepackage{multirow}
\usepackage{bm}
\usepackage{rotating}

% but we need to include txfonts to use \jmath now
\usepackage[T1]{fontenc}
\usepackage{txfonts}

\title{}
\begin{document}
\listoffigures
\listoftables



""")
    return fh



def write_latex_footer(fh):
    fh.write(r'\end{document}'+"\n\n")



    
def make_latex_document(odir,PI,methods=["TI","TI3","BAR","MBAR","IEXP","DEXP"]): #,tequil,tmax):
    #import tianalysis as tia
    import os

    #methods = ["TI","TI3","BAR","MBAR","IEXP","DEXP"]
    ncols = 2+len(methods)

    if not os.path.exists(odir):
        os.mkdir(odir)
        
    latex = write_latex_header( os.path.join(odir,"Main.tex") )
    
    summary = ddict( lambda: ddict( str ) )

    skips=[]
    
    for compound in sorted(PI):

        dloc = PI[compound]["bio"]
        if len( [ trial.tequil for step in dloc for trial in dloc[step] ] ) == 0:
            print(("%s bio is empty; skipping"%(compound)))
            skips.append(compound)
            continue
        tequil = min( [ trial.tequil for step in dloc for trial in dloc[step] ] )
        tmax   = max( [ trial.tmax  for step in dloc for trial in dloc[step] ] )
        biotrans = Transform( "%s.bio"%(compound), dloc, methods=methods )
        biotrans.analyze(tequil=tequil,tmax=tmax,plotdir=odir)
        if "ref" not in PI[compound]:
            PI[compound]["ref"] = None
        if PI[compound]["ref"] is None:
            reftrans = None
            reltrans = RelativeTransform( "%s.bio-ref"%(compound), biotrans, reftrans )
            #reltrans = biotrans
            #name = biotrans.name.replace("_","\_")
        else:
            dloc = PI[compound]["ref"]
            tequil = min( [ trial.tequil for step in dloc for trial in dloc[step] ] )
            tmax   = max( [ trial.tmax  for step in dloc for trial in dloc[step] ] )
            reftrans = Transform( "%s.ref"%(compound), dloc, methods=methods )
            reftrans.analyze(tequil=tequil,tmax=tmax,plotdir=odir)
            try:
                f = reftrans.full
            except:
                reftrans=None
            reltrans = RelativeTransform( "%s.bio-ref"%(compound), biotrans, reftrans )
        name = reltrans.name.replace("_","\_")
        
        latex.write(r"""
\clearpage
\pagebreak
\begin{table*}
\caption{%s}
{\small
\begin{tabular}{l %s}
\hline
"""%(name, " ".join( [ "r" for i in range(ncols-1) ] ) ) )
        
        latex.write( "%40s"%("Calculation") )
        for method in methods:
            latex.write(" & %18s"%(method))
        latex.write(" & max $\\tau$")
        latex.write("\\\\\n\\hline")

        for ifrac in [ 0, 3, 4, 6 ]:
            idx = 11-ifrac
            latex.write("\\multicolumn{%i}{c}{Excluding the first %i/12 of the simulations} \\\\\n"%(ncols,ifrac))
            s = "%40s & %s \\\\\n"%( name," & ".join( collect_ene_vals( reltrans.reverse[idx], methods ) ) )
            latex.write( s )
            summary[compound][ifrac] = s
            
            for model_name in reltrans.stages:
                model = reltrans.stages[model_name]
                latex.write( "%40s & %s \\\\\n"%( model.name.replace("_","\_")," & ".join( collect_ene_vals( model.reverse[idx], methods ) ) ) )
                if len(model.stages) > 1:
                    for stage_name in model.stages:
                        stage = model.stages[stage_name]
                        latex.write( "%40s & %s \\\\\n"%( stage.name.replace("_","\_")," & ".join( collect_ene_vals( stage.reverse[idx], methods ) ) ) )

        latex.write(r"""
\hline
\end{tabular}
}
\end{table*}
""")


        for model_name in reltrans.stages:
            model = reltrans.stages[model_name]
            if len(model.stages) > 0:
                for stage_name in model.stages:
                    stage = model.stages[stage_name]
                    if len(stage.trials) > 1 or len(model.stages) > 1:


                        latex.write(r"""

\begin{table*}
\caption{%s}
{\small
\begin{tabular}{l %s}
\hline
"""%(stage.name.replace("_","\_"), " ".join( [ "r" for i in range(ncols-1) ] ) ) )
        
                        latex.write( "%40s"%("Calculation") )
                        for method in methods:
                            latex.write(" & %18s"%(method))
                        latex.write(" & max $\\tau$")
                        latex.write("\\\\\n\\hline")

                        
                        for ifrac in [ 0, 3, 4, 6 ]:
                            idx = 11-ifrac
                            latex.write("\\multicolumn{%i}{c}{Excluding the first %i/12 of the simulations} \\\\\n"%(ncols,ifrac))
                            for trial in stage.trials:
                                latex.write( "%40s & %s \\\\\n"%( trial.name.replace("_","\_")," & ".join( collect_ene_vals( trial.reverse[idx], methods ) ) ) )

                        latex.write(r"""
\hline
\end{tabular}
}
\end{table*}
""")

        
        reltrans.write_plot(odir,latex)

    
        
    for ifrac in [ 0, 3, 4, 6 ]:
        idx = 11-ifrac

        if ifrac == 0:
            latex.write(r"""
\clearpage
\pagebreak
""")
            
        latex.write(r"""
\begin{table*}
\caption{Summary, excluding the first %i/12 of the simulations}
{\small
\begin{tabular}{l %s}
\hline
"""%(ifrac," ".join( [ "r" for i in range(ncols-1) ] ) ) )
        
        latex.write( "%40s"%("Calculation") )
        for method in methods:
            latex.write(" & %18s"%(method))
        latex.write(" & max $\\tau$")
        latex.write("\\\\\n\\hline")

        for compound in sorted(PI):
            if compound not in skips:
                s = "%40s & %s \\\\\n"%( reltrans.name.replace("_","\_")," & ".join( collect_ene_vals( reltrans.reverse[idx], methods ) ) )
                latex.write( summary[compound][ifrac] )
            
        latex.write(r"""
\hline
\end{tabular}
}
\end{table*}
""")

    write_latex_footer(latex)



    





                
    
if __name__ == "__main__":

    import argparse
    import fnmatch
    import os
    import glob
    import re
    import sys

    
    parser = argparse.ArgumentParser \
             ( formatter_class=argparse.RawDescriptionHelpFormatter,
               description="""
               """,
               epilog="""
               """)

    
    parser.add_argument \
        ("-i","--inpdir",
         help="directory containing dvdl and efep dat files",
         type=str,
         default="",
         required=False )

    
    parser.add_argument \
        ("-o","--outdir",
         help="directory to write dvdl and time-series plots",
         type=str,
         default="",
         required=False )

    
    parser.add_argument \
        ("-l","--label",
         help="only process dat files of the form dvdl_LABEL_tlam.dat and efep_LABEL_tlam_plam.dat",
         type=str,
         default="",
         required=False )

    
    parser.add_argument \
        ("-t","--tequil",
         help="equilibration time to exclude from analysis (ps)",
         type=int,
         default=0,
         required=False )

    
    parser.add_argument \
        ("-d","--dtime",
         help="delta-time used when generating the time-series (ps)",
         type=int,
         default=0,
         required=False )

    parser.add_argument \
        ("--no-ti",
         help="don't run TI calculation",
         action='store_true',
         required=False )

    parser.add_argument \
        ("--no-bar",
         help="don't run BAR calculation",
         action='store_true',
         required=False )

    parser.add_argument \
        ("--no-mbar",
         help="don't run MBAR calculation",
         action='store_true',
         required=False )

    parser.add_argument \
        ("--no-iexp",
         help="don't run IEXP calculation",
         action='store_true',
         required=False )

    parser.add_argument \
        ("--no-dexp",
         help="don't run DEXP calculation",
         action='store_true',
         required=False )


    

    
    args = parser.parse_args()

    dvdls = []
    efeps = []
    slabel = ""
    if len(args.label) > 0:
        slabel = "%s_"%(args.label)
    if len( args.inpdir ) == 0:
        for root, dirnames, filenames in os.walk('.'):
            for filename in fnmatch.filter(filenames, "dvdl_%s*.dat"%(slabel)):
                dvdls.append(os.path.join(root, filename))
            for filename in fnmatch.filter(filenames, "efep_%s*.dat"%(slabel)):
                efeps.append(os.path.join(root, filename))
    else:
        for f in glob.glob( os.path.join(args.inpdir,"dvdl_%s*.dat"%(slabel)) ):
            dvdls.append(f)
        for f in glob.glob( os.path.join(args.inpdir,"efep_%s*.dat"%(slabel)) ):
            efeps.append(f)

    dir_lab_lam = ddict( lambda: ddict( list ) )
    for f in dvdls:
        dname,fname = os.path.split(f)
        m = re.match( r"dvdl(.*)_([0-9\.]+?)\.dat", fname )
        if m:
            label = m.group(1)
            if len(label) > 0:
                if label[0] == "_":
                    label = label[1:]
            dir_lab_lam[dname][ label ].append( m.group(2) )
    for f in efeps:
        dname,fname = os.path.split(f)
        m = re.match( r"efep(.*)_([0-9\.]+?)_([0-9\.]+?)\.dat", fname )
        #m = re.match( r"efep_([0-9\.]+?)_([0-9\.]+?)\.dat", fname )
        if m:
            label = m.group(1)
            if len(label) > 0:
                if label[0] == "_":
                    label = label[1:]
            dir_lab_lam[dname][ label ].append( m.group(2) )
            dir_lab_lam[dname][ label ].append( m.group(3) )

    print ("DATA FILES LOCATED:\n")
    for d in sorted(dir_lab_lam):
        for lab in sorted(dir_lab_lam[d]):
            dir_lab_lam[d][lab] = [ float(lam) for lam in set(dir_lab_lam[d][lab]) ]
            dir_lab_lam[d][lab].sort()
            print(("directory:",d))
            print(("label:    ",lab))
            print(("lambdas:  ",dir_lab_lam[d][lab]))
            print ("")

    jobs = []
    methods=[]
    for d in sorted(dir_lab_lam):
        for lab in sorted(dir_lab_lam[d]):
            lams = dir_lab_lam[d][lab]
            job  = AlchemicalTransform(lams,directory=d,label=lab, methods=methods)

            dvdl_autocor=True
            if args.no_ti:
                job.ti=False
                job.ti3=False
                dvdl_autocor=False
            if args.no_bar:
                job.bar=False
            if args.no_mbar:
                job.mbar=False
            if args.no_iexp:
                job.iexp=False
            if args.no_dexp:
                job.dexp=False

            if not job.has_dvdl():
                print(("WARNING: directory='%s' label='%s'"%(d,lab)))
                print ("Some dvdl data is missing; skipping TI and TI3 analysis and dvdl) autocorrelations")
                job.ti=False
                job.ti3=False
                dvdl_autocor=False
                

            ene,dvdl,taus = job.get_energies(tequil=args.tequil,dvdl_autocor=dvdl_autocor)
            #taus     = job.get_autocorrelation_times(tequil=args.tequil,dvdl_autocor=dvdl_autocor)
            if len(jobs) == 0:
                methods = job.get_methods()
                print ("")
                job.write_header( sys.stdout )
            job.write_row(sys.stdout,ene,taus)

            fname = os.path.join( args.outdir, "tau_%s.dat"%(job.label) )
            fh = open(fname,"w")
            for k in range(job.nlam):
                lam = job.lams[k]
                fh.write("%.8f %19.3f\n"%(lam,taus[k]))
            fh.close()
            
            if job.ti or job.ti3:
                fname = os.path.join( args.outdir, "dvdl_%s.dat"%(job.label) )
                fh = open(fname,"w")
                for k in range(job.nlam):
                    lam = job.lams[k]
                    val,err = dvdl[k]
                    fh.write("%.8f %19.10e %19.10e\n"%(lam,val,1.96*err))
                fh.close()

            if job.bar:
                fname = os.path.join( args.outdir, "bar_%s.dat"%(job.label) )
                fh = open(fname,"w")
                for k in range(job.nlam):
                    lam = job.lams[k]
                    val = ene["BAR"][2][k]
                    fh.write("%.8f %19.10e\n"%(lam,val))
                fh.close()
                
            if job.mbar:
                fname = os.path.join( args.outdir, "mbar_%s.dat"%(job.label) )
                fh = open(fname,"w")
                for k in range(job.nlam):
                    lam = job.lams[k]
                    val = ene["MBAR"][2][k]
                    fh.write("%.8f %19.10e\n"%(lam,val))
                fh.close()
                
            if job.iexp:
                fname = os.path.join( args.outdir, "iexp_%s.dat"%(job.label) )
                fh = open(fname,"w")
                for k in range(job.nlam):
                    lam = job.lams[k]
                    val = ene["IEXP"][2][k]
                    fh.write("%.8f %19.10e\n"%(lam,val))
                fh.close()

            if job.dexp:
                fname = os.path.join( args.outdir, "dexp_%s.dat"%(job.label) )
                fh = open(fname,"w")
                for k in range(job.nlam):
                    lam = job.lams[k]
                    val = ene["DEXP"][2][k]
                    fh.write("%.8f %19.10e\n"%(lam,val))
                fh.close()
                
            jobs.append(job)


            

    if args.dtime > 0:
        for job in jobs:
            print ("")
            fwd = job.get_forward_timeseries(tequil=args.tequil,dtime=args.dtime,verbose=True)
            times = [t for t in sorted(fwd)]
            for method in methods:
                fname = os.path.join( args.outdir, "fwdGvsT_%s_%s.dat"%(job.label,method) )
                fh = open(fname,"w")
                for time in times:
                    val,err,x = fwd[time][method]
                    fh.write("%10.3f %15.8f %15.8f\n"%(time/1000.,val,1.96*err))
                fh.write("\n")
            fh.close()
            print ("")
            rev = job.get_reverse_timeseries(tequil=args.tequil,dtime=args.dtime,verbose=True)
            times = [t for t in sorted(rev)]
            for method in methods:
                fname = os.path.join( args.outdir, "revGvsT_%s_%s.dat"%(job.label,method) )
                fh = open(fname,"w")
                for time in times:
                    val,err,x = rev[time][method]
                    fh.write("%10.3f %15.8f %15.8f\n"%(time/1000.,val,1.96*err))
                fh.write("\n")
            fh.close()
            print ("")
            seg = job.get_segmented_timeseries(tequil=args.tequil,dtime=args.dtime,verbose=True)
            times = [t for t in sorted(seg)]
            for method in methods:
                fname = os.path.join( args.outdir, "segGvsT_%s_%s.dat"%(job.label,method) )
                fh = open(fname,"w")
                for time in times:
                    val,err,x = seg[time][method]
                    fh.write("%10.3f %15.8f %15.8f\n"%(time/1000.,val,1.96*err))
                fh.write("\n")
            fh.close()
            # xmgrace -settype xydy fwdGvsT_unlabeled_TI3.dat revGvsT_unlabeled_TI3.dat -settype xy segGvsT_unlabeled_TI3.dat -pexec 's0 errorbar size 0.5; s1 errorbar size 0.5; s0 errorbar riser linestyle 0; s1 errorbar riser linestyle 0; s2 line linestyle 0; s2 symbol 1; s2 symbol size 0.3333; s2 symbol fill 1; s0 line linewidth 2; s1 line linewidth 2'
            
    exit(0)
    
    nlam=12
    lams = [ float(i)/float(nlam-1.) for i in range(nlam) ]

    tequil = 10000
    dtime  = 10000
    inpdir = "/run/media/giese/homedir/piscti/protein/scti/mbar/data"
    label  = "sctest2"
    outdir = "./"
    
    t = AlchemicalTransform(lams,directory=inpdir,label=label, methods=methods)
    t.generic_run( tequil, dtime=dtime, outputdir=outdir )


    
