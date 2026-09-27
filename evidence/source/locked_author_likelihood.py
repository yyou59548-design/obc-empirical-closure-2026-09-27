import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import sys,pickle,yaml,re
import numpy as np
from scipy.linalg import solve_discrete_lyapunov, sqrtm
from scipy.special import ndtri,gammaln,betaln
from scipy.optimize import brentq

ROOT=Path('/mnt/data/obc_repro/OBC_frontier_run')
SCRIPTS=ROOT/'scripts'
if str(SCRIPTS) not in sys.path: sys.path.insert(0,str(SCRIPTS))
import model_engine as engine

META=ROOT/'reference/previous/raw/rank_spreads_exo_ztrend_SW_BAA_6420_ninit_0_meta.npz'


def _lhs_standard_normals(rng, size, dim):
    """Exact Chaospy 3.3.9 LHS+Normal construction in legacy np.random.RandomState order.

    Chaospy create_latin_hypercube_samples(order,dim): random(dim,order), then
    one permutation(order) per dimension; MvNormal.inv applies C @ ndtri(q).
    Returns array with shape (dim,)+size, before econsieve moveaxis/shuffle.
    """
    if isinstance(size,int): size_tuple=(size,)
    else: size_tuple=tuple(size)
    order=int(np.prod(size_tuple))
    u=rng.random_sample(order*dim).reshape(dim,order)
    for j in range(dim):
        perm=rng.permutation(order)
        u[j]=(perm+u[j])/order
    z=ndtri(u)
    return z.reshape((dim,)+size_tuple)


def _mvn_lhs(rng, mean, cov, size):
    mean=np.asarray(mean,float); cov=np.asarray(cov,float)
    dim=len(mean)
    try:
        C=np.linalg.cholesky(cov)
    except np.linalg.LinAlgError:
        C=np.real(sqrtm(cov))
    z=_lhs_standard_normals(rng,size,dim)
    flat=z.reshape(dim,-1)
    out=(C@flat + mean[:,None]).reshape(z.shape)
    # econsieve: np.moveaxis(res,0,res.ndim-1); np.random.shuffle(res)
    out=np.moveaxis(out,0,-1)
    rng.shuffle(out)
    return out


def _historical_logpdf(x,mean,cov,cond=1e-9):
    # exact grgrlib 2021 PSD logpdf convention
    s,u=np.linalg.eigh(cov)
    eps=cond*np.max(np.abs(s))
    d=s[s>eps]
    spinv=np.array([0.0 if abs(v)<=eps else 1.0/v for v in s])
    U=u*np.sqrt(spinv)[None,:]
    dev=np.asarray(x)-np.asarray(mean)
    maha=np.sum(np.square(dev@U),axis=-1)
    return float(-0.5*(len(d)*np.log(2*np.pi)+np.sum(np.log(d))+maha))


def _inv_gamma_spec(mu,sigma):
    def f(nu):
        return np.log(2*mu**2)-np.log((sigma**2+mu**2)*(nu-2))+2*(gammaln(nu/2)-gammaln((nu-1)/2))
    # brent is numerically equivalent for the density parameters at our tolerance
    nu=brentq(f,2+1e-12,1e12,xtol=1e-13)
    return (sigma**2+mu**2)*(nu-2),nu


class LockedAuthorLikelihood:
    """Author-2021 TEnKF likelihood with a full-horizon CRN bank sliced by endpoint.

    Full-sample behavior matches the archived 6420 protocol:
    N=350, archived filter_R, seed=1, Chaospy-3.3.9 Latin-hypercube Gaussian
    draws, nonlinear OBC transition, historical 2021 filter algebra.

    Endpoint calls slice one pre-generated 1964Q1--2019Q4 random bank, ensuring
    strict prefix identity. Only data endpoint changes.
    """
    def __init__(self,N=350,seed=1,l_max=4,k_max=16):
        with np.load(META,allow_pickle=True) as z:
            raw=z['yaml_raw'].item(); data=pickle.loads(z['data'].item())
            self.R=z['filter_R'].copy(); self.archive_mode=z['mode_x'].copy()
            self.archive_mode_f=float(z['mode_f']); self.archive_obs=[str(x) for x in z['obs']]
            self.archive_prior_names=[str(x) for x in z['prior_names']]
            self.archive_bounds=np.asarray(z['prior_bounds'],float).T
        self.raw=raw; self.data=data; self.N=int(N);self.seed=int(seed);self.l_max=int(l_max);self.k_max=int(k_max)
        # historical parser order: declaration shocks and declaration observables
        raw_clean=re.sub(r'@ ?\n',' ',raw.replace('^','**').replace(';',''))
        raw_clean=raw_clean.replace('\n ~ ','\n - ').replace('\n  ~ ','\n  - ').replace('   ~ ','   - ')
        y0=yaml.safe_load(raw_clean)
        self.shocks=list(y0['declarations']['shocks'])
        self.obs=list(y0['declarations']['observables'])
        yy,m=engine.parse_raw(raw,shock_order=self.shocks)
        # model_engine parser otherwise follows equation-dict order; restore historical declaration observable order
        # reuse the exact Variable class instantiated by parser
        V=type(m['observables'][0])
        m['observables']=[V(k) for k in self.obs]
        self.yy,self.m=yy,m
        self.funs,_=engine.build_lambdas(m,yy)
        self.pnames=list(yy['calibration']['parameters'])
        self.priors=yy['estimation']['prior']; self.prior_names=list(self.priors)
        assert self.prior_names==self.archive_prior_names
        self.ix=np.array([self.pnames.index(k) for k in self.prior_names],int)
        self.base=np.array([yy['calibration']['parameters'][k] for k in self.pnames],float)
        self.bounds=np.array([self.priors[k][1:3] for k in self.prior_names],float)
        self.Y=data[self.obs].to_numpy(dtype=float)
        self.periods=data.index.to_period('Q')
        self.endpoints=['2007Q4','2008Q4','2009Q4','2011Q4','2019Q4']
        self.ends={e:int(np.where(self.periods==e)[0][0])+1 for e in self.endpoints}
        # fixed full-horizon historical CRN bank, then endpoints are prefixes
        rng=np.random.RandomState(self.seed)
        T=len(self.Y)
        self.mus=_mvn_lhs(rng,np.zeros(len(self.obs)),self.R,(T,self.N))
        # eps covariance depends on parameters, so preserve the parameter-free standard normal LHS bank.
        # Historical MvNormal with diagonal Q is C @ z; generate z with same RNG calls/shuffles.
        # To preserve exact RNG path we generate using identity covariance here.
        self.eps_z=_mvn_lhs(rng,np.zeros(len(self.shocks)),np.eye(len(self.shocks)),(T,self.N))
        # Initial-state draw is parameter-dependent; preserve standard normal LHS bank and its historical RNG position.
        # Its dimension is fixed at 25 for this model.
        # Need dimension from a cheap baseline system.
        s0=engine.gen_system(self.m,self.yy,self.funs,self.base,l_max=self.l_max,k_max=self.k_max)
        self.nx=s0.dimq-s0.neps
        self.x_z=_mvn_lhs(rng,np.zeros(self.nx),np.eye(self.nx),self.N)
        self._ig={}
        for name,spec in self.priors.items():
            kind,mu,sd=spec[3:]
            if kind=='inv_gamma_dynare': self._ig[name]=_inv_gamma_spec(float(mu),float(sd))
    def full_from_estimated(self,x):
        p=self.base.copy();p[self.ix]=np.asarray(x,float);return p
    def estimated_from_full(self,p): return np.asarray(p,float)[self.ix]
    def logprior_estimated(self,x,enforce_bounds=False):
        x=np.asarray(x,float)
        if enforce_bounds and (np.any(x<self.bounds[:,0]) or np.any(x>self.bounds[:,1])): return -np.inf
        val=0.0
        for name,a in zip(self.prior_names,x):
            kind,mu,sd=self.priors[name][3:]; mu=float(mu);sd=float(sd)
            if kind=='normal': val += -.5*((a-mu)/sd)**2-np.log(sd)-.5*np.log(2*np.pi)
            elif kind=='beta':
                aa=(1-mu)*mu**2/sd**2-mu; bb=aa*(1/mu-1)
                if a<=0 or a>=1:return -np.inf
                val+=(aa-1)*np.log(a)+(bb-1)*np.log1p(-a)-betaln(aa,bb)
            elif kind=='gamma':
                if a<=0:return -np.inf
                scale=sd**2/mu;shape=mu/scale
                val+=(shape-1)*np.log(a)-a/scale-gammaln(shape)-shape*np.log(scale)
            elif kind=='inv_gamma_dynare':
                if a<=0:return -np.inf
                ss,nu=self._ig[name]
                val+=np.log(2)-gammaln(nu/2)-nu/2*(np.log(2)-np.log(ss))-(nu+1)*np.log(a)-.5*ss/a**2
            else: raise NotImplementedError(kind)
        return float(val)
    def _system(self,p): return engine.gen_system(self.m,self.yy,self.funs,p,l_max=self.l_max,k_max=self.k_max)
    @staticmethod
    def _root_for_mvn(P):
        try:return np.linalg.cholesky(P)
        except np.linalg.LinAlgError:return np.real(sqrtm(P))
    def evaluate_full(self,p,upto=None,return_path=False):
        p=np.asarray(p,float); s=self._system(p); neps=s.neps
        pm,qm,pt,qt,bm,bt=s.precalc
        # historical reduced-form initialization uses qmat[1,0]
        F=qm[1,0][:-neps,:-neps]; E=qm[1,0][:-neps,-neps:]
        sig=np.array([p[self.pnames.index(str(self.yy['calibration']['covariances'][k]))] for k in self.shocks],float)
        Qsh=np.diag(sig**2); Qstate=E@Qsh@E.T
        P=solve_discrete_lyapunov(F.T,Qstate)
        Cx=self._root_for_mvn(P)
        # x_z has shape N,nx after econsieve moveaxis/shuffle
        X=(self.x_z@Cx.T).T  # dimx,N
        # historical eps sample for diagonal Q: z @ chol(Q).T ; eps_z is T,N,neps
        eps=self.eps_z*sig[None,None,:]
        T=len(self.Y) if upto is None else (self.ends[upto] if isinstance(upto,str) else int(upto))
        I1=np.ones(self.N); I2=np.eye(self.N)-np.outer(I1,I1)/self.N
        Yens=np.empty((len(self.obs),self.N))
        ll=np.empty(T); flags=np.empty(T,dtype=int)
        forecasts=np.empty((T,len(self.obs)))
        # local transition loop, matrices follow historical t_func_jit
        for t,z in enumerate(self.Y[:T]):
            flagcnt=0
            for i in range(self.N):
                state=X[:,i]; shock=eps[t,i]
                qin=np.hstack((state,shock))
                l,k,flag=engine.find_lk(bm,bt,s.x_bar,qin)
                flagcnt += int(flag!=0)
                pp=pm[l,k]@qin+pt[l,k]
                qq=qm[l,k][:-neps]@qin+qt[l,k][:-neps]
                X[:,i]=qq
                Yens[:,i]=s.hx[0]@pp+s.hx[1]@qq+s.hx[2]
            Xbar=X@I2; Ybar=Yens@I2; ZZ=np.outer(z,I1)
            S=np.cov(Yens)+self.R
            forecasts[t]=Yens.mean(axis=1)
            ll[t]=_historical_logpdf(z-forecasts[t],np.zeros(len(self.obs)),S)
            # exact historical algebra
            X += Xbar@Ybar.T@np.linalg.inv((self.N-1)*S)@(ZZ-Yens-self.mus[t].T)
            flags[t]=flagcnt
        if return_path:return np.cumsum(ll),flags,forecasts,P
        return float(ll.sum())
    def logposterior_estimated(self,x,upto=None,enforce_bounds=False):
        lp=self.logprior_estimated(x,enforce_bounds=enforce_bounds)
        if not np.isfinite(lp):return -np.inf
        return lp+self.evaluate_full(self.full_from_estimated(x),upto=upto)

if __name__=='__main__':
    import json,time
    L=LockedAuthorLikelihood()
    x=L.archive_mode
    t=time.time(); lp=L.logprior_estimated(x); ll=L.evaluate_full(L.full_from_estimated(x)); dt=time.time()-t
    print(json.dumps({'archive_mode_f':L.archive_mode_f,'logprior':lp,'loglikelihood':ll,'reconstructed_logposterior':lp+ll,'difference':lp+ll-L.archive_mode_f,'seconds':dt},indent=2))

# Fast exact-algebra implementation: numba only accelerates transition evaluation.
try:
    from numba import njit
    @njit(cache=True)
    def _select_lk_fast(bm,bt,xbar,q):
        lm=bt.shape[1];km=bt.shape[2];l=0;k=0;flag=0
        while np.dot(bm[2,l,0],q)+bt[2,l,0]>xbar:
            l+=1
            if l==lm:break
        if l<lm:
            found=False
            for ll in range(lm):
                for kk in range(1,km):
                    if ll and np.dot(bm[0,ll,kk],q)+bt[0,ll,kk]<xbar:continue
                    if ll>1 and np.dot(bm[1,ll,kk],q)+bt[1,ll,kk]<xbar:continue
                    if np.dot(bm[4,ll,kk],q)+bt[4,ll,kk]<xbar:continue
                    if np.dot(bm[2,ll,kk],q)+bt[2,ll,kk]>xbar:continue
                    if kk>1 and np.dot(bm[3,ll,kk],q)+bt[3,ll,kk]>xbar:continue
                    l=ll;k=kk;found=True;break
                if found:break
            if not found:
                flag=1;l=0;k=0;loop=False
                while True:
                    k+=1
                    if not loop:loop=np.dot(bm[4,0,k],q)+bt[4,0,k]<xbar
                    elif np.dot(bm[4,0,k],q)+bt[4,0,k]>xbar:break
                    if k==km-1:flag=2;break
        else:l=0
        return l,k,flag

    @njit(cache=True)
    def _transition_fast(X,E,pm,qm,pt,qt,bm,bt,Hp,Hq,hc,xbar,neps):
        N=X.shape[0]; nstate=X.shape[1]; nobs=len(hc)
        Xn=np.empty_like(X); Z=np.empty((N,nobs)); flags=np.zeros(N,np.int64)
        for i in range(N):
            qin=np.empty(nstate+neps)
            qin[:nstate]=X[i]; qin[nstate:]=E[i]
            l,k,flag=_select_lk_fast(bm,bt,xbar,qin)
            pp=pm[l,k]@qin+pt[l,k]
            qq=qm[l,k][:-neps]@qin+qt[l,k][:-neps]
            Xn[i]=qq; Z[i]=Hp@pp+Hq@qq+hc; flags[i]=flag
        return Xn,Z,flags
except Exception:
    _transition_fast=None


def _evaluate_full_fast(self,p,upto=None,return_path=False):
    p=np.asarray(p,float); s=self._system(p); neps=s.neps
    pm,qm,pt,qt,bm,bt=s.precalc
    F=qm[1,0][:-neps,:-neps]; E=qm[1,0][:-neps,-neps:]
    sig=np.array([p[self.pnames.index(str(self.yy['calibration']['covariances'][k]))] for k in self.shocks],float)
    P=solve_discrete_lyapunov(F.T,E@np.diag(sig**2)@E.T)
    Cx=self._root_for_mvn(P)
    X=self.x_z@Cx.T  # N, dimx
    eps=self.eps_z*sig[None,None,:]
    T=len(self.Y) if upto is None else (self.ends[upto] if isinstance(upto,str) else int(upto))
    ll=np.empty(T);flags=np.empty(T,dtype=int);forecasts=np.empty((T,len(self.obs)))
    I1=np.ones(self.N); I2=np.eye(self.N)-np.outer(I1,I1)/self.N
    for t,z in enumerate(self.Y[:T]):
        X,Z,fl=_transition_fast(X,eps[t],pm,qm,pt,qt,bm,bt,s.hx[0],s.hx[1],s.hx[2],s.x_bar,neps)
        # Preserve the historical matrix-centering algebra literally.  This
        # matters at OBC branch boundaries: algebraically equivalent mean
        # subtraction can change roundoff enough to select a different branch.
        Xt=X.T; Yt=Z.T
        Xbar=Xt@I2; Ybar=Yt@I2; ZZ=np.outer(z,I1)
        S=np.cov(Yt)+self.R
        zm=np.mean(Yt,axis=1)
        ll[t]=_historical_logpdf(z-zm,np.zeros(len(self.obs)),S); forecasts[t]=zm
        Xt += Xbar@Ybar.T@np.linalg.inv((self.N-1)*S)@(ZZ-Yt-self.mus[t].T)
        X=Xt.T
        flags[t]=int(np.count_nonzero(fl))
    if return_path:return np.cumsum(ll),flags,forecasts,P
    return float(ll.sum())

LockedAuthorLikelihood.evaluate_full_fast=_evaluate_full_fast
LockedAuthorLikelihood.evaluate_full=_evaluate_full_fast

# Optional whole-filter JIT retaining the historical mixed centering conventions.
try:
    @njit(cache=True)
    def _logpdf_psd_numba(dev,S,cond=1e-9):
        s,u=np.linalg.eigh(S)
        mx=0.0
        for i in range(len(s)):
            a=abs(s[i])
            if a>mx: mx=a
        eps=cond*mx; k=0; logdet=0.0; maha=0.0
        for j in range(len(s)):
            if s[j]>eps:
                k+=1; logdet+=np.log(s[j])
                proj=0.0
                for i in range(len(dev)): proj += dev[i]*u[i,j]
                maha += proj*proj/s[j]
        return -0.5*(k*np.log(2*np.pi)+logdet+maha)

    @njit(cache=True)
    def _filter_loop_numba(Yobs,X,eps,mus,R,pm,qm,pt,qt,bm,bt,Hp,Hq,hc,xbar,neps,I2):
        T=Yobs.shape[0]; N=X.shape[0]; nobs=Yobs.shape[1]; nx=X.shape[1]
        ll=np.empty(T); flags=np.empty(T,np.int64); forecasts=np.empty((T,nobs))
        I1=np.ones(N)
        for t in range(T):
            X,Z,fl=_transition_fast(X,eps[t],pm,qm,pt,qt,bm,bt,Hp,Hq,hc,xbar,neps)
            Xt=X.T; Yt=Z.T
            Xbar=Xt@I2; Ybar=Yt@I2
            # np.cov(Yt) convention: subtract arithmetic row means, divide by N-1.
            Sm=np.empty((nobs,nobs))
            means=np.empty(nobs)
            for a in range(nobs):
                ss=0.0
                for i in range(N): ss+=Yt[a,i]
                means[a]=ss/N
            for a in range(nobs):
                for b in range(nobs):
                    ss=0.0
                    for i in range(N): ss+=(Yt[a,i]-means[a])*(Yt[b,i]-means[b])
                    Sm[a,b]=ss/(N-1)+R[a,b]
            for a in range(nobs): forecasts[t,a]=means[a]
            dev=np.empty(nobs)
            for a in range(nobs): dev[a]=Yobs[t,a]-means[a]
            ll[t]=_logpdf_psd_numba(dev,Sm)
            # A = inv((N-1)*S); update = Xbar Ybar' A (ZZ-Yt-mus')
            A=np.linalg.inv((N-1)*Sm)
            G=Xbar@Ybar.T@A
            Innov=np.empty((nobs,N))
            for a in range(nobs):
                for i in range(N): Innov[a,i]=Yobs[t,a]-Yt[a,i]-mus[t,i,a]
            Xt = Xt + G@Innov
            X=Xt.T
            cnt=0
            for i in range(N):
                if fl[i]!=0: cnt+=1
            flags[t]=cnt
        return ll,flags,forecasts,X

    def evaluate_full_jit(self,p,upto=None,return_path=False):
        p=np.asarray(p,float); s=self._system(p); neps=s.neps
        pm,qm,pt,qt,bm,bt=s.precalc
        F=qm[1,0][:-neps,:-neps]; E=qm[1,0][:-neps,-neps:]
        sig=np.array([p[self.pnames.index(str(self.yy['calibration']['covariances'][k]))] for k in self.shocks],float)
        P=solve_discrete_lyapunov(F.T,E@np.diag(sig**2)@E.T)
        Cx=self._root_for_mvn(P); X=self.x_z@Cx.T; eps=self.eps_z*sig[None,None,:]
        T=len(self.Y) if upto is None else (self.ends[upto] if isinstance(upto,str) else int(upto))
        I2=np.eye(self.N)-np.ones((self.N,self.N))/self.N
        ll,flags,forecasts,_=_filter_loop_numba(self.Y[:T],X,eps[:T],self.mus[:T],self.R,pm,qm,pt,qt,bm,bt,s.hx[0],s.hx[1],s.hx[2],s.x_bar,neps,I2)
        if return_path:return np.cumsum(ll),flags,forecasts,P
        return float(ll.sum())
    LockedAuthorLikelihood.evaluate_full_jit=evaluate_full_jit
except Exception:
    pass
