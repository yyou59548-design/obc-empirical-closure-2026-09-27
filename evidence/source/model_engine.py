import sys, types, numpy as np, pandas as pd, yaml, sympy, scipy.linalg as sl, itertools, pickle
from pathlib import Path
from sympy.printing.str import StrPrinter
from sympy.core.cache import clear_cache

StrPrinter._print_TSymbol = lambda self, x: x.__str__()
class Parameter(sympy.Symbol):
    def __init__(self,name,exp_date=0): super().__init__(); self.name=name
    def __repr__(self): return self.name
class TSymbol(sympy.Symbol):
    def __init__(self,name,**args):
        super().__init__(); date=args.get('date',0); exp=args.get('exp_date',0)
        self._assumptions_orig['date']=date; self._assumptions0 += ('date',date),
        self._assumptions_orig['exp_date']=exp; self._assumptions0 += ('exp_date',exp),
        self._mhash=None; self.__hash__()
    def __call__(self,lead):
        newdate=int(self.date)+int(lead); newname=str(self.name); clear_cache(); return self.__class__(newname,date=newdate)
    @property
    def date(self): return self.assumptions0['date']
    @property
    def exp_date(self): return self.assumptions0['exp_date']
    def _hashable_content(self): return (self.name,str(self.date),str(self.exp_date))
    def __getstate__(self): return {}
    def class_key(self): return (2,0,self.name,self.date)
    @property
    def lag(self): return self.date
    def __str__(self): return self.name if self.lag==0 else self.name+f'({self.lag})'
class Variable(TSymbol):
    @property
    def fortind(self): return 'v_'+self.name if self.date<=0 else 'v_E'+self.name
    def __str__(self): return super().__str__() if self.exp_date==0 else 'E['+str(self.exp_date)+']'+super().__str__()
    def __repr__(self): return self.__str__()
    __sstr__=__str__
class Shock(TSymbol):
    @property
    def fortind(self): return 'e_'+self.name if self.date<=0 else 'e_E'+self.name
class Equation(sympy.Equality):
    def __new__(cls,lhs,rhs,name=None): return super(sympy.Equality,cls).__new__(cls,lhs,rhs)
    @property
    def set_eq_zero(self): return self.lhs-self.rhs
    @property
    def variables(self): return [v for v in self.atoms() if isinstance(v,Variable)]

def parse_raw(raw, shock_order=None):
    mtxt=raw.replace('^','**').replace(';','')
    # mimic parser cleanup
    import re
    mtxt=re.sub(r'@ ?\n',' ',mtxt)
    mtxt=mtxt.replace('\n ~ ','\n - ').replace('\n  ~ ','\n  - ').replace('   ~ ','   - ')
    yy=yaml.safe_load(mtxt)
    dec=yy['declarations']; cal=yy['calibration']
    var_ordering=[Variable(v) for v in dec['variables']]
    par_ordering=[Parameter(v) for v in cal['parameters']]
    if shock_order is None: shock_order=list(cal['covariances'])
    shk_ordering=[Shock(v) for v in shock_order]
    other=[Parameter(v) for v in cal.get('parafunc',{})]
    ctx={s.name:s for s in var_ordering+par_ordering+shk_ordering+other}
    obs_equations={k:v for k,v in yy['equations']['observables'].items()}
    # constraint
    raw_const=yy['equations']['constraint'][0]
    lhs,rhs=raw_const.split('=')
    c_var=Variable(lhs.strip())
    lhs=eval(lhs,{'__builtins__':{}},ctx); rhs=eval(rhs,{'__builtins__':{}},ctx)
    if isinstance(lhs,(int,float)): lhs=sympy.sympify(lhs)
    if isinstance(rhs,(int,float)): rhs=sympy.sympify(rhs)
    const_eq=Equation(lhs,rhs)
    # model equations
    equations=[]
    raw_equations=yy['equations']['model']
    for eq in raw_equations:
        if '=' in eq: lhs,rhs=eq.split('=',1)
        else: lhs,rhs=eq,'0'
        lhs=eval(lhs,{'__builtins__':{}},ctx); rhs=eval(rhs,{'__builtins__':{}},ctx)
        if isinstance(lhs,(int,float)): lhs=sympy.sympify(lhs)
        if isinstance(rhs,(int,float)): rhs=sympy.sympify(rhs)
        equations.append(Equation(lhs,rhs))
    # arbitrary lags of exo shocks exactly as current parser
    it=itertools.chain.from_iterable
    all_shocks=[list(eq.atoms(Shock)) for eq in equations]
    max_lag_exo={}
    for s in shk_ordering:
        dates=[i.date for i in it(all_shocks) if i.name==s.name]
        max_lag_exo[s]=min(dates)
        all_shocks=[list(eq.atoms(Shock)) for eq in equations]  # reset because chain exhausted issue workaround
    # Above source creates a fresh chain comprehension each s in actual list comp. recompute robustly
    max_lag_exo={s:min([i.date for eq in equations for i in eq.atoms(Shock) if i.name==s.name]) for s in shk_ordering}
    for s in shk_ordering:
        if abs(max_lag_exo[s])>0:
            var_s=Variable(s.name+'_VAR'); var_ordering.append(var_s); equations.append(Equation(var_s,s))
            subs1=[s(-i) for i in np.arange(1,abs(max_lag_exo[s])+1)]
            subs2=[var_s(-i) for i in np.arange(1,abs(max_lag_exo[s])+1)]
            equations=[eq.subs(dict(zip(subs1,subs2))) for eq in equations]
    # endo >1 lags/leads
    all_vars=[list(eq.atoms(Variable)) for eq in equations]
    max_lead={v:max([i.date for eq in equations for i in eq.atoms(Variable) if i.name==v.name]) for v in var_ordering}
    max_lag={v:min([i.date for eq in equations for i in eq.atoms(Variable) if i.name==v.name]) for v in var_ordering}
    subs_dict={}; old_var=var_ordering[:]
    for v in old_var:
        for i in np.arange(2,abs(max_lag[v])+1):
            var_l=Variable(v.name+'_LAG'+str(i-1)); var_l_1=v(-1) if i==2 else Variable(v.name+'_LAG'+str(i-2),date=-1)
            subs_dict[Variable(v.name,date=-i)]=var_l(-1); var_ordering.append(var_l); equations.append(Equation(var_l,var_l_1))
        for i in np.arange(2,abs(max_lead[v])+1):
            var_l=Variable(v.name+'_LEAD'+str(i-1)); var_l_1=v(+1)
            subs_dict[Variable(v.name,date=+i)]=var_l(+1); var_ordering.append(var_l); equations.append(Equation(var_l,var_l_1))
    equations=[eq.subs(subs_dict) for eq in equations]
    # evaluate observables with current context (sum/range available)
    ctx2=dict(ctx);ctx2.update({'sum':np.sum,'range':range})
    obs_sym={k:eval(v,{'__builtins__':{}},ctx2) for k,v in obs_equations.items()}
    return yy, {'var_ordering':var_ordering,'par_ordering':par_ordering,'shk_ordering':shk_ordering,'other_para':other,'para_func':cal['parafunc'],'const_var':c_var,'const_eq':const_eq,'perturb_eq':equations,'observables':[Variable(k) for k in obs_equations], 'obs_equations':obs_sym}

def build_lambdas(m,yy):
    params=m['par_ordering']; other=m['other_para']; sub_var=m['var_ordering']; slist=m['shk_ordering']
    fvarl=[v(+1) for v in sub_var];lvarl=[v(-1) for v in sub_var]
    subs={};subs.update({v:0 for v in sub_var});subs.update({v(1):0 for v in sub_var});subs.update({v(-1):0 for v in sub_var})
    no_var=len(sub_var);evar=len(slist);ovar=len(m['observables'])
    bb=sympy.zeros(1,no_var+no_var);bb_PSI=sympy.zeros(1,evar)
    AA=sympy.zeros(no_var-1,no_var);BB=sympy.zeros(no_var-1,no_var);CC=sympy.zeros(no_var-1,no_var);PSI=sympy.zeros(no_var-1,evar)
    bb_var=[x for x in m['const_eq'].atoms(Variable) if x.date<=0]; full_var=sub_var+lvarl
    for v in bb_var: bb[full_var.index(v)] = -m['const_eq'].set_eq_zero.diff(v).subs(subs)
    for s in m['const_eq'].atoms(Shock): bb_PSI[slist.index(s)] = -m['const_eq'].set_eq_zero.diff(s).subs(subs)
    for i,eq in enumerate(m['perturb_eq']):
        for v in [x for x in eq.atoms(Variable) if x.date>0]: AA[i,fvarl.index(v)]=eq.set_eq_zero.diff(v).subs(subs)
        for v in [x for x in eq.atoms(Variable) if x.date==0]: BB[i,sub_var.index(v)]=eq.set_eq_zero.diff(v).subs(subs)
        for v in [x for x in eq.atoms(Variable) if x.date<0]: CC[i,lvarl.index(v)]=eq.set_eq_zero.diff(v).subs(subs)
        for s in eq.atoms(Shock): PSI[i,slist.index(s)] = -eq.set_eq_zero.diff(s).subs(subs)
    ZZ0=sympy.zeros(ovar,no_var);ZZ1=sympy.zeros(ovar,1);vlist=sub_var
    for i,obs in enumerate(m['observables']):
        eq=m['obs_equations'][str(obs)]; ZZ1[i,0]=eq.subs(subs)
        for v in [x for x in eq.atoms(Variable) if x.date>=0]: ZZ0[i,vlist.index(v)] = eq.diff(v).subs(subs)
    allsyms=params+other
    funs={k:sympy.lambdify([allsyms],obj,modules='numpy') for k,obj in {'AA':AA,'BB':BB,'CC':CC,'PSI':PSI,'bb':bb,'bb_PSI':bb_PSI,'ZZ0':ZZ0,'ZZ1':ZZ1}.items()}
    return funs, (AA,BB,CC,PSI,bb,bb_PSI,ZZ0,ZZ1)

def parafunc_values(yy,p):
    pnames=list(yy['calibration']['parameters']);ctx={n:float(v) for n,v in zip(pnames,p)};exprs=yy['calibration']['parafunc'];vals={};pending=dict(exprs)
    safe={'exp':np.exp,'log':np.log,'sqrt':np.sqrt,'abs':abs,'np':np}
    for _ in range(len(pending)+5):
        prog=False
        for k,e in list(pending.items()):
            try: vals[k]=float(eval(str(e).replace('^','**'),{'__builtins__':{}},{**safe,**ctx,**vals}));del pending[k];prog=True
            except NameError: pass
        if not pending:break
        if not prog:raise RuntimeError(pending)
    return [vals[k] for k in exprs],vals

def fast0(A,mode=-1,tol=1e-8):
    con=np.abs(A)<tol
    return con if mode==-1 else con.all(axis=mode)
def ouc(x,y):
    out=np.empty_like(x,dtype=bool);nz=y!=0;out[~nz]=True;out[nz]=np.abs(x[nz]/y[nz])>1.;return out
def klein(A,B,nstates,force=False):
    SS,TT,alp,bet,Q,Z=sl.ordqz(A,B,sort='ouc');out=ouc(alp,bet)
    if nstates!=out.sum() and not force:raise ValueError(f'BK {nstates} != {out.sum()}')
    S11=SS[:nstates,:nstates];T11=TT[:nstates,:nstates];Z11=Z[:nstates,:nstates];Z21=Z[nstates:,:nstates]
    return np.real_if_close(Z21@np.linalg.inv(Z11)).astype(float), np.real_if_close(Z11@np.linalg.inv(S11)@T11@np.linalg.inv(Z11)).astype(float)

def preprocess(S,T,V,W,h,fq1,fp1,fq0,omg,lam,x_bar,l_max,k_max):
    S=np.array(S,float,copy=True);T=np.array(T,float,copy=True);V=np.array(V,float,copy=True);W=np.array(W,float,copy=True);h=np.array(h,float,copy=True)
    dimp,dimq=omg.shape;lm=l_max+1;km=k_max+1
    T22i=np.linalg.inv(T[dimq:,dimq:]);T[dimq:]=T22i@T[dimq:];S[dimq:]=T22i@S[dimq:]
    W22i=np.linalg.inv(W[dimq:,dimq:]);W[dimq:]=W22i@W[dimq:];V[dimq:]=W22i@V[dimq:];h[dimq:]=W22i@h[dimq:]
    pmat=np.empty((lm,km,dimp,dimq));qmat=np.empty((lm,km,dimq,dimq));pterm=np.empty((lm,km,dimp));qterm=np.empty((lm,km,dimq));pmat[0,0]=omg;pterm[0,0]=0;qmat[0,0]=lam;qterm[0,0]=0
    for l in range(lm):
      for k in range(km):
       if k or l:
        ll=max(l-1,0);kk=k if l else max(k-1,0);A=S if l else V;B=T if l else W;c=np.zeros(dimq) if l else h[:dimq]
        inv=np.linalg.inv(A[:dimq,:dimq]+A[:dimq,dimq:]@pmat[ll,kk]);qmat[l,k]=inv@B[:dimq,:dimq];qterm[l,k]=inv@(c-A[:dimq,dimq:]@pterm[ll,kk])
        cc=np.zeros(dimp) if l else h[dimq:];dum=A[dimq:,:dimq]+A[dimq:,dimq:]@pmat[ll,kk];pterm[l,k]=dum@qterm[l,k]+A[dimq:,dimq:]@pterm[ll,kk]-cc;pmat[l,k]=dum@qmat[l,k]-B[dimq:,:dimq]
    bmat=np.empty((5,lm,km,dimq));bterm=np.empty((5,lm,km))
    for l in range(lm):
      for k in range(km):
       lm0=np.eye(dimq);xi=np.zeros(dimq)
       for ss in range(l+k+1):
        ll=max(l-ss,0);kk=max(min(k,k+l-ss),0);y2r=fp1@pmat[ll,kk]+fq1@qmat[ll,kk]+fq0;cr=fp1@pterm[ll,kk]+fq1@qterm[ll,kk]
        if ss==0:bmat[0,l,k]=y2r@lm0;bterm[0,l,k]=cr+y2r@xi
        if ss==l-1:bmat[1,l,k]=y2r@lm0;bterm[1,l,k]=cr+y2r@xi
        if ss==l:bmat[2,l,k]=y2r@lm0;bterm[2,l,k]=cr+y2r@xi
        if ss==l+k-1:bmat[3,l,k]=y2r@lm0;bterm[3,l,k]=cr+y2r@xi
        if ss==l+k:bmat[4,l,k]=y2r@lm0;bterm[4,l,k]=cr+y2r@xi
        lm0=qmat[ll,kk]@lm0;xi=qmat[ll,kk]@xi+qterm[ll,kk]
    return pmat,qmat,pterm,qterm,bmat,bterm

def gen_system(m,yy,funs,p,l_max=3,k_max=17,force=False):
    class S:pass
    s=S();s.shocks=[str(x) for x in m['shk_ordering']];s.neps=len(s.shocks);s.const_var=m['const_var'];ov,vals=parafunc_values(yy,p);ppar=list(map(float,p))+ov;s.x_bar=vals['x_bar']
    AA0=np.array(funs['AA'](ppar),float);BB0=np.array(funs['BB'](ppar),float);CC0=np.array(funs['CC'](ppar),float);DD0=-np.array(funs['PSI'](ppar),float);fbc=np.array(funs['bb'](ppar),float).flatten();fd0=-np.array(funs['bb_PSI'](ppar),float).flatten();vv0=np.array([v.name for v in m['var_ordering']],object)
    fb0=-fbc[:len(vv0)];fc0=-fbc[len(vv0):];ZZ0=np.array(funs['ZZ0'](ppar),float);ZZ1=np.array(funs['ZZ1'](ppar),float).squeeze()
    # exact normalization from gen_sys
    c_arg=list(vv0).index(str(s.const_var));fc0=-fc0/fb0[c_arg];fb0=-fb0/fb0[c_arg]
    inall=(~fast0(AA0,0)) & (~fast0(CC0,0))
    if np.any(inall):
      n=int(inall.sum());vv0=np.hstack((vv0,[v+'_lag' for v in vv0[inall]]));AA0=np.pad(AA0,((0,n),(0,n)));BB0=np.pad(BB0,((0,n),(0,n)));CC0=np.pad(CC0,((0,n),(0,n)));DD0=np.pad(DD0,((0,n),(0,0)));fb0=np.pad(fb0,(0,n));fc0=np.pad(fc0,(0,n));ZZ0=np.pad(ZZ0,((0,0),(0,n)));BB0[-n:,-n:]=np.eye(n);BB0[-n:,:-n][:,inall]=-np.eye(n);CC0[:,-n:]=CC0[:,:-n][:,inall];CC0[:,:-n][:,inall]=0
    s.primitive_A=AA0.copy(); s.primitive_B=BB0.copy(); s.primitive_C=CC0.copy(); s.primitive_D=DD0.copy(); s.primitive_names=vv0.copy()
    de=len(s.shocks);AA0=np.pad(AA0,((0,de),(0,de)));BB0=sl.block_diag(BB0,np.eye(de));CC0=np.block([[CC0,DD0],[np.zeros((de,AA0.shape[1]))]]);fb0=np.pad(fb0,(0,de));fc0=-np.hstack((fc0,fd0))
    inq=(~fast0(CC0,0)) | (~fast0(fc0));inp=((~fast0(AA0,0)) | (~fast0(BB0,0))) & (~inq);dimq=int(inq.sum());dimp=int(inp.sum());zp=ZZ0[:,inp[:-de]];zq=ZZ0[:,inq[:-de]];zc=ZZ1
    AA=np.pad(AA0,((0,1),(0,0)));BBU=np.vstack((BB0,fb0));CCU=np.vstack((CC0,fc0));BBR=np.pad(BB0,((0,1),(0,0)));CCR=np.pad(CC0,((0,1),(0,0)));BBR[-1,list(vv0).index(str(s.const_var))]=-1;fb0[list(vv0).index(str(s.const_var))]=0
    PU=-np.hstack((BBU[:,inq],AA[:,inp]));MU=np.hstack((CCU[:,inq],BBU[:,inp]));PR=-np.hstack((BBR[:,inq],AA[:,inp]));MR=np.hstack((CCR[:,inq],BBR[:,inp]));gg=np.pad([float(s.x_bar)],(dimp+dimq-1,0));R,Q=sl.rq(MU.T);MU=R.T;PU=Q@PU;R,Q=sl.rq(MR.T);MR=R.T;PR=Q@PR;gg=Q@gg
    omg,lam=klein(PU,MU,dimq,force=force);fq0=fc0[inq];fp1=fb0[inp];fq1=fb0[inq];prec=preprocess(PU,MU,PR,MR,gg,fq1,fp1,fq0,omg,lam,s.x_bar,l_max,k_max)
    s.sys=(omg,lam,s.x_bar);s.precalc=prec;s.hx=(zp,zq,zc);s.dimq=dimq;s.dimp=dimp;s.dimeps=de;s.dimy=dimp+dimq;s.inq=inq;s.inp=inp;s.vv=np.hstack((vv0[inp[:-de]],vv0[inq[:-de]]));s.ppar=ppar;s.observable_names=[str(x) for x in m['observables']]
    return s

def find_lk(bmat,bterm,x_bar,q):
    _,lm,km=bterm.shape
    def chk(si,l,k):return bmat[si,l,k]@q+bterm[si,l,k]
    l=k=0;flag=0
    while chk(2,l,0)-x_bar>0:
      l+=1
      if l==lm:break
    if l<lm:
      found=None
      for ll in range(lm):
       for kk in range(1,km):
        if ll and chk(0,ll,kk)-x_bar<0:continue
        if ll>1 and chk(1,ll,kk)-x_bar<0:continue
        if chk(4,ll,kk)-x_bar<0:continue
        if chk(2,ll,kk)-x_bar>0:continue
        if kk>1 and chk(3,ll,kk)-x_bar>0:continue
        found=(ll,kk);break
       if found:break
      if found:l,k=found
      else:
       flag=1;l=k=0;loop=False
       while True:
        k+=1
        if not loop:loop=chk(4,0,k)-x_bar<0
        elif chk(4,0,k)-x_bar>0:break
        if k==km-1:flag=2;break
    else:l=0
    return l,k,flag

def step(s,state,shock,linear=False,get_obs=True):
    omg,lam,xbar=s.sys;pmat,qmat,pterm,qterm,bmat,bterm=s.precalc;dimp,dimq=omg.shape;de=s.dimeps
    qs=np.hstack((state[-dimq+de:],shock))
    l,k,flag=(0,0,0) if linear else find_lk(bmat,bterm,xbar,qs)
    p=pmat[l,k]@qs+pterm[l,k]
    # pydsge.tools.t_func slices the last `dimeps` q rows: shock states are supplied anew each period
    q=qmat[l,k][:-de]@qs+qterm[l,k][:-de]
    po=s.hx[0]@p+s.hx[1]@q+s.hx[2]
    return np.hstack((p,q)),(l,k),flag,p,q

def simulate(s,init,resid,linear=False):
    state=np.array(init,float);X=[state];LK=[];flags=[];full=[]
    for e in resid:
      state,lk,fl,p,q=step(s,state,e,linear=linear);X.append(state);LK.append(lk);flags.append(fl);full.append((p,q))
    return np.array(X),np.array(LK),np.array(flags),full

if __name__=='__main__':
    meta='/mnt/data/dsge6420/rank_spreads_exo_ztrend_SW_BAA_6420_ninit_0_meta.npz';res='/mnt/data/dsge6420/rank_spreads_exo_ztrend_SW_BAA_6420_ninit_0_res.npz'
    with np.load(meta,allow_pickle=False) as z:raw=z['yaml_raw'].item();data=pickle.loads(z['data'].item())
    # preserve historical model shock order from declarations, matching stored residuals
    y0=yaml.safe_load(raw);shock_order=y0['declarations']['shocks']
    yy,m=parse_raw(raw,shock_order);funs,mats=build_lambdas(m,yy)
    print('nvars',len(m['var_ordering']),'neq',len(m['perturb_eq']),'shocks', [str(x) for x in m['shk_ordering']])
    print('vars', [str(x) for x in m['var_ordering']])
    A,B,C,PSI,*_=mats
    print('lag cols',[(str(m['var_ordering'][j]),sum(1 for x in C[:,j] if x!=0)) for j in range(C.shape[1]) if any(x!=0 for x in C[:,j])])
    with np.load(res,allow_pickle=False) as z:pars=z['pars'];init=z['init'];rr=z['resid']
    p=np.median(pars,0);s=gen_system(m,yy,funs,p)
    print('dims dimp dimq dimy, state-input',s.dimp,s.dimq,s.dimy,len(s.observable_names)+s.dimq-s.dimeps,'xbar',s.x_bar)
    print('inp',s.vv[:s.dimp]);print('inq',s.vv[s.dimp:])
    X,LK,flags,full=simulate(s,np.median(init,0),np.median(rr,0))
    print('X',X.shape,'flags',np.unique(flags,return_counts=True),'LK top')
    u,cnt=np.unique(LK,axis=0,return_counts=True);print(sorted(zip(map(tuple,u),cnt),key=lambda x:-x[1])[:20])
    print('data',data.shape, data.index[0],data.index[-1])
    # compare observation components in X against data (state contains obs then q-no-shocks)
    obs=X[:,:len(s.observable_names)]
    print('RMSE obs vs data',np.sqrt(np.nanmean((obs-np.array(data))**2,axis=0)))
    print('corr', [np.corrcoef(obs[:,j],np.array(data)[:,j])[0,1] for j in range(obs.shape[1])])
