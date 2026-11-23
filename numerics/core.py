"""Portable extraction of the audited original-Hamiltonian numerical core.

Provenance: audit/verify_prefactor.py, used by correction/numerics.py and the
stage-3 validation.  Only the old output-only predicted/Cardoso fields were
removed. Boundary recurrence, propagation, extraction, and MP recurrence are
unchanged. No asymptotic transition coefficient enters the numerical problem.

The binary64 method uses mpmath only for endpoint/boundary preparation; its
``dps`` argument is not the precision of the DOP853 propagation. Reported drift,
last terms, and solver tolerances are diagnostics, not certified error bounds.
"""
import numpy as np
from scipy.integrate import solve_ivp
import mpmath as mp


def jost_coefficients(delta, eta, order, dps=80):
    """b=e^(-i*s/(2*delta))*psi_out, b=(f,g), f(0)=0,g(0)=1.
    The power series variable is z=1/tau, at positive infinity.
    -i*delta*z^3*b_z = z*(h+sqrt(1+tau^2)*I)*b.
    """
    with mp.workdps(dps):
        d, a = mp.mpf(str(delta)), mp.mpf(str(eta))
        f, g = [mp.mpc(0)]*(order+2), [mp.mpc(0)]*(order+2)
        g[0], f[1] = mp.mpf(1), -mp.mpf('0.5')
        e = {2*j:mp.binomial(mp.mpf('0.5'),j) for j in range(1,order//2+2)}
        for m in range(1,order+1):
            known = -1j*d*(m-1)*f[m-1]
            known -= sum(e[k]*f[m+1-k] for k in e if k<=m+1)
            known -= 1j*a*d/2*sum((-1)**(j-1)*g[m-2*j]
                                   for j in range(1,m//2+1))
            known /= 2
            other = sum(e[k]*g[m+2-k] for k in e if 4<=k<=m+2)
            other -= 1j*a*d/2*sum((-1)**(j-1)*f[m+1-2*j]
                                   for j in range(1,(m+1)//2+1))
            g[m] = (known+other)/(-1j*d*m)
            f[m+1] = known-g[m]/2
        return f[:order+1],g[:order+1]



def sval(t):
    return t*np.sqrt(1+t*t)+np.arcsinh(t)



def scattering(delta, eps, T=8., order=32, rtol=2e-12, c=1., dps=80,
               real_axis=False, both_sides=False):
    eta = 1-eps
    f,g = jost_coefficients(delta,eta,order,dps)
    height = 0. if real_axis else np.pi/2-c*delta
    if height < 0: raise ValueError('choose c*delta < pi/2')
    X = float(sval(T))
    with mp.workdps(dps):
        target=mp.mpc(X,height)
        fun=lambda t:t*mp.sqrt(1+t*t)+mp.asinh(t)-target
        endpoint=mp.findroot(fun,mp.mpc(T,height/(2*T)))
        y0=mp.findroot(lambda y:y*mp.sqrt(1-y*y)+mp.asin(y)-height,
                       (mp.mpf('0.1'),mp.mpf('0.999999')), solver='anderson') if height else mp.mpf(0)
        center=complex(1j*y0)
        b0=np.array([complex(mp.polyval(list(reversed(v)),1/endpoint))
                     for v in (f,g)])
        tail_last=max(float(abs(v[-1]/endpoint**order)) for v in (f,g))
        endpoint=complex(endpoint)

    def rhs(x,Y):
        t,u,v=Y
        E=np.sqrt(1+t*t)
        a=eta*delta/(2*(1+t*t))
        return np.array([1/(2*E),
            -1j*((t+E)*u+(1+1j*a)*v)/(2*delta*E),
            -1j*((1-1j*a)*u+(-t+E)*v)/(2*delta*E)])

    def propagate(x0, z0, b):
        sol=solve_ivp(rhs,(x0,0),np.r_[z0,b],method='DOP853',
                      rtol=rtol,atol=rtol/100)
        if not sol.success: raise RuntimeError(sol.message)
        Y=sol.y[:,-1]
        drift=abs(Y[0]-center)
        # Correct the small numerical inverse-map drift by one Taylor step.
        deriv=rhs(0,Y)[1:]*(2*np.sqrt(1+Y[0]**2))
        b=Y[1:]+deriv*(center-Y[0])
        return b,drift,sol.nfev

    bout,drift,nfev=propagate(X,endpoint,b0)
    if both_sides:
        bin_,drift2,n2=propagate(-X,-endpoint.conjugate(),b0.conjugate()[::-1])
        W=bin_[0]*bout[1]-bin_[1]*bout[0]
        drift=max(drift,drift2);nfev+=n2
    else:
        W=abs(bout[1])**2-abs(bout[0])**2
    scale=np.exp((np.pi-2*height)/delta)
    C=float(abs(W)**2*scale)
    return dict(delta=delta,eps=eps,eta=eta,T=T,order=order,rtol=rtol,c=c,
                dps=dps,real_axis=real_axis,both_sides=both_sides,
                C=C,Wronskian_real=float(np.real(W)),Wronskian_imag=float(np.imag(W)),
                center_drift=float(drift),last_boundary_term=tail_last,nfev=nfev)



def scattering_mp(delta, eps, T=6, order=80, dps=65, degree=54, c=1):
    """Arbitrary-precision Taylor propagation of the ORIGINAL rational TDSE.

    Multiply TDSE by 1+tau^2 to obtain polynomial coefficients (degree 3).
    Every Taylor coefficient then needs only four preceding coefficients.
    This propagates psi directly, independently of solve_ivp and its ODE.
    """
    with mp.workdps(dps):
        d,e,eta=mp.mpf(str(delta)),mp.mpf(str(eps)),1-mp.mpf(str(eps))
        h=mp.pi/2-c*d
        S=lambda t:t*mp.sqrt(1+t*t)+mp.asinh(t)
        X=S(mp.mpf(T))
        t=mp.findroot(lambda t:S(t)-X-1j*h,mp.mpc(T,h/(2*T)))
        f,g=jost_coefficients(d,eta,order,dps)
        y=[mp.polyval(list(reversed(v)),1/t)*mp.exp(1j*X/(2*d)) for v in (f,g)]
        x=X;steps=0; rejected=0;maxerr=mp.mpf(0)
        tol=mp.power(10,-dps+12)
        while x>0:
            p0,p1=1+t*t,2*t
            ds=[t*p0,1+3*t*t,3*t,mp.mpf(1)]
            os=[p0,2*t,mp.mpf(1),mp.mpf(0)]
            D=[-1j*v/d for v in ds]
            B=[-1j*v/d for v in os]
            cf=[y]
            for n in range(degree):
                u=v=mp.mpc(0)
                for k in range(min(3,n)+1):
                    fu,fv=cf[n-k]
                    u+=D[k]*fu+B[k]*fv
                    v+=B[k]*fu-D[k]*fv
                u+=eta/2*cf[n][1]-p1*n*cf[n][0]
                v-=eta/2*cf[n][0]+p1*n*cf[n][1]
                if n:
                    u-=(n-1)*cf[n-1][0]
                    v-=(n-1)*cf[n-1][1]
                cf.append([u/(p0*(n+1)),v/(p0*(n+1))])
            dx=min(2*d,x)
            while True:
                xn=x-dx
                tn=mp.findroot(lambda z:S(z)-xn-1j*h,t-dx/(2*mp.sqrt(1+t*t)))
                dt=tn-t
                err=max(abs(cf[n][j]*dt**n) for n in range(degree-3,degree+1) for j in (0,1))
                if err<tol:break
                dx/=2;rejected+=1
            y=[mp.polyval([v[j] for v in reversed(cf)],dt) for j in (0,1)]
            t,x=tn,xn;steps+=1;maxerr=max(err,maxerr)
        W=abs(y[1])**2-abs(y[0])**2
        C=W**2*mp.exp(2*c)
        return dict(delta=str(d),eps=str(e),T=T,order=order,dps=dps,degree=degree,c=c,
                    C=mp.nstr(C,dps-8),steps=steps,rejected=rejected,
                    max_last_terms=mp.nstr(maxerr,8),
                    center_drift=mp.nstr(abs(mp.re(t)),8))

