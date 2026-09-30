"""One-file independent reproduction of the ideal pipe transient.

Run `python reproduce.py`. This is not a valve-opening or pipe-rating model.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

R=json.loads((Path(__file__).parent/"reference"/"result.json").read_text())
L=float(R["length_m"])
RHO=float(R["rho_kg_m3"])
A=float(R["wave_speed_m_s"])
V0=float(R["initial_speed_m_s"])
Z=RHO*A
END=float(R["window_s"])
ROUND=2*L/A


def h(t,tc):
    """Drop in prescribed outlet fluid speed, expressed in pressure units."""
    t=np.asarray(t,dtype=float)
    return Z*V0*np.where(t<0,0.,np.minimum(1.,t/tc))


def launched(t,tc):
    t=np.asarray(t,dtype=float)
    value=np.zeros_like(t)
    for n in range(9):
        value += (-1. if n%2 else 1.)*h(t-n*ROUND,tc)
    return value


def valve_exact(t,tc):
    t=np.asarray(t,dtype=float)
    return launched(t,tc)-launched(t-ROUND,tc)


def exact_peak(tc):
    # A finite sum of time-shifted linear ramps has no interior peak except
    # a plateau. Thus every extremum is attained at a shifted breakpoint.
    points={0.,END}
    for n in range(9):
        for q in (n*ROUND,tc+n*ROUND):
            if 0<=q<=END: points.add(round(q,12))
    ts=np.array(sorted(points))
    values=valve_exact(ts,tc)
    return float(np.max(values)),float(np.min(values))


def moc(tc,cells=400,wrong_reservoir=False):
    dx=L/cells
    dt=dx/A
    n=int(round(END/dt))
    p=np.zeros(cells+1)
    u=np.zeros(cells+1)
    pv=np.zeros(n+1)
    reservoir_u=np.zeros(n+1)
    reservoir_p=np.zeros(n+1)
    for k in range(1,n+1):
        rplus=p+Z*u
        rminus=p-Z*u
        pn=np.empty_like(p);un=np.empty_like(u)
        pn[1:-1]=(rplus[:-2]+rminus[2:])/2
        un[1:-1]=(rplus[:-2]-rminus[2:])/(2*Z)
        if wrong_reservoir:
            # Executed negative control: rigid wall instead of fixed pressure.
            un[0]=0.
            pn[0]=rminus[1]
        else:
            pn[0]=0.
            un[0]=-rminus[1]/Z
        un[-1]=-V0*min(1.,k*dt/tc)
        pn[-1]=rplus[-2]-Z*un[-1]
        p,u=pn,un
        pv[k]=p[-1]
        reservoir_u[k]=u[0]
        reservoir_p[k]=p[0]
    return np.arange(n+1)*dt,pv,reservoir_u,reservoir_p


def close(a,b,tol,message):
    if abs(a-b)>tol:
        raise AssertionError(f"{message}: {a} vs {b}")


def run():
    peaks={}
    for name,tc in R["closure_s"].items():
        tc=float(tc)
        exact_max,exact_min=exact_peak(tc)
        ts,pv,ru,rp=moc(tc)
        close(float(pv.max()),exact_max,0.005*exact_max,name+" grid/series peak")
        close(float(pv.min()),exact_min,0.005*Z,name+" grid/series minimum")
        close(exact_max,float(R["positive_valve_peaks_pa"][name]),500.,name+" frozen positive peak")
        close(exact_min,float(R["negative_valve_minima_pa"][name]),500.,name+" frozen minimum")
        if np.max(np.abs(rp))>1e-8:
            raise AssertionError("fixed-pressure reservoir violated")
        arr=ts[np.flatnonzero(np.abs(ru)>1e-7)[0]]
        reflected=ts[np.flatnonzero(np.abs(pv-h(ts,tc))>1e-5)[0]]
        close(arr,float(R["reservoir_arrival_s"]),ts[1]*1.1,name+" reservoir arrival")
        close(reflected,float(R["valve_reflection_return_s"]),ts[1]*1.1,name+" valve return")
        peaks[name]=exact_max
        print(f"{name}: +{exact_max/1e6:.2f} MPa; minimum {exact_min/1e6:+.2f} MPa")
    ratio=peaks["fast"]/peaks["slow"]
    close(ratio,float(R["fast_to_slow_peak_ratio"]),.001,"frozen ratio")
    bad_ts,bad_pv,_,bad_rp=moc(float(R["closure_s"]["fast"]),cells=100,
                               wrong_reservoir=True)
    if np.max(np.abs(bad_rp))<.1*Z or np.max(np.abs(bad_pv-valve_exact(bad_ts,.05)))<.1*Z:
        raise AssertionError("wrong reservoir reflection escaped negative control")
    print("PASS: characteristic grid, analytic waves, frozen peaks, arrivals, wrong-sign control")


if __name__=="__main__":
    run()
