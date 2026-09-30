# A flow stop and its returning pressure wave

This standalone package reproduces one ideal horizontal pipe comparison.
The pipe is 100 m long, starts with a uniform 1 m/s water speed and has a
fixed-pressure upstream reservoir. The downstream boundary is a prescribed
**fluid-speed** ramp to zero over either 0.05 or 0.80 s. It is not a model of
valve aperture, valve hydraulics or a specific installed pressure system.

Run with Python and NumPy:

```sh
python -m pip install -r requirements.txt
python reproduce.py
```

`reproduce.py` reruns the one-dimensional inviscid characteristic equations,
independently constructs the exact traveling-wave series, searches the
piecewise-linear peak over 0–1.20 s, checks the frozen `reference/result.json`,
and exits nonzero on a failed check. The source is intentionally a single
readable script. The photograph and video editing are not inputs to this
calculation.

The positive valve-end pressure rise is 1.00 MPa for the 0.05 s stop and
0.25 MPa for the 0.80 s stop. These peaks need not occur simultaneously.
Pressure changes are relative to an assumed 2.00 MPa absolute baseline.
In the fast case the outlet fluid speed is zero after 0.05 s but the returning
wave reaches it at 0.20 s and its pressure perturbation later becomes negative.
The slow case's pressure perturbation returns to zero once its prescribed
flow stop finishes; persistent post-stop variation is not claimed for it.

The stipulated effective wave speed is 1000 m/s. No friction, real valve
pressure–flow law, nonlinear pipe compliance, junction, vapor cavity or
measured pressure history is represented. The reported difference is **not**
a real pipe safety rating or operating recommendation.

Background on characteristic water-hammer modelling: [NASA's technical
report](https://ntrs.nasa.gov/citations/20040076876) and [USC's research
dissertation](https://cee.usc.edu/wp-content/uploads/2017/02/Dissertation_Kim2012.pdf).
