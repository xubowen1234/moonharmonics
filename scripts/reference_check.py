"""Independent SciPy/SVD oracle; no generated fixtures are production implementation."""
import argparse, json, math, os, pathlib, subprocess
import numpy as np
from scipy.special import sph_harm_y, roots_legendre, eval_legendre
from scipy.linalg import lstsq

def basis(degree, theta, phi):
    result=[]
    for l in range(degree+1):
        for m in range(-l,l+1):
            y=sph_harm_y(l,abs(m),theta,phi)
            result.append(y.real if m==0 else math.sqrt(2)*(y.real if m>0 else y.imag))
    return np.asarray(result).T

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--target',default='wasm-gc')
    parser.add_argument('--output',type=pathlib.Path)
    args=parser.parse_args()
    root=pathlib.Path(__file__).resolve().parents[1]
    run=subprocess.run(['moon','run','tools/reference','--target',args.target],cwd=root,text=True,capture_output=True,encoding='utf-8')
    if run.returncode: raise RuntimeError(run.stderr)
    records=[json.loads(line) for line in run.stdout.splitlines() if line.startswith('{')]
    assert len(records)==29, f'unexpected reference output count {len(records)}'
    metrics={};comparisons=0
    def check(name,actual,expected,tolerance):
        nonlocal comparisons
        a=np.asarray(actual);e=np.asarray(expected)
        assert a.shape==e.shape,(name,a.shape,e.shape)
        assert np.all(np.isfinite(a)),name
        error=float(np.max(np.abs(a-e)))
        metrics[name]=max(metrics.get(name,0.0),error)
        comparisons+=a.size
        assert error<=tolerance, (name,error,tolerance)
    nodes,weights=roots_legendre(17)
    theta=np.repeat(np.arccos(nodes),33)
    phi=np.tile(np.arange(33)*2*math.pi/33,17)
    matrix=basis(16,theta,phi)
    coefficients=np.sin(np.arange(1,290))/np.arange(1,290)
    quadrature_weights=np.repeat(weights*2*math.pi/33,33)
    for record in records:
        kind=record['kind'];v=record['value']
        if kind=='basis':
            t,p,values=v
            check(kind,values,basis(64,t,p),2e-11)
        elif kind=='quadrature':
            x,w=roots_legendre(256)
            check('quadrature_nodes',v[0],x,3e-15)
            check('quadrature_weights',v[1],w,5e-14)
        elif kind=='synthesis':
            check(kind,v,matrix@coefficients,2e-12)
        elif kind=='analysis':
            x=np.sin(theta)*np.cos(phi);y=np.sin(theta)*np.sin(phi);z=np.cos(theta)
            samples=np.exp(.5*z)+x*y+.2*x
            check(kind,v,matrix.T@(quadrature_weights*samples),2e-12)
        elif kind=='gradient':
            t,p,gx,gy,gz=v
            dt=0.;dp=0.
            for l in range(17):
                for m in range(-l,l+1):
                    _,deriv=sph_harm_y(l,abs(m),t,p,diff_n=1)
                    real_deriv=deriv.real if m==0 else math.sqrt(2)*(deriv.real if m>0 else deriv.imag)
                    dt+=coefficients[l*l+l+m]*real_deriv[0]
                    dp+=coefficients[l*l+l+m]*real_deriv[1]/math.sin(t)
            expected=[math.cos(t)*math.cos(p)*dt-math.sin(p)*dp,math.cos(t)*math.sin(p)*dt+math.cos(p)*dp,-math.sin(t)*dt]
            check(kind,[gx,gy,gz],expected,2e-10)
        elif kind=='lambertian':
            x,w=roots_legendre(40);t=(x+1)/2
            expected=[math.pi*np.dot(w,t*eval_legendre(l,t)) for l in range(65)]
            check(kind,v,expected,3e-13)
        elif kind=='fit':
            x,w=roots_legendre(6)
            t=np.repeat(np.arccos(x),11);p=np.tile(np.arange(11)*2*math.pi/11,6)
            b=basis(2,t,p);i=np.arange(len(t));weight=(i%7+1).astype(float)
            values=2+np.sin(t)*np.cos(p)-3*np.sin(t)*np.sin(p)+np.cos(t)+.01*np.sin(3*i)
            expected,_,rank,_=lstsq(b*np.sqrt(weight)[:,None],values*np.sqrt(weight),lapack_driver='gelsd')
            assert rank==9
            rms=math.sqrt(np.dot(weight,(b@expected-values)**2)/sum(weight))
            check('fit_coefficients',v[0],expected,2e-12)
            check('fit_weighted_rms',v[1],rms,2e-13)
        else: raise AssertionError(kind)
    report={'oracle':'SciPy 1.16.3 spherical harmonics, quadrature and SVD least squares','target':args.target,'records':len(records),'scalar_comparisons':int(comparisons),'maximum_absolute_errors':metrics,'passed':True}
    text=json.dumps(report,indent=2)
    if args.output:args.output.write_text(text+'\n',encoding='utf-8')
    print(text)
if __name__=='__main__': main()
