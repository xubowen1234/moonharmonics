"""Independent SciPy/SVD oracle; no generated fixtures are production implementation."""
import argparse, json, math, os, pathlib, subprocess
import numpy as np
from scipy.special import sph_harm_y, roots_legendre, eval_legendre
from scipy.linalg import lstsq, expm
from scipy.spatial.transform import Rotation

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
    assert len(records)==56, f'unexpected reference output count {len(records)}'
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
        elif kind in ('dh_analysis','dh_synthesis'):
            degree=16;n=degree+1
            theta=np.repeat(np.arange(2*n)*math.pi/(2*n),2*n)
            phi=np.tile(np.arange(2*n)*math.pi/n,2*n)
            dh_basis=basis(degree,theta,phi)
            dh_weights=(2/n)*np.sin(np.arange(2*n)*math.pi/(2*n))*np.asarray([
                sum(math.sin((2*k+1)*j*math.pi/(2*n))/(2*k+1) for k in range(n))
                for j in range(2*n)
            ])
            sample_weights=np.repeat(dh_weights*math.pi/n,2*n)
            if kind=='dh_synthesis':
                check(kind,v,dh_basis@coefficients,3e-12)
            else:
                x=np.sin(theta)*np.cos(phi);y=np.sin(theta)*np.sin(phi);z=np.cos(theta)
                samples=np.exp(.5*z)+x*y+.2*x
                check(kind,v,dh_basis.T@(sample_weights*samples),3e-12)
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
        elif kind in ('so3_encoding','so3_field'):
            a,b,g=v[0][:3]
            rotation=(Rotation.from_rotvec([0,0,a])*Rotation.from_rotvec([0,b,0])*Rotation.from_rotvec([0,0,g]))
            def angles(cartesian):
                return np.arctan2(np.hypot(cartesian[...,0],cartesian[...,1]),cartesian[...,2]), np.arctan2(cartesian[...,1],cartesian[...,0])
            if kind=='so3_encoding':
                t,p=v[0][3:]
                d=rotation.apply([math.sin(t)*math.cos(p),math.sin(t)*math.sin(p),math.cos(t)])
                rt,rp=angles(d)
                check(kind,v[1],basis(64,rt,rp),5e-10)
                t2,p2=v[2]
                check('so3_direction',[math.sin(t2)*math.cos(p2),math.sin(t2)*math.sin(p2),math.cos(t2)],d,3e-14)
            else:
                samples=np.asarray(v[2]);t=samples[:,0];p=samples[:,1]
                points=np.column_stack([np.sin(t)*np.cos(p),np.sin(t)*np.sin(p),np.cos(t)])
                it,ip=angles(rotation.inv().apply(points))
                canonical=np.sin(np.arange(1,4226))/np.arange(1,4226)
                check('so3_field_values',samples[:,2],basis(64,it,ip)@canonical,5e-11)
                check('so3_field_output_basis',samples[:,2],basis(64,t,p)@np.asarray(v[1]),5e-11)
                # Independent dense Pade expm (not the library's stepped Taylor action).
                l=64;m=np.arange(-l,l+1);band=canonical[l*l:(l+1)**2]
                ac=np.zeros(2*l+1,dtype=complex);ac[l]=band[l]
                for k in range(1,l+1):
                    ac[l+k]=(band[l+k]-1j*band[l-k])/math.sqrt(2)
                    ac[l-k]=(-1.)**k*np.conj(ac[l+k])
                edge=.5*np.sqrt((l-m[:-1])*(l+m[:-1]+1))
                generator=np.diag(edge,1)-np.diag(edge,-1)
                ac=np.exp(-1j*m*a)*(expm(b*generator)@(np.exp(-1j*m*g)*ac))
                expected=np.empty(2*l+1);expected[l]=ac[l].real
                for k in range(1,l+1):
                    expected[l+k]=math.sqrt(2)*ac[l+k].real
                    expected[l-k]=-math.sqrt(2)*ac[l+k].imag
                check('so3_band64_pade',v[1][l*l:(l+1)**2],expected,2e-13)
                actual_energy=[np.dot(np.asarray(v[1])[j*j:(j+1)**2],np.asarray(v[1])[j*j:(j+1)**2]) for j in range(65)]
                expected_energy=[np.dot(canonical[j*j:(j+1)**2],canonical[j*j:(j+1)**2]) for j in range(65)]
                check('so3_band_energy',actual_energy,expected_energy,2e-13)
        elif kind=='convention':
            ni,pi,oi=map(int,v[0]);canonical=np.sin(np.arange(1,82))/np.arange(1,82)
            external=np.cos(np.arange(1,82))/np.arange(1,82)
            exported=np.empty(81);imported=np.empty(81)
            for l in range(9):
                order=list(range(-l,l+1)) if oi==0 else ([0]+[m for k in range(1,l+1) for m in (k,-k)] if oi==1 else list(range(l+1))+list(range(-1,-l-1,-1)))
                factor=1 if ni==0 else math.sqrt(4*math.pi/(1 if ni==1 else 2*l+1))
                for offset,m in enumerate(order):
                    i=l*l+l+m;j=l*l+offset;sign=(-1.)**abs(m) if pi else 1
                    exported[j]=canonical[i]/(factor*sign)
                    imported[i]=external[j]*factor*sign
            check('convention_export',v[1],exported,3e-15)
            check('convention_import',v[2],imported,3e-15)
        else: raise AssertionError(kind)
    report={'oracle':'SciPy 1.16.3 spherical harmonics, spatial rotations, Pade matrix exponential, quadrature and SVD least squares','target':args.target,'records':len(records),'scalar_comparisons':int(comparisons),'maximum_absolute_errors':metrics,'passed':True}
    text=json.dumps(report,indent=2)
    if args.output:args.output.write_text(text+'\n',encoding='utf-8')
    print(text)
if __name__=='__main__': main()
