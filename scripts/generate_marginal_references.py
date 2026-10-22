"""Independent high-precision reference data; mpmath is not a runtime/test dependency."""
from pathlib import Path
import hashlib,json,sys
import mpmath as mp

assert mp.__version__=='1.3.0'
maximum=sys.float_info.max;tiny=float.fromhex('0x0.0000000000001p-1022')
cases=[]
def add(name,kind,lower,upper,nominal,probability,sigma=None,precisions=(100,180)):
    cases.append((name,kind,lower,upper,nominal,probability,sigma,precisions))
for p in (.01,.25,.5,.75,.99):
    add(f'log_wide_{p}','log-uniform',tiny,maximum,1.,p,precisions=(400,500))
    add(f'normal_bounded_{p}','truncated-normal',-1.,1.,0.,p,1.)
    add(f'normal_shifted_{p}','truncated-normal',1.,4.,2.,p,.5)
for p in (tiny,.25,.75,float.fromhex('0x1.fffffffffffffp-1')):
    add(f'normal_full_{p.hex()}','truncated-normal',-maximum,maximum,0.,p,1.,(400,500))
for p in (.25,.5,.75):
    add(f'normal_narrow_{p}','truncated-normal',0.,tiny,0.,p,maximum,(2200,2400))
rows=[]
for name,kind,lower,upper,nominal,p,sigma,precisions in cases:
    answers=[];zeros=[]
    for digits in precisions:
        with mp.workdps(digits):
            a,b,m,u=map(mp.mpf,(lower,upper,nominal,p))
            if kind=='log-uniform':q=mp.exp((1-u)*mp.log(a)+u*mp.log(b))
            else:
                s=mp.mpf(sigma);alpha=(a-m)/(s*mp.sqrt(2));beta=(b-m)/(s*mp.sqrt(2))
                q=m+s*mp.sqrt(2)*mp.erfinv(mp.erf(alpha)+u*(mp.erf(beta)-mp.erf(alpha)))
            answers.append(float(q).hex());zeros.append(q==0)
    assert answers[0]==answers[1],(name,answers)
    exact_zero=(kind=='truncated-normal' and nominal==0 and lower==-upper and p==.5)
    row={'case':name,'distribution':kind,'lower_hex':lower.hex(),'upper_hex':upper.hex(),
         'nominal_hex':nominal.hex(),'probability_hex':p.hex(),
         'sigma_hex':None if sigma is None else sigma.hex(),
         'expected_float_hex':answers[1],'exact_zero_by_symmetry':exact_zero,
         'positive_value_rounds_to_zero':float.fromhex(answers[1])==0 and not exact_zero,
         'reference_decimal_precisions':list(precisions)}
    rows.append(row)
dest=Path('tests/fixtures/marginal_quantile_reference.json');dest.parent.mkdir(parents=True,exist_ok=True)
report={'scope':'Independent mpmath centered-erf inverse and log/exp reference values, exact binary float inputs; two precisions must agree on the exported binary value. This is finite reference evidence, not a universal rounding proof or physical validation.',
        'generator':'mpmath 1.3.0','cases':rows}
dest.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'reference_cases':len(rows),'sha256':hashlib.sha256(dest.read_bytes()).hexdigest()}))
