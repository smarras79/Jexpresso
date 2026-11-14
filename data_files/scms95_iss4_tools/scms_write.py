import numpy as np
d=np.load("scms_iss/profile.npz"); p,T,Td,qv=d['p'].copy(),d['T'].copy(),d['Td'].copy(),d['qv'].copy()
Rd,cp,g=287.04,1004.,9.80616; k=Rd/cp
TK=T+273.15; th=TK*(1000/p)**k
# 1) bottom: the surface temperature sensor is sun-heated below ~107 m (p > 1003.9 hPa), so theta there is
#    replaced by the mixed-layer value.
i107=np.argmin(abs(p-1003.9))
th[:i107]=th[i107]
# 2) light 3-point smoothing of theta and the mixing ratio (pixel quantisation), keeping the end points
def sm(a): b=a.copy(); b[1:-1]=0.25*a[:-2]+0.5*a[1:-1]+0.25*a[2:]; return b
th=sm(th); qv=sm(qv)
# 3) mixed-layer moisture: the measured surface layer (20 g/kg at the ground) cannot be kept with a mixed theta
#    (theta_v would fall with height: statically unstable), and dropping it (v1) left the layer 0.4-0.5 g/kg too
#    dry. Use a uniform mixing ratio equal to the archive analysis' lowest-50-mb mean (17.1 g/kg) over that
#    layer, blended linearly in p back to the measured profile over the next 30 hPa.
R_ML=17.1e-3; p_ml=p[0]-50.0; p_bl=p_ml-30.0
a=np.clip((p-p_bl)/(p_ml-p_bl),0.0,1.0)
qv=a*R_ML+(1-a)*qv
TK=th*(p/1000)**k
Tv=TK*(1+0.608*qv/(1+qv)); z=np.zeros_like(p)
for i in range(1,len(p)): z[i]=z[i-1]+Rd*0.5*(Tv[i]+Tv[i-1])/g*np.log(p[i-1]/p[i])
thv=Tv*(1000/p)**k; low=z<1500
dthv=np.diff(thv[low])
print(f"theta_v below 1.5 km: min step {dthv.min():+.3f} K ({'stable' if dthv.min()>=-1e-3 else 'UNSTABLE step at z=%.0f m' % z[np.argmin(dthv)]})")
m=p>=p[0]-50
print(f"lowest-50-mb mean: theta={th[m].mean():.2f} K (archive 301.4), r={qv[m].mean()*1e3:.2f} g/kg (archive 17.1)")
# 4) the model reads the moisture column as specific humidity qv = qt: convert from mixing ratio
q=qv/(1+qv)
out="scms_iss/LasherTrapp2001_SCMS_ISS4_22Jul1995_1545_v3.dat"
with open(out,"w") as f:
    for zi,pi,ti,qi in zip(z,p,TK-273.15,q*1e3):
        f.write(f"{zi:.1f} {pi*100:.1f} {ti:.3f} {qi:.4f} 0.0 0.0\n")
print(out, len(p),"levels, z top",round(z[-1]), "m  (column 4: specific humidity, g/kg)")
