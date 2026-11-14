import numpy as np
tc=np.load("scms_iss/tcurve.npy"); dc=np.load("scms_iss/dcurve.npy")
Rd,Rv,cp,g,Lv=287.04,461.5,1004.,9.80616,2.5e6; eps=Rd/Rv
def es(Tc): return 611.2*np.exp(17.67*Tc/(Tc+243.5))      # Bolton, Pa
# common pressure grid: from surface (max p with both curves) to 100 hPa, using the pixel rows of the T curve
ps=min(tc[:,0].max(),dc[:,0].max())
o1=np.argsort(tc[:,0]); o2=np.argsort(dc[:,0])
p=np.unique(np.round(np.concatenate([tc[:,0],dc[:,0]]),2)); p=np.sort(p[(p<=ps)&(p>=100)])[::-1]
T=np.interp(p,tc[o1,0],tc[o1,1]); Td=np.interp(p,dc[o2,0],dc[o2,1])
Td=np.minimum(Td,T)
e=es(Td); qv=eps*e/(p*100-e)                              # kg/kg (mixing ratio)
TK=T+273.15; Tv=TK*(1+0.608*qv/(1+qv))
z=np.zeros_like(p)
for i in range(1,len(p)): z[i]=z[i-1]+Rd*0.5*(Tv[i]+Tv[i-1])/g*np.log(p[i-1]/p[i])
th=TK*(1000/p)**(Rd/cp)
# 50 hPa mean layer (pressure-weighted ~ uniform in p)
m=p>=ps-50
print(f"surface p={ps:.1f} hPa T={T[0]:.2f}C Td={Td[0]:.2f}C;  50-mb mean theta={th[m].mean():.2f} K (archive 301.4), w={qv[m].mean()*1e3:.2f} g/kg (archive 17.1)")
np.savez("scms_iss/profile.npz",p=p,T=T,Td=Td,qv=qv,z=z)
for pp in (1000,900,850,700,600,500,400,300,250,200,150,100):
    i=np.argmin(abs(p-pp)); print(f"  p={p[i]:6.1f} z={z[i]:7.0f} m T={T[i]:6.1f} Td={Td[i]:6.1f} qv={qv[i]*1e3:6.2f} g/kg RH={e[i]/es(T[i])*100:5.1f}%")
