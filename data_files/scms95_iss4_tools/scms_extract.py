from PIL import Image; import numpy as np, sys
im=np.array(Image.open("scms_iss/i7221545.gif").convert("RGB")).astype(int)
# pressure calibration: least squares on gridlines
yl=np.array([109,226,293,341,379,410,435,458,478,496,503.]); pl=np.array([100,200,300,400,500,600,700,800,900,1000,1050.])
b,a=np.polyfit(np.log(pl),yl,1)
P=lambda y: np.exp((y-a)/b)
S=-1.168; X0=372.7; DX=6.13            # px per degC
T=lambda x,y: (x-(X0+S*(y-496)))/DX
print(f"p fit: y = {a:.2f} + {b:.2f} ln p ; max resid {np.abs(a+b*np.log(pl)-yl).max():.2f} px")
def curve(rgb):
    msk=(im==rgb).all(2); out=[]
    for y in range(105,505):
        xs=np.where(msk[y,118:580])[0]+118
        if len(xs): out.append((P(y), np.median([T(x,y) for x in xs]), len(xs)))
    return np.array(out)
tc=curve([0,255,255]); dc=curve([255,255,0])
np.save("scms_iss/tcurve.npy",tc); np.save("scms_iss/dcurve.npy",dc)
def at(c,p): 
    o=np.argsort(c[:,0]); return np.interp(p,c[o,0],c[o,1])
for p in (1000,950,900,850,700,500,400,300,200,150):
    print(f"p={p:5d}  T={at(tc,p):6.2f} C   Td={at(dc,p):6.2f} C")
print("bottom T point:",tc[np.argmax(tc[:,0])], " bottom Td:",dc[np.argmax(dc[:,0])])
