import ephem
import ebf
import numpy as np
import matplotlib.pyplot  as plt
import os
import asctab

__version__ = "0.0.1"


def bitison(x,bitno):
    return x&(2**bitno)!=0

def bitisoff(x,bitno):
    return x&(2**bitno)==0

def bit_on(x,bitno):
    return x|(2**bitno)

def bit_off(x,bitno):
    return x&(~(2**bitno))

def bit_toggle(x,bitno):
    return x^(2**bitno)

def angular_area(angle):
    return (1-np.cos(np.radians(angle)))*64800/(np.pi)

# def mjd2djd(x):
#     return x-15019.5

# def djd2mjd(x):
#     return x+15019.5

def radec_decimal(ra,dec):
    a=ephem.Equatorial(0.0,0.0)
    if (hasattr(ra,"__len__") == False):
        a.from_radec(ra,dec)
        temp=np.degrees(a.to_radec())
        return [temp[0],temp[1]]
    else:
        rao=np.zeros(len(ra),dtype='float64')
        deco=np.zeros(len(ra),dtype='float64')
        for i in range(len(ra)):
            a.from_radec(ra[i],dec[i])
            temp=np.degrees(a.to_radec())
            rao[i]=temp[0]
            deco[i]=temp[1]
        return rao,deco


def angle_znorm(angle):
    if (hasattr(angle,"__len__") == False):
        if angle>180.0:
            angle=angle-360.0
    else:
        angle=np.array(angle)
        ind=np.where(angle>180.0)[0]
        angle[ind]=angle[ind]-360.0
    return angle

def angdist_xyz(x1,y1,z1,x2,y2,z2):
    return 2*np.arcsin(np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)/2.0)*180.0/np.pi

def angsep(l1,b1,l2,b2):
    x1,y1,z1=lbr2xyz(l1,b1)
    x2,y2,z2=lbr2xyz(l2,b2)
    return np.degrees(2*np.arcsin(np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)/2.0))

def sspring(date=None):
    observer=ephem.Observer()
    observer.long, observer.lat = '149.9356', '-31.2733'
    observer.date = '2014/01/01 12:00:00' 
    if date != None:
        observer.date = date
    else:
        # 12 noon at siding spring
        observer.date =  observer.date-np.degrees(observer.long)/360.0 

    return observer

def transit_ra(observer,date,duration):
    observer.date=date
    observer.date+=(0.5*duration/24.0)
    return np.degrees(observer.radec_of(np.pi,np.pi/2)[0].norm)


def lbr2xyz(l,b,r=1.0):
#    dtor=np.pi/180.0
    l=np.radians(l)
    b=np.radians(b)
    return [r*np.cos(b)*np.cos(l),r*np.cos(b)*np.sin(l),r*np.sin(b)]

def xyz2lbr(x,y,z):
#    dtor=np.pi/180.0
    rc=np.sqrt(x*x+y*y)
    l=np.arctan2(y,x)
    b=np.arctan(z/rc)
    r=np.sqrt(x*x+y*y+z*z)
    l=np.degrees(l)
    b=np.degrees(b)
    return [l,b,r]


def vxyz2lbr(px,py,pz,vx,vy,vz):
    r=np.sqrt(px*px+py*py+pz*pz)
    ind=np.where(r == 0.0)[0]
    r[ind]=1.0
    px=px/r
    py=py/r
    pz=pz/r
    rc=np.sqrt(px*px+py*py)
    ind=np.where(rc == 0.0)[0]
    rc[ind]=1.0
    tm_00=-py/rc; tm_01=-pz*px/rc; tm_02= rc*px/rc
    tm_10= px/rc; tm_11=-pz*py/rc; tm_12= rc*py/rc
    tm_20= 0.0  ; tm_21= rc  ; tm_22= pz
    vl=(vx*tm_00+vy*tm_10+vz*tm_20)
    vb=(vx*tm_01+vy*tm_11+vz*tm_21)
    vr=(vx*tm_02+vy*tm_12+vz*tm_22)
    return vl,vb,vr


def vlbr2xyz(l,b,r,vl,vb,vr):
    l=np.radians(l)
    b=np.radians(b)
    tm_00=-np.sin(l) ; tm_01=-np.sin(b)*np.cos(l) ; tm_02= np.cos(b)*np.cos(l)
    tm_10= np.cos(l) ; tm_11=-np.sin(b)*np.sin(l) ; tm_12= np.cos(b)*np.sin(l)
    tm_20= 0.0       ; tm_21= np.cos(b)           ; tm_22= np.sin(b)
    vx=vl*tm_00+vb*tm_01+vr*tm_02
    vy=vl*tm_10+vb*tm_11+vr*tm_12
    vz=vl*tm_20+vb*tm_21+vr*tm_22
    # px=r*np.cos(b)*np.cos(l)
    # py=r*np.cos(b)*np.sin(l)
    # pz=r*np.sin(b)
    return vx,vy,vz

def radius(x,y,z):
    return np.sqrt(x*x+y*y+z*z)

def rotate_xy(th,x,y):
    x1=x*1.0 
    y1=y*1.0
    th=np.radians(th)
    x=x1*np.cos(th)+y1*np.sin(th)
    y=-x1*np.sin(th)+y1*np.cos(th)
    return [x,y]

def plotcirc1(l,b,angle,norm=False,fill=False,fc='y',ec=None,nsize=64):
    th=np.linspace(0,2*np.pi,nsize)
    angle=np.radians(angle)
    l1=np.array(l)
    b1=np.array(b)
    if l1.ndim==0:
        l1=np.array([l1])
    if b1.ndim==0:
        b1=np.array([b1])
    angle1=np.zeros(l1.size,dtype='float64')
    angle1[:]=angle
    fc1=np.zeros(l1.size,dtype='S2')
    fc1[:]=fc
    for i in range(l1.size):
        xc=np.cos(th)*np.sin(angle1[i])
        yc=np.sin(th)*np.sin(angle1[i])
        zc=np.cos(angle1[i])
        [zc,xc]=rotate_xy(-(90.0-b1[i]),zc,xc)
        [xc,yc]=rotate_xy(-l1[i],xc,yc)
        [lc,bc,rc]=xyz2lbr(xc,yc,zc)
        if norm == True:
            lc=lc%360.0
        d=np.max(lc)-np.min(lc)
        if d>358.0:
            ind=np.where(lc<180.0)
            lc[ind]=lc[ind]+360.0
            plt.gca().add_patch(plt.Polygon(np.transpose(np.array([lc,bc])),fill=fill,fc=fc1[i],ec=ec))
            lc=lc-360.0

        plt.gca().add_patch(plt.Polygon(np.transpose(np.array([lc,bc])),fill=fill,fc=fc1[i],ec=ec))
#        plt.draw()



def plotcirc(l,b,angle,norm=False,ls=None,lw=None):
    l1=np.array(l)
    b1=np.array(b)
    if l1.ndim==0:
        l1=np.array([l1])
    if b1.ndim==0:
        b1=np.array([b1])
    angle1=np.zeros(l1.size,dtype='float64')
    angle1[:]=np.radians(angle)

    th=np.linspace(0,2*np.pi,64)
    if ls == None:
        ls='k-'
    if lw == None:
        lw=1.0

    for i in range(l1.size):
        xc=np.cos(th)*np.sin(angle1[i])
        yc=np.sin(th)*np.sin(angle1[i])
        zc=np.cos(angle1[i])
        [zc,xc]=rotate_xy(-(90.0-b1[i]),zc,xc)
        [xc,yc]=rotate_xy(-l1[i],xc,yc)
        [lc,bc,rc]=xyz2lbr(xc,yc,zc)
        if norm == True:
            lc=lc%360.0
        d=np.max(lc)-np.min(lc)
        if d>358.0:
            ind=np.where(lc<180.0)
            lc[ind]=lc[ind]+360.0
            plt.plot(lc,bc,ls,lw=lw)
            lc=lc-360.0
        plt.plot(lc,bc,ls,lw=lw)
    

def make_patch(l,b,radius,nsize):
    phi=(np.random.ranf(nsize)-0.5)*180.0*2
    temp=1-np.cos(np.radians(radius))
    theta=np.degrees(np.arccos(1-np.random.ranf(nsize)*temp))
    [x,y,z]=lbr2xyz(phi,90.0-theta)
    [z,x]=rotate_xy(-(90.0-b),z,x)
    [x,y]=rotate_xy(-l,x,y)
    return xyz2lbr(x,y,z)

def move_patch(l1,b1,l2,b2,l3,b3):
    [x,y,z]=lbr2xyz(l1,b1)
    [z,x]=rotate_xy(-(b2-b3),z,x)
    [x,y]=rotate_xy(-(l3-l2),x,y)
    [l1,b1,r1]=xyz2lbr(x,y,z)
    return [l1,b1]


def plot_spatch(l,b,lc,bc,radius):
    [x,y,z]=lbr2xyz(l,b,1.0)
    [x,y]=rotate_xy(lc,x,y)
    [z,x]=rotate_xy((90.0-bc),z,x)
    plt.clf()
    plt.plot(np.degrees(x),np.degrees(y),'k.')
    plt.axis([-radius,radius,-radius,radius])

def kernel_density2(r,r0):
#    fd=[0.75000113,0.63661975,0.59683102,0.60792706,0.66492018,0.77403678]
    r0=r0*r0
    return np.sum(1-r*r/r0)*0.63/r0


def delta_airmass(observer,date,duration,ra,dec):
    """
    Compute change in air mass for a given date, duration and ra, dec
    """
    return get_airmass(observer,date,duration,ra,dec)

def get_airmass(observer,date,duration,ra,dec):
    """
    Compute change in air mass for a given date, duration and ra, dec
    """
    transitra=transit_ra(observer,date,duration)
    midpt=np.abs(ra-transitra)
    track_angle=360.0*duration/24.0
    lat=np.degrees(observer.lat)
    radiff2=midpt+track_angle*0.5
    radiff1=np.clip(midpt-track_angle*0.5,0,360.0)
    return [np.abs(airmass(radiff1,dec,lat)-airmass(radiff2,dec,lat)),transitra]

def airmass(radiff,dec,lat=-31.2733):
    th0=np.radians(90.0-lat)
    dec=np.radians(dec)
    alpha=np.radians(180.0-np.abs(radiff))
    cosz=(-np.cos(dec)*np.cos(alpha)*np.sin(th0)+np.sin(dec)*np.cos(th0))
    return (1.002432*cosz*cosz+0.148386*cosz+0.0096467)/(cosz**3+0.149864*cosz*cosz+0.0102963*cosz+0.000303978)

# def delta_airmass(radiff1,radiff2,dec,lat=-31.2733):
#     return np.abs(airmass(radiff1,dec,lat)-airmass(radiff2,dec,lat))

# def delta_airmass_track(midpt,track_angle,dec,lat=-31.2733):
#     radiff2=midpt+track_angle*0.5
#     radiff1=np.clip(midpt-track_angle*0.5,0,360.0)
#     return np.abs(airmass(radiff1,dec,lat)-airmass(radiff2,dec,lat))


def ephem_object(ra,dec):
    temp = ephem.readdb("GalahF,f|S,00:00:00,00:00:00,12.00,2000")
    temp._ra+=np.radians(ra)
    temp._dec+=np.radians(dec)
    return temp

def gal2ecl(l, b,epoch='2000'):
    ra,dec=gal2equ(l,b)
    return equ2ecl(ra,dec)

def ecl2gal(l, b,epoch='2000'):
    ra,dec=ecl2equ(l,b)
    return equ2gal(ra,dec)


def gal2equ(l, b,epoch='2000'):   
   # al_gp = np.radians(192.85948) 
   # de_gp = np.radians(27.12825)
   # l_cp = np.radians(122.932e0) 
   g=ephem.Galactic(0.0,np.pi/2.0,epoch=str(epoch))
   [al_gp,de_gp]=g.to_radec()
   g.from_radec(0.0,np.pi/2.0)
   l_cp=g.lon.znorm
   
   l = np.radians(l) 
   b = np.radians(np.clip(b,-89.999,89.999))
   
   al = al_gp + np.arctan2(np.cos(b) * np.sin(l_cp - l), np.cos(de_gp) * np.sin(b) - np.sin(de_gp) * np.cos(b) * np.cos(l_cp - l))   
   temp = np.sin(al - al_gp) * (np.sin(de_gp) * np.sin(b) + np.cos(de_gp) * np.cos(b) * np.cos(l_cp - l)) / (np.cos(b) * np.sin(l_cp - l))   
   de = np.degrees(np.arctan(temp))   
   al = np.degrees(al)%360.0    
   al=al%360.0

   return [al,de]


def equ2gal(al, de,epoch='2000'):   
   g=ephem.Galactic(0.0,np.pi/2.0,epoch=str(epoch))
   [al_gp,de_gp]=g.to_radec()
   g.from_radec(0.0,np.pi/2.0)
   l_cp=g.lon.znorm
   al = np.radians(al) 
   de = np.radians(de)

   l=l_cp-np.arctan2(np.cos(de)*np.sin(al-al_gp),np.cos(de_gp)*np.sin(de)-np.sin(de_gp)*np.cos(de)*np.cos(al-al_gp))
   temp=np.sin(l_cp-l)*(np.sin(de_gp)*np.sin(de)+np.cos(de_gp)*np.cos(de)*np.cos(al-al_gp))/(np.cos(de)*np.sin(al-al_gp))
   b=np.degrees(np.arctan(temp))
   l=np.degrees(l)%360.0
   
   return [l,b]


def equ2gal_pm(al, de, mua,mud,epoch='2000'):   
   g=ephem.Galactic(0.0,np.pi/2.0,epoch=str(epoch))
   [al_gp,de_gp]=g.to_radec()
   g.from_radec(0.0,np.pi/2.0)
   l_cp=g.lon.znorm
   al = np.radians(al) 
   de = np.radians(de)

   l=l_cp-np.arctan2(np.cos(de)*np.sin(al-al_gp),np.cos(de_gp)*np.sin(de)-np.sin(de_gp)*np.cos(de)*np.cos(al-al_gp))
   temp=np.sin(l_cp-l)*(np.sin(de_gp)*np.sin(de)+np.cos(de_gp)*np.cos(de)*np.cos(al-al_gp))/(np.cos(de)*np.sin(al-al_gp))
   b=np.degrees(np.arctan(temp))
   l=np.degrees(l)%360.0
   c1=np.sin(de_gp)*np.cos(de)-np.cos(de_gp)*np.sin(de)*np.cos(al-al_gp)
   c2=np.cos(de_gp)*np.sin(al-al_gp)
   cosb=np.sqrt(c1*c1+c2*c2)
   mul=(c1*mua+c2*mud)/cosb
   mub=(-c2*mua+c1*mud)/cosb   
   return [l,b,mul,mub]


def ecl2equ(l, b,epoch='2000'):   
   g=ephem.Ecliptic(0.0,np.pi/2.0,epoch=str(epoch))
   [al_gp,de_gp]=g.to_radec()
   g.from_radec(0.0,np.pi/2.0)
   l_cp=g.lon.znorm
   
   l = np.radians(l) 
   b = np.radians(b)
   
   al = al_gp + np.arctan2(np.cos(b) * np.sin(l_cp - l), np.cos(de_gp) * np.sin(b) - np.sin(de_gp) * np.cos(b) * np.cos(l_cp - l))   
   temp = np.sin(al - al_gp) * (np.sin(de_gp) * np.sin(b) + np.cos(de_gp) * np.cos(b) * np.cos(l_cp - l)) / (np.cos(b) * np.sin(l_cp - l))   
   de = np.degrees(np.arctan(temp))   
   al = np.degrees(al)%360.0    
   al=al%360.0

   return [al,de]


def equ2ecl(al, de,epoch='2000'):   
   g=ephem.Ecliptic(0.0,np.pi/2.0,epoch=str(epoch))
   [al_gp,de_gp]=g.to_radec()
   g.from_radec(0.0,np.pi/2.0)
   l_cp=g.lon.znorm
   al = np.radians(al) 
   de = np.radians(de)

   l=l_cp-np.arctan2(np.cos(de)*np.sin(al-al_gp),np.cos(de_gp)*np.sin(de)-np.sin(de_gp)*np.cos(de)*np.cos(al-al_gp))
   temp=np.sin(l_cp-l)*(np.sin(de_gp)*np.sin(de)+np.cos(de_gp)*np.cos(de)*np.cos(al-al_gp))/(np.cos(de)*np.sin(al-al_gp))
   b=np.degrees(np.arctan(temp))
   l=np.degrees(l)%360.0
   
   return [l,b]


def stat(x):
    return [np.min(x),np.max(x),np.mean(x),np.sqrt(np.var(x))]


def uniqueid(x):
    ind=np.argsort(x)
    x=x[ind]
    flag=np.zeros(x.size,dtype='int32')
    for i in range(1,x.size):
        if x[i] == x[i-1]:
            flag[i]=1
    ind1=np.where(flag == 0)[0]
    return ind[ind1]


def distmod(r):
    return 5.0*np.log10(100.0*r)

def jk2vmag(j,ks):
    return ks+2.0*((j-ks)+0.14)+0.382*np.exp((j-ks-0.2)*2)

def equ2hor(ra,dec,observer,date):
    observer.date=date
    if (hasattr(ra,"__len__") == False):
        g=ephem_object(ra,dec)
        g.compute(observer)
        lon=np.degrees(g.az.znorm)
        lat=np.degrees(g.alt.znorm)
    else:
        ra=np.array(ra)
        dec=np.array(dec)
        lon=ra.copy()
        lat=dec.copy()
        for i in range(len(ra)):
            g=ephem_object(ra[i],dec[i])
            g.compute(observer)
            lon[i]=np.degrees(g.az.znorm)
            lat[i]=np.degrees(g.alt.znorm)        
    return [lon,lat]

def rad(x,y,z):
    return np.sqrt(x*x+y*y+z*z)

# def gal2equ(lon,lat,epoch='2000'):
#     g=ephem.Galactic(0.0,0.0,epoch=str(epoch))
#     if (hasattr(lon,"__len__") == False):
#         lon=np.array([lon])
#         lat=np.array([lat])
#         status=True
#     else:
#         status=False 
#     if (type(lon) in (tuple,list)):
#         lon=np.array(lon)
#         lat=np.array(lat)
        
#     ra=lon.copy()
#     dec=lat.copy()
#     for i in range(len(ra)):
#         g.set(np.radians(lon[i]),np.radians(lat[i]))
#         ra[i],dec[i]=np.degrees(g.to_radec())
#     if status:
#         ra=ra[0]
#         dec=dec[0]
#     return [ra,dec]

# def ecl2equ(lon,lat,epoch='2000'):
#     g=ephem.Ecliptic(0.0,0.0,epoch=str(epoch))
#     if (hasattr(lon,"__len__") == False):
#         lon=np.array([lon])
#         lat=np.array([lat])
#         status=True
#     else:
#         status=False    
#     if type(lon) in (tuple,list):
#         lon=np.array(lon)
#         lat=np.array(lat)
#     ra=lon.copy()
#     dec=lat.copy()
#     for i in range(len(ra)):
#         g.set(np.radians(lon[i]),np.radians(lat[i]))
#         ra[i],dec[i]=np.degrees(g.to_radec())
#     if status:
#         ra=ra[0]
#         dec=dec[0]
#     return [ra,dec]

# def equ2gal(ra,dec,epoch='2000'):
#     g=ephem.Galactic(0.0,0.0,epoch=str(epoch))
#     if (hasattr(ra,"__len__") == False):
#         ra=np.array([ra])
#         dec=np.array([dec])
#         status=True
#     else:
#         status=False    
#     if type(ra) in (tuple,list):
#         print 'here'
#         ra=np.array(ra)
#         dec=np.array(dec)
#     lon=ra.copy()
#     lat=dec.copy()
#     for i in range(len(ra)):
#         g.from_radec(np.radians(ra[i]),np.radians(dec[i]))
#         lon[i]=np.degrees(g.lon.znorm)
#         lat[i]=np.degrees(g.lat.znorm)
#     if status:
#         lon=lon[0]
#         lat=lat[0]
#     return [lon,lat]


# def equ2ecl(ra,dec,epoch='2000'):
#     g=ephem.Ecliptic(0.0,0.0,epoch=str(epoch))
#     if (hasattr(ra,"__len__") == False):
#         ra=np.array([ra])
#         dec=np.array([dec])
#         status=True
#     else:
#         status=False    
#     if type(ra) in (tuple,list):
#         ra=np.array(ra)
#         dec=np.array(dec)
#     lon=ra.copy()
#     lat=dec.copy()
#     for i in range(len(ra)):
#         g.from_radec(np.radians(ra[i]),np.radians(dec[i]))
#         lon[i]=np.degrees(g.lon.znorm)
#         lat[i]=np.degrees(g.lat.znorm)
#     if status:
#         lon=lon[0]
#         lat=lat[0]
#     return [lon,lat]


def histogram(a, **kwargs):
    [h,e]=np.histogram(a,**kwargs)
    d=np.digitize(a,e)
    ind=np.where((d>0)&(d<e.size))[0]
    ind=ind[np.argsort(d[ind])]
    ri=np.zeros(e.size+ind.size,dtype='int64')
    ri[0]=e.size        
    ri[1:e.size]=np.cumsum(h)+e.size        
    ri[e.size:]=ind
    return h,e,ri

def hist_indices(ri,i):
    if i<ri[0]:
        return ri[ri[i]:ri[i+1]]
    else:
        return np.array([],dtype='int64')

def create_ri(d,bins=None):
    if bins == None:
        bins=np.max(d)+1
    h,e,ri=histogram(d,bins=bins,range=[0,bins])
    return ri


def spring_next_transit(date,ra,dec):
    s=sspring(date)
    b=ephem_object(ra,dec)
    return s.next_transit(b)



def cross_match(src,quant,suffix='',prefix=''):
    mydict=dict(zip(src,np.arange(src.size)))
    ind=np.zeros(quant.size,dtype='int64')
    if len(suffix+prefix)>0:
        for i,temp in enumerate(quant):
            ind[i]=mydict.get(prefix+temp+suffix,-1)
    else:
        for i,temp in enumerate(quant):
            ind[i]=mydict.get(temp,-1)
    return ind

def tocsv(data,filename,basekey=None,keylist=None,delimiter=', '):
    if type(data) == dict:
        with open(filename,'w') as fp:
            if keylist==None:
                keylist=data.keys()
            if basekey == None:
                nsize=data[keylist[0]].size
            else:
                nsize=data[basekey].size        
            keylist=[key for key in keylist if data[key].size==nsize] 
            # s=str(keylist)
            # s=s[1:-1].replace("'",'')+'\n'
            s=delimiter.join([str(key) for key in keylist])+'\n'
            fp.write(s)
            for i in range(data[keylist[0]].size):
                s=', '.join([str(data[key][i]) for key in keylist])+'\n'
                fp.write(s)
    else:
        with open(filename,'w') as fp:
            if keylist==None:
                # s=str(data.dtype.names)
                # s=s[1:-1].replace("'",'')+'\n'
                s=delimiter.join([str(key) for key in data.dtype.names])+'\n'
                fp.write(s)
                for temp in data:
                    s=delimiter.join([str(temp[key]) for key in data.dtype.names])+'\n'
                    fp.write(s)
                # for temp in data:
                #     s=str(temp)
                #     s=s[1:-1].replace("'",'')+'\n'
                #     fp.write(s)
            else:
                # s=str(keylist)
                # s=s[1:-1].replace("'",'')+'\n'
                s=delimiter.join([str(key) for key in keylist])+'\n'
                fp.write(s)
                for i in range(data[keylist[0]].size):
                    s=delimiter.join([str(data[key][i]) for key in keylist])+'\n'
                    fp.write(s)
    print 'Written file:',filename



def ebfread(file1,keys):
#    keys1=[path+key for key in keys]
    data={}
    for key in keys:
        data[key.rsplit('/',1)[1]]=ebf.read(file1,key)
    return data

# def advance_radec(ra,dec,pmra,pmdec,years,actual=True):
#     """
#     Advance the RA, Dec by the given proper motions
#     Parameters
#     ----------
#     ra:    degrees
#     dec:   degrees
#     pmra:  degrees/year
#     pmdec: degress/year
#     years: year
#     actual: If True then it means pmra= cos(DEC)*dRA/dt 
#     """
#     if actual:
#         ra=(ra+pmra*advance/cos(np.radians(dec)))%360.0
#     else:
#         ra=(ra+pmra*advance)%360.0
#     dec=dec+pmdec*advance
#     if dec > 90.0:
#         dec=90-(dec-90.0)
#     elif dec < -90.0:
#         dec=-90-(dec+90.0)
#     return [ra,dec]

# def gal2ecl1(l, b,epoch='2000'):   
#    # al_gp = np.radians(192.85948) 
#    # de_gp = np.radians(27.12825)
#    # l_cp = np.radians(122.932e0) 
#    g=ephem.Galactic(0.0,np.pi/2.0,epoch=str(epoch))
#    [al_gp,de_gp]=g.to_radec()

#    temp=ephem.Ecliptic(0.0,0.0,epoch=str(epoch))
#    ephem.from_radec(al_gp,de_gp,epoch=str(epoch))
   
#    g.from_radec(0.0,np.pi/2.0)
#    l_cp=g.lon.znorm
   
#    l = np.radians(l) 
#    b = np.radians(b)
   
#    al = al_gp + np.arctan2(np.cos(b) * np.sin(l_cp - l), np.cos(de_gp) * np.sin(b) - np.sin(de_gp) * np.cos(b) * np.cos(l_cp - l))   
#    temp = np.sin(al - al_gp) * (np.sin(de_gp) * np.sin(b) + np.cos(de_gp) * np.cos(b) * np.cos(l_cp - l)) / (np.cos(b) * np.sin(l_cp - l))   
#    de = np.degrees(np.arctan(temp))   
#    al = np.degrees(al)%360.0    
#    al=al%360.0

#    return [al,de]

def read(file1,schema=None,catalog=None,key=None,quant=None):
    if file1.endswith('.ebf'):
        if ebf.containsKey(file1,'/data'):
            data1=ebf.read(file1,'/data')
            data1=npstruct2dict(data1)
        else:
            data1=ebf.read(file1,'/')
    elif file1.endswith('.fit') or file1.endswith('.fits'):
        data1=fitsread(file1)
    else:
        data1=asctab.read(file1)

    if catalog != None:
        schema='/work1/sharma/dbm/schemas/cdef_'+catalog+'.txt'
    if schema != None:
        l=asctab.read(schema)
        if 'alias' in l.keys():
            l['ukey']=l['alias']
        s=set(data1.keys())
        data2={}
        for i in range(len(l['key'])):
            if l['key'][i] in s:
                data2[l['ukey'][i]]=data1[l['key'][i]]
        data1=data2
    if key != None:
        if quant != None: 
            return filter_dict(data1,cross_match(data1['key'],quant))
        else:
            raise RuntimeError('FOr cross match, quant should not be None')
    else:
        return data1



def fitsread(filename,ext=1):
    import pyfits
    data1=np.array(pyfits.getdata(filename,ext))
    data={}
    for x in data1.dtype.descr:
        data[x[0].lower()]=data1[x[0]]
    return data


def str2num(x,mytype=np.float64):
    status=0
    try:        
        y=np.array(x,dtype=mytype)
    except ValueError:
        x=np.array(x)    
        if x.ndim==0:
            status=1
            x=x.reshape((1,))
        if int(x.dtype.descr[0][1].rsplit('S')[1]) < 3:
            x=np.array(x,dtype='S3')                    
        x=np.char.lower(np.char.strip(x))
        x[np.where((x=='')|(x=='null'))[0]]='nan'
        try:
            y=np.array(x,dtype=mytype)            
        except ValueError:        
            y=np.zeros(x.size,dtype=mytype)            
            if (y.dtype.char != 'f')|(y.dtype.char != 'd'):
                y1=np.array(x,dtype=np.float64)
                ind=np.where(np.isfinite(y1)==True)[0]
                y[ind]=mytype(y1[ind])
                ind=np.where(np.isfinite(y1)==False)[0]
                y[ind]=np.iinfo(y.dtype).max

    if status==1:
        return y[0]
    else:
        return y


def dict2npstruct(data,basekey=None,keylist=None,sortby=None,dtype=None):    
    if keylist==None:
        keylist=data.keys()
    if basekey == None:
        nsize=data[keylist[0]].size
    else:
        nsize=data[basekey].size        

    if dtype==None:    
        dt=[]
        for key in keylist:
            if data[key].size == nsize:
                dt.append((key,data[key].dtype))
    else:
        dt=dtype

    data1=None
    if len(dt)>0:
        data1=np.zeros(nsize,dtype=dt)
        print data.keys()
        for key in data1.dtype.names:
            data1[key.lower()]=data[key]

    return data1

# def dict2npstruct(data,basekey=None,keylist=None):    
#     if keylist==None:
#         keylist=data.keys()
#     if basekey == None:
#         nsize=data[keylist[0]].size
#     else:
#         nsize=data[basekey].size        

#     dt=[]
#     for key in keylist:
#         if data[key].size == nsize:
#             dt.append((key,data[key].dtype))

#     data1=None
#     if len(dt)>0:
#         data1=np.zeros(nsize,dtype=dt)
#         for key in data1.dtype.names:
#             data1[key.lower()]=data[key]

#     return data1


def cross_match_angular(tree,ra,dec,radius):
    index=np.zeros(ra.size,dtype='int64')-1
    dist=np.zeros(ra.size,dtype='float64')
    for i in range(ra.size):
        [ind1,r,ind2]=tree.ngbsearch_angular(ra[i],dec[i],radius)
        if ind2.size>0:            
            index[i]=ind2[0]
            dist[i]=r[0]
    return index,dist


def filter_dict(data,ind,basekey=None,keylist=None,aliaslist=None):
    if keylist==None:
        if type(data)==dict:
            keylist=data.keys()
        else:
            keylist=data.dtype.names
    if aliaslist==None:
        aliaslist=keylist

    ind1=np.array(ind)
    ind1=ind1[np.where(ind>=0)[0]]

    if basekey == None:
        nsize=data[keylist[0]].size
    else:
        nsize=data[basekey].size        
    data1={}
    for i in range(len(keylist)):
#        print key,data[key].size,nsize
        if data[keylist[i]].size == nsize:
            if ind == None:
                data1[aliaslist[i]]=data[keylist[i]]
            else:
                data1[aliaslist[i]]=data[keylist[i]][ind1]
    return data1


def undef_values(dtype):
    if dtype != np.dtype:
        dtype=np.dtype(dtype)
    t=''
    if dtype.kind == 'i':
        t=np.iinfo(dtype).max
    elif dtype.kind == 'u':
        t=np.iinfo(dtype).max
    elif dtype.kind == 'f':
        t=np.inf
    return t

def print_utd(data,catalog=None):        
    if type(data) == str:
        file1=data
        if (file1.endswith('.fits'))or(file1.endswith('.fit')): 
            import pyfits
            data=pyfits.getdata(file1,1)
        elif file1.endswith('.ebf'):
            if ebf.containsKey(file1,'/data'):
                data=ebf.read_ind(file1,'/data',[0,1,2])
            else:
                data={}
                for key in ebf.keys(file1,'/'):
                    data[key]=ebf.read_ind(file1,'/'+key,[0,1,2])
        else:
            data=asctab.read(file1,struct=True)

    if type(data) == dict:
        keys=data.keys()
    else:
        keys=data.dtype.names
    

    if catalog != None:
        if catalog != None:
            catfile='/work1/sharma/dbm/schemas/cdef_'+catalog+'.txt'
            # if os.path.isfile(catfile):
            #     raise RuntimeError('Catalog definition aleady exists')
            f=open(catfile,'w')
            f.write('{0:16s}, {1:8s}, {2:4s}\n'.format('key','datatype','ukey'))
            for key in keys:
                f.write('{0:16s}, {1:8s}, {2:16s}\n'.format(key.lower(),data[key].dtype.str[1:],catalog+'_'+key.lower()))
            f.close()
    else:
        print '{0:16s}, {1:8s}'.format('key','datatype')
        for key in keys:
            print '{0:16s}, {1:8s}'.format(key.lower(),data[key].dtype.str[1:])



def t_size(data1):
    return data1[data1.keys()[0]].size

def t_addcol(data1,key,dt,value=None):
    nsize=data1[data1.keys()[0]].size
    data1[key]=np.zeros(nsize,dtype=dt)
    if value != None:
        data1[key][:]=value
    

def t_unique(data1,pkey):
    x,ind=np.unqiue(data1['pkey'],return_index=True)
    if ind.size<data1[pkey].size:
        t_filter(data1,ind)

def t_where(data1,condition):
    t_filter(data1,np.where(condition)[0])

def t_filter(data1,ind):
    for key in data1.keys():
        data1[key]=data1[key][ind]
    

def t_union(data1,data2,pkey,cols_or=None,cols_update=None):
    indl,indr,imdrnm=autil.cross_match(data1[pkey],data2[pkey],split=True)
    for key in cols_or:
        data1[key][indl]=data1[key][indl]|data2[key][indr]
    for key in cols_update:
        data1[key][indl]=data2[key][indr]
    for key in data1.keys():
        data1[key]=np.concantenate([data1[key],data2[key][indrnm]])
    

def t_ljoin(data1,data2,pkey1,pkey2=None):
    if pkey2 == None:
        pkey2=pkey1

    indl,indr,indrnm=autil.cross_match(data2[pkey2],data1[pkey1],split=True)

    if type(data2) == dict:
        s2=data2.keys()
    else:
        s2=data2.dtype.names
    if type(data1) == dict:
        s1=set(data1.keys())
    else:
        s1=set()

    for tkey in s2:
        if tkey not in s1:
            temp=np.zeros(data[pkey[0]].size,dtype=datat[tkey].dtype)
            temp[indrnm]=undef_values(temp.dtype)
            temp[indr]=datat[tkey][indl]
            data1[tkey]=temp
                
