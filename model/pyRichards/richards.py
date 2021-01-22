import numpy as np
import scipy.sparse as sparse
import time
import numba

class richards:

 def __init__(self,nhru,nhru2,nsoil): #laura

  self.theta = np.zeros((nhru2,nsoil))
  self.thetar = np.zeros(nhru2)
  self.thetas = np.zeros(nhru2)
  self.b = np.zeros(nhru2)
  self.satpsi = np.zeros(nhru2)
  self.ksat = np.zeros(nhru2)
  self.dem = np.zeros(nhru)
  self.demhband = np.zeros(nhru2) #laura added
  self.dem1hband=np.zeros(nhru2) #laura added
  self.slope = np.zeros(nhru2)
  #self.hand = np.zeros(nhru)
  self.area = np.zeros(nhru2)
  self.dz = np.zeros((nhru2,nsoil))
  self.hdiv = np.zeros((nhru2,nsoil))
  self.m = np.zeros(nhru2)

  #Initialize the width array
  self.width = {}#[] #laura
  self.I = {}#[]#laura

  return

 def calculate_soil_moisture_potential(self,il):
  
  eps = 0.01
  theta = self.theta[:,il]
  thetar = self.thetar
  thetas = self.thetas
  b = self.b
  satpsi = self.satpsi #meters
  m = (theta <= (1+eps)*thetar)
  theta[m] = (1+eps)*thetar[m]
  #psi = 1000.0*np.ones(theta.shape) #limit
  #m = theta > thetar
  #psi[m] = satpsi[m]*((theta[m] - thetar[m])/(thetas[m] - thetar[m]))**(-b[m])
  #psi[m] = satpsi[m]*((theta[m] - thetar[m])/(thetas[m] - thetar[m]))**(-b[m])

  with np.errstate(invalid='ignore'):
    psi = satpsi*((theta-thetar)/(thetas-thetar))**-b

  return psi

 def calculate_hydraulic_conductivity(self,psi,il):
 
  af = 1.0 #safe
  #sdz = np.cumsum(self.dz,axis=1)-self.dz/2.
  #df = np.exp(-self.m[:,np.newaxis]/sdz)[:,il]
  #Ksat_x = af*df*self.ksat[:] #lateral saturated hydraulic conductivity (multiply times anisotropy factor) [m/s]
  Ksat_x = af*self.ksat[:] #lateral saturated hydraulic conductivity (multiply times anisotropy factor) [m/s]
  with np.errstate(invalid='ignore'):
   K_x = Ksat_x*(psi/self.satpsi)**(-2-3/self.b)

  return K_x

 def calculate_transmissivity(self,psi,ztop,zbot):
  
  #af = 1.0  #safe
  af = 2.0
  m = np.copy(self.m)
  #m[:] = 1000.0
  Ksat_x = af*self.ksat[:] #lateral saturated hydraulic conductivity (multiply times anisotropy factor) [m/s]
  with np.errstate(invalid='ignore', divide='ignore'):
   K_x = Ksat_x*np.true_divide(psi,self.satpsi)**(-2-np.true_divide(3.,self.b))
   K_x[~np.isfinite(K_x)] = np.nan
   #Calculate transmissivity at top layer (exponential decay)
   Ttop = m*K_x*np.exp(-ztop/m)
   #Calculate transmissivity at bottom of layer (exponential decay)
   Tbot = m*K_x*np.exp(-zbot/m)
   T = Ttop - Tbot
    
  return T

 def calculate_hydraulic_head(self,psi,depth):
  
  h = self.dem - depth - psi

  return h

 #def calculate_divergence_dense(self,h,K_x,dz):
 #@numba.jit(nopython=True,cache=True)
 def calculate_divergence_dense(self,h,T):
 
  #tic = time.time()
  dh = h[:,np.newaxis] - h[np.newaxis,:]
  w = self.w
  dx = self.dx
  area = self.area 
  That = np.true_divide((2*T[:,np.newaxis]*T[np.newaxis,:]),(T[:,np.newaxis] + T[np.newaxis,:]))
  That[~np.isfinite(That)] = np.nan
  #[mm/s] = [mm/m]*[m/s]*[m]/[m]*[m]*[m]/[m2]
  calc_div = -1000.0*That*np.true_divide(dh,dx)*np.true_divide(w,area) # mm/s
  calc_div[~np.isfinite(calc_div)] = np.nan
  #print('calc_div',time.time() - tic)

  return calc_div

 #def calculate_divergence_sparse(self,h,K_x,dz):
 def calculate_divergence_sparse(self,h,T):

  #Define the boolean matrix (connections or not?)
  I = self.I
  #Calculate dh
  h1 = (I != 0).multiply(sparse.csr_matrix(h))
  dh = h1.T - h1
  #Calculate dx
  d1 = (I != 0).multiply(sparse.csr_matrix(self.dem))
  dx = d1.T - d1#**2 + self.dx**2)
  dx = dx.power(2)
  dx.data += self.dx**2
  dx = dx.power(0.5)
  #dx = (np.abs(d1.T - d1)**2 + self.dx**2)**0.5
  #dx = (np.abs(d1.T - d1)**2 + self.dx**2)**0.5
  #dx = (np.abs(self.dem[:,np.newaxis] - self.dem[np.newaxis,:])**2 + self.dx**2)**0.5
  #Calculate the effective hydraulic conductivity 
  #k1 = (I != 0).multiply(sparse.csr_matrix(K_x))
  #n = 2*k1.T.multiply(k1)
  #d = k1.T+k1
  #Khat = n.multiply(d.power(-1))
  t1 = (I != 0).multiply(sparse.csr_matrix(T))
  n = 2*t1.T.multiply(t1)
  d = t1.T+t1
  That = n.multiply(d.power(-1)).tocsr()
  print(That.count_nonzero,self.width.count_nonzero)
  #Calculate the flux
  #[m/s] = [m/s]*[m]/[m]*[m]*[m]/[m2]
  print(That.multiply(dh).shape)
  print(That.multiply(dh).multiply(self.width).count_nonzero)
  print(1.0/self.area)
  print(That.multiply(dh).multiply(self.width).multiply(1.0/self.area).count_nonzero)
  print(That.multiply(dh).multiply(self.width).multiply(1.0/self.area).multiply(dx.power(-1)).count_nonzero)
  return -That.multiply(dh).multiply(self.width).multiply(1.0/self.area).multiply(dx.power(-1)).multiply(1000) #mm/s
  #return -Khat.multiply(dh).multiply(self.width).multiply(dz/self.dx/self.area).multiply(1000) #mm/s

 def update(self):

  #Iterate per layer
  for il in range(self.theta.shape[1]):
   #Calculate soil moisture potential
   psi = self.calculate_soil_moisture_potential(il)
   zbot = np.sum(self.dz[:,0:il+1],axis=1)
   ztop = zbot - self.dz[:,il]
   T = self.calculate_transmissivity(psi,ztop,zbot)
   #Calculate hydraulic head
   h = self.calculate_hydraulic_head(psi,ztop)
   #Calculate the divergence
   q = self.calculate_divergence_dense(h,T)
   self.hdiv[:,il] = np.sum(q,axis=0) #mm/s  

  return

 def update_numba(self):

  theta = self.theta
  dz = self.dz
  hdiv = self.hdiv
  thetar = self.thetar
  thetas = self.thetas
  b = self.b
  satpsi = self.satpsi
  m = self.m
  ksat = self.ksat
  #hand = self.dem
  hand = self.dem1hband
  w = self.w
  dx = self.dx
  area = self.area
  ncsbasins=self.ncsbasins #laura, number of characteristic subbasins
  #self.hdiv[:] = update_workhorse(theta,dz,hdiv,thetar,thetas,b,satpsi,m,ksat,hand,w,dx,area)
  #Compute divergence independently per characteristic subbasin, laura
  aux=0
  div=np.empty(self.hdiv.shape)
  for bas in range(1,ncsbasins+1):
   w_bas=w['Basin%s' %bas]
   dx_bas=dx['Basin%s' %bas]
   init=aux
   fin=aux+w_bas.shape[0]
   div[init:fin,:]=update_workhorse(theta[init:fin,:],dz[init:fin,:],hdiv[init:fin,:],thetar[init:fin],thetas[init:fin],b[init:fin],satpsi[init:fin],m[init:fin],ksat[init:fin],hand[init:fin],w_bas,dx_bas,area[init:fin])
   aux=fin
  self.hdiv=div

  return
 
@numba.jit(nopython=True,cache=True)
def update_workhorse(theta,dz,hdiv,thetar,thetas,b,satpsi,m,ksat,hand,w,dx,area):

 #Iterate per layer
 for il in range(theta.shape[1]):
  #Calculate soil moisture potential
  psi = calculate_soil_moisture_potential(il,theta,thetar,thetas,b,satpsi)
  zbot = np.sum(dz[:,0:il+1],axis=1)
  ztop = zbot - dz[:,il]
  T = calculate_transmissivity(psi,ztop,zbot,m,ksat,satpsi,b)
  #Calculate hydraulic head
  h = calculate_hydraulic_head(hand,psi,ztop)
  #Calculate the divergence
  q = calculate_divergence(h,T,w,dx,area)
  hdiv[:,il] = np.sum(q,axis=0) #mm/s'''

 return hdiv

@numba.jit(nopython=True,cache=True)
def calculate_soil_moisture_potential(il,theta,thetar,thetas,b,satpsi):
  
 eps = 0.01
 theta = theta[:,il]
 m = (theta <= (1+eps)*thetar)
 theta[m] = (1+eps)*thetar[m]
 psi = satpsi*((theta-thetar)/(thetas-thetar))**-b

 return psi

@numba.jit(nopython=True,cache=True)
def calculate_transmissivity(psi,ztop,zbot,m,ksat,satpsi,b):
  
 af = 1.0#10.0#2.0
 Ksat_x = af*ksat #lateral saturated hydraulic conductivity (multiply times anisotropy factor) [m/s]
 K_x = Ksat_x*np.true_divide(psi,satpsi)**(-2-np.true_divide(3.,b))
 #Calculate transmissivity at top layer (exponential decay)
 Ttop = m*K_x*np.exp(-ztop/m)
 #Calculate transmissivity at bottom of layer (exponential decay)
 Tbot = m*K_x*np.exp(-zbot/m)
 T = Ttop - Tbot
  
 return T

@numba.jit(nopython=True,cache=True)
def calculate_hydraulic_head(hand,psi,depth):
  
 h = hand - depth - psi

 return h

@numba.jit(nopython=True,cache=True)
def calculate_divergence(h,T,w,dx,area):
 
 #Calculate dh
 #dh = h[:,np.newaxis] - h[np.newaxis,:]
 dh = calculate_dh(h)
 #Calculate That
 #That = np.true_divide((2*T[:,np.newaxis]*T[np.newaxis,:]),(T[:,np.newaxis] + T[np.newaxis,:]))
 That = calculate_That(T)
 #That[~np.isfinite(That)] = np.nan
 #[mm/s] = [mm/m]*[m/s]*[m]/[m]*[m]*[m]/[m2]
 calc_div = -1000.0*That*dh/dx*w/area # mm/s
 #calc_div[~np.isfinite(calc_div)] = np.nan

 return calc_div

@numba.jit(nopython=True,cache=True)
def calculate_dh(h):

 dh = np.zeros((h.size,h.size))
 for i in range(h.size):
  for j in range(h.size):
   dh[i,j] = h[i] - h[j]
   #dh = h[:,np.newaxis] - h[np.newaxis,:]

 return dh

@numba.jit(nopython=True,cache=True)
def calculate_That(T):

 That = np.zeros((T.size,T.size))
 for i in range(T.size):
  for j in range(T.size):
   That[i,j] = (2*T[i]*T[j])/(T[i] + T[j])
   #That[i,j] = np.true_divide((2*T[:,np.newaxis]*T[np.newaxis,:]),(T[:,np.newaxis] + T[np.newaxis,:]))

 return That
