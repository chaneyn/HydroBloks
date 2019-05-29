import numpy as np
import matplotlib.pyplot as plt
import time
import netCDF4 as nc
import scipy.sparse as sparse
import scipy.sparse.linalg
import pickle
import numba
import sys

@numba.jit(nopython=True,cache=True)
def calculate_hydraulic_radius(A,P,W,A1):

 Rh = np.zeros(A.shape[0])
 for i in range(A.shape[0]):
  A0 = 0.0
  P0 = 0.0
  for j in range(A.shape[1]):
   W1 = W[i,j]
   if A[i,j] > A1[i]:break
   if A[i,j+1] == 0.0:break
   A0 = A[i,j]
   P0 = P[i,j]
  #Calculate the height above the segment
  h = (A1[i] - A0)/W1
  #Calculate the wetted perimeter
  P1 = P0 + 2*h + W1
  #Calculate hydraulic radius
  Rh[i] = A1[i]/P1

 return Rh

#cid = 7100454060
#cid = 7100456120
#cid = 7100454070
cid = 7100456120
#Read in the stream network information
file = '/home/nc153/soteria/projects/hydroblocks_inter_catchment/run_marion_county/input_data/catchments/%d/input_file_enhanced_test.nc' % cid
fp = nc.Dataset(file)
grp = fp['stream_network']
dbc = {}
for var in grp.variables:
 dbc[var] = grp[var][:]
ucids = np.unique(dbc['cid'])
fp.close()

#Read in the other info
odbc = {}
for ucid in ucids:
 file = '/home/nc153/soteria/projects/hydroblocks_inter_catchment/run_marion_county/input_data/catchments/%d/input_file_enhanced.nc' % ucid
 fp = nc.Dataset(file)
 grp = fp['stream_network']
 odbc[ucid] = {}
 for var in ['channelid','cid']:
  odbc[ucid][var] = grp[var][:]

#Read in the boundary condition data to assimilate
bcdata = {}
for ucid in ucids:
 file = 'workspace/%s.pck' % ucid
 bcdata[ucid] = pickle.load(open(file,'rb'))

#Read in the runoff output
file = '/home/nc153/soteria/projects/hydroblocks_inter_catchment/run_marion_county/output_data/%s/output/2002-01-01.nc' % cid
fp = nc.Dataset(file)
runoff = fp['data']['runoff'][:]

#Read in reach/hru area
file = '../../../hydroblocks_inter_catchment/run_marion_county/input_data/catchments/%d/routing_info.pck' % cid
db = pickle.load(open(file,'rb'))['reach_hru_area']

#Read in the reach/hand database and remap to new topology (THERE COULD BE AN ORDERING PROBLEM)
hdb = {}
for ucid in ucids:
 file = '/home/nc153/soteria/projects/hydroblocks_inter_catchment/run_marion_county/input_data/catchments/%d/routing_info.pck' % ucid
 ohdb = pickle.load(open(file,'rb'))['reach_cross_section']
 for var in ohdb:
  if var not in hdb:hdb[var] = np.zeros((dbc['topology'].size,ohdb[var].shape[1]))
  m = dbc['cid'] == ucid
  hdb[var][m,:] = ohdb[var][dbc['channelid'][m],:]

#Create a matrix
reach2hru = np.zeros((np.sum(odbc[cid]['cid']==cid),runoff.shape[1]))
for reach in db:
 for hru in db[reach]:
  reach2hru[reach-1,hru] = db[reach][hru]
  #reach2hru[reach-1,hru] = db[reach][hru]
reach2hru = sparse.csr_matrix(reach2hru)
area = np.sum(reach2hru,axis=0)

#Construct mapping of other catchment network to that of the current network
mapping_ucid = {}
for ucid in ucids:
 mapping_ucid[ucid] = {'cid':[],'ocid':[]}
 idxs = np.where(dbc['cid'] == ucid)[0]
 for idx in idxs:
  m = (odbc[ucid]['cid'] == ucid) & (odbc[ucid]['channelid'] == dbc['channelid'][idx])
  mapping_ucid[ucid]['cid'].append(idx)
  mapping_ucid[ucid]['ocid'].append(np.where(m)[0][0])
 #Convert to arrays
 for var in mapping_ucid[ucid]:
  mapping_ucid[ucid][var] = np.array(mapping_ucid[ucid][var])

#Assemble the connectivity array (Best to define this reordering in the database creation)
corg = np.arange(dbc['topology'].size)
cdst = dbc['topology'][:]
m = cdst != -1
nc = cdst.size
cmatrix = sparse.coo_matrix((np.ones(cdst[m].size),(corg[m],cdst[m])),shape=(nc,nc),dtype=np.float32)
cmatrix = cmatrix.tocsr().T

#Assemble mapping
mapping = np.zeros(np.sum(dbc['cid']==cid)).astype(np.int64)
for i in range(mapping.size):
 mapping[i] = np.where((dbc['cid'] == cid) & (dbc['channelid'] == i))[0][0]

#Identify headwaters
cup = -1*np.ones(cdst.size)
for id in corg:
 n = np.sum(cdst == id)
 if n > 0:cup[id] = 1

#import copy
c_length = dbc['length'][:]
c_slope = dbc['slope'][:]
c_width = dbc['width'][:]
c_n = dbc['manning'][:]
Ainit = np.zeros(c_length.size)
A0 = np.copy(Ainit)
A1 = np.copy(Ainit)

#Filler
dt = 10800 #s
#tmax = 100*3600*24
tmax = dt*250#runoff.shape[0]
#maxc = 10.0 #m/s
nt = int(tmax/dt)
#Define initial conditions
Qinit = np.zeros(nc)
Q0 = Qinit[:]
qin = np.zeros(c_length.size)
qout = np.zeros(c_length.size)

out = {'Q':[],'A':[],'qin':[],'qout':[]}
dif0 = -9999
max_niter = 10
for t in range(nt):
 print(cid,t)
 A0_org = np.copy(A0)
 #Q0_org = np.copy(Q0)
 qin[:] = 0.0
 #Update initial conditions using upstream information
 if t > 0:
  for ucid in ucids:
   if ucid == cid:continue
   A0_org[mapping_ucid[ucid]['cid']] = bcdata[ucid]['A'][t-1,mapping_ucid[ucid]['ocid']]
 #Update the lateral inputs/outputs
 for ucid in ucids:
  if ucid == cid:continue
  qin[mapping_ucid[ucid]['cid']] = bcdata[ucid]['qin'][t,mapping_ucid[ucid]['ocid']]

 for it in range(max_niter):
  #Compute inflows
  qin[mapping] = reach2hru.dot(runoff[t,:]/1000.0/dt)/c_length[mapping] #m/s
  qin[qin < 0] = 0.0
  #Determine hydraulic radius
  rh = calculate_hydraulic_radius(hdb['A'],hdb['P'],hdb['W'],A0)
  #Determine velocity
  u = rh**(2.0/3.0)*c_slope**0.5/c_n
  #Fill non-diagonals
  LHS = cmatrix.multiply(-dt*u)
  #Fill diagonal
  LHS.setdiag(c_length + dt*u)
  #Set right hand side
  RHS = c_length*A0_org + dt*qin*c_length - dt*qout*c_length
  #Ax = b
  A1 = scipy.sparse.linalg.spsolve(LHS,RHS,use_umfpack=True)
  #QC
  A0[A0 < 0] = 0.0
  dif1 = np.mean(np.abs(A0 - A1))
  if (dif1 < 10**-10) | (it == max_niter-1):
   #Reset A0
   A0[:] = A1[:]
   h = A0/c_length #rectangular channel
   #Determine hydraulic radius
   rh = calculate_hydraulic_radius(hdb['A'],hdb['P'],hdb['W'],A0)
   #Determine velocity
   u = rh**(2.0/3.0)*c_slope**0.5/c_n
   #Calculate Q1
   Q1 = A0*u
   dif0 = -9999
   break
  else:
   #Reset A0
   A0[:] = A1[:]
   dif0 = dif1 

 #Append to output
 out['Q'].append(np.copy(Q1))
 out['A'].append(np.copy(A1))
 out['qin'].append(np.copy(qin))
 out['qout'].append(np.copy(qout))
for var in out:
 out[var] = np.array(out[var])
m = cdst == -1
#print(np.sum(np.trapz(out['Q'][:,m].T,dx=dt)),"m3")
dVh = -np.sum(c_length*np.diff(out['A'],axis=0),axis=1)
dVh += np.sum(c_length*dt*out['qin'],axis=1)[1:]
dVh -= np.sum(c_length*dt*out['qout'],axis=1)[1:]
dVQ = dt*np.sum(out['Q'][1:,m],axis=1)
print(np.sum(dVh),np.sum(dVQ))
print(np.sum(dVQ)/np.sum(area))
#plt.plot(out['Q'][:,dbc['cid']==cid])
#plt.plot(bcdata[cid]['Q'][:,odbc[cid]['cid']==cid])
#A0_org[mapping_ucid[ucid]['cid']] = bcdata[ucid]['A'][t-1,mapping_ucid[ucid]['ocid']]
Q1 = out['Q'][:,mapping_ucid[cid]['cid']]
Q0 = bcdata[cid]['Q'][:,mapping_ucid[cid]['ocid']]
maxe = 0
for b in range(Q1.shape[1]):
 if 100*np.mean(np.abs(Q1[:,b]-Q0[:,b]))/(np.max(Q0[:,b]) - np.min(Q0[:,b])) > maxe:
  maxe = 100*np.mean(np.abs(Q1[:,b]-Q0[:,b]))/(np.max(Q0[:,b]) - np.min(Q0[:,b]))
 print('b',b,100*np.mean(np.abs(Q1[:,b]-Q0[:,b]))/(np.max(Q0[:,b]) - np.min(Q0[:,b])))
print(maxe)
plt.plot(Q1[:,0])
plt.plot(Q0[:,0])
plt.show()
exit()
plt.subplot(121)
plt.plot(out['Q'][:,:])
#plt.plot(np.log10(out['Q'][:,0]))
plt.legend(['C1','C2'])
plt.subplot(122)
plt.plot(out['A'][:,:])
#plt.plot(out['A'][:,0])
plt.show()
#pickle.dump(out,open('workspace/%s.pck' % cid,'wb'))
