import os
import glob
import json
import netCDF4 as nc
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
import scipy.sparse as sparse
import scipy.sparse.linalg
import pickle
import numba
import time
from mpi4py import MPI

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

#Get some general info
dir = os.getcwd()

#Determine communication information
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()

#Define the catchment
rdir = '/home/nc153/soteria/projects/hydroblocks_inter_catchment/regions/ohio_basin'
tmp = glob.glob('%s/input_data/domain/*' % rdir)
tmp.remove('%s/input_data/domain/domain_database.pck' % rdir)
cids = []
for cid in tmp:
 cids.append(int(cid.split('/')[-1]))
cid = cids[rank]

#Determine catchments that rely on this catchment (ids to send to)
scids = []
for ucid in cids:
 if ucid == cid:continue
 file = '%s/input_data/domain/%d/input_file_enhanced_test.nc' % (rdir,ucid)
 fp = nc.Dataset(file)
 grp = fp['stream_network']
 ucids = np.unique(grp['cid'][:])
 fp.close()
 if cid in ucids:scids.append(ucid)

#Read in the stream network information
file = '%s/input_data/domain/%d/input_file_enhanced_test.nc' % (rdir,cid)
fp = nc.Dataset(file)
grp = fp['stream_network']
dbc = {}
for var in grp.variables:
 dbc[var] = grp[var][:]
ucids = np.unique(dbc['cid'])
fp.close()

#Define ids to receive from
rcids = list(ucids)
rcids.remove(cid)

#Determine the rank of all the cid
ranks = []
for ucid in ucids:
 if ucid == cid:continue
 ranks.append(cids.index(ucid))

#Read in the other info
odbc = {}
for ucid in ucids:
 file = '%s/input_data/domain/%d/input_file_enhanced_test.nc' % (rdir,ucid)
 fp = nc.Dataset(file)
 grp = fp['stream_network']
 odbc[ucid] = {}
 for var in ['channelid','cid']:
  odbc[ucid][var] = grp[var][:]

#Read in the runoff output
file = '%s/output_data/%s/output/2002-01-01.nc' % (rdir,cid)
fp = nc.Dataset(file)
runoff = fp['data']['runoff'][:]

#Read in reach/hru area
file = '%s/input_data/domain/%d/routing_info.pck' % (rdir,cid)
db = pickle.load(open(file,'rb'))['reach_hru_area']

#Read in the reach/hand database and remap to new topology (THERE COULD BE AN ORDERING PROBLEM)
hdb = {}
for ucid in ucids:
 file = '%s/input_data/domain/%d/routing_info.pck' % (rdir,ucid)
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
reach2hru = sparse.csr_matrix(reach2hru)
area = np.sum(reach2hru,axis=0)

#Bring in all the pertaining mapping_ucid
odb = {}
for ucid in scids:
 odb[ucid] = pickle.load(open('input/%s.pck' % ucid,'rb'))
for ucid in [cid,]:
 odb[ucid] = pickle.load(open('input/%s.pck' % ucid,'rb'))

#Assemble the connectivity array (Best to define this reordering in the database creation)
corg = np.arange(dbc['topology'].size)
cdst = dbc['topology'][:]
m = cdst != -1
nc = cdst.size
cmatrix = sparse.coo_matrix((np.ones(cdst[m].size),(corg[m],cdst[m])),shape=(nc,nc),dtype=np.float32)
cmatrix = cmatrix.tocsr().T.tocsr()

#Assemble mapping
mapping = np.zeros(np.sum(dbc['cid']==cid)).astype(np.int64)
for i in range(mapping.size):
 m = (dbc['cid'] == cid) & (dbc['channelid'] == i)
 mapping[i] = np.where(m)[0][0]

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
Ainit[:] = 0.001
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
stime = 0.0
rtime = 0.0
for t in range(nt):
 print(cid,t)
 A0_org = np.copy(A0)
 qin[:] = 0.0
 #Compute inflows
 qin[mapping] = reach2hru.dot(runoff[t,:]/1000.0/dt)/c_length[mapping] #m/s
 qin[qin < 0] = 0.0
 #Everyone hold here before continuing
 comm.Barrier()
 tic = time.time()
 #Assemble data to send
 recv = {}
 #Send data (Send to everyone as a first pass)
 for ucid in scids:
  db = {'qin':qin[odb[ucid]['mapping_ucid'][cid]['ocid']],
        'A0':A0[odb[ucid]['mapping_ucid'][cid]['ocid']]}
  comm.send(db,dest=list(cids).index(ucid),tag=11)
 stime += time.time() - tic
 tic = time.time()
 #Receive data
 for ucid in rcids:
  if ucid not in recv:recv[ucid] = {}
  db = comm.recv(source=list(cids).index(ucid),tag=11)
  for var in db:
   recv[ucid][var] = db[var]
 '''for var in ['qin','A0']:
  for ucid in scids:
   #if ucid == cid:continue
   if var == 'qin':data = qin[odb[ucid]['mapping_ucid'][cid]['ocid']]
   if var == 'A0':data = A0[odb[ucid]['mapping_ucid'][cid]['ocid']]
   comm.send(data,dest=list(cids).index(ucid),tag=11)
   #comm.Send(data,dest=list(cids).index(ucid),tag=13)
  #Receive data
  for ucid in rcids:
   if ucid not in recv:recv[ucid] = {}
   recv[ucid][var] = comm.recv(source=list(cids).index(ucid),tag=11)
   #recv[ucid][var] = np.empty(odb[cid]['mapping_ucid'][ucid]['ocid'].size,np.float64)
   #comm.Recv(recv[ucid][var], source=list(cids).index(ucid), tag=13)'''
 rtime += time.time() - tic
 #Update initial conditions using upstream information
 if t > 0:
  for ucid in ucids:
   if ucid == cid:continue
   A0_org[odb[cid]['mapping_ucid'][ucid]['cid']] = recv[ucid]['A0'][:]#[mapping_ucid[ucid]['ocid']]
 #Update the lateral inputs/outputs
 for ucid in ucids:
  if ucid == cid:continue
  qin[odb[cid]['mapping_ucid'][ucid]['cid']] = recv[ucid]['qin'][:]#[mapping_ucid[ucid]['ocid']]

 for it in range(max_niter):
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
  A1 = scipy.sparse.linalg.spsolve(LHS.tocsr(),RHS,use_umfpack=True)
  #A1 = scipy.sparse.linalg.spsolve(LHS,RHS,use_umfpack=True)
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
print(cid,'stime',stime/nt)
print(cid,'rtime',rtime/nt)
#print(np.sum(dVh),np.sum(dVQ))
#print(np.sum(dVQ)/np.sum(area))
pickle.dump(out,open('../workspace/%s_parallel.pck' % cid,'wb'))
