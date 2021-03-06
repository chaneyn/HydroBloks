import os
import numpy as np
import netCDF4 as nc
import h5py
import datetime
import time
import sys
import scipy.sparse as sparse
import pickle

#def assign_string(nelem,pstring):
def assign_string(dtype,pstring):

 #print(str(dtype))
 nelem = int(str(dtype)[2:])
 #nelem = len)
 tmp = np.chararray(shape=(nelem,))
 tmp[:] = ' '
 tmp2 = []
 for letter in pstring:
  #print(letter)
  tmp2.append(letter)
 tmp[0:len(tmp2)] = tmp2

 return tmp

def initialize(info):

 HB = HydroBlocks(info)

 return HB

class HydroBlocks:

 def __init__(self,info):

  #Read in general information
  self.general_information(info)

  #Initialize Noah-MP
  print("Initializing Noah-MP")
  self.initialize_noahmp()

  #Initialize subsurface module
  print("Initializing subsurface module")
  self.initialize_subsurface()

  #Initialize human water use module
  #print("Initializing Human Water Management")
  self.initialize_hwu(info)

  #Other metrics
  self.dE = 0.0
  self.r = 0.0
  self.dr = 0.0
  self.et = 0.0
  self.etran = 0.0
  self.esoil = 0.0
  self.ecan = 0.0
  self.prcp = 0.0
  self.q = 0.0
  self.errwat = 0.0
  self.erreng = 0.0
  self.dzwt0 = np.zeros(self.nhru,dtype=np.float32)
  self.beg_wb = np.zeros(self.nhru,dtype=np.float32)
  self.end_wb = np.zeros(self.nhru,dtype=np.float32)
  self.itime = 0

  #Restart from initial conditions?
  self.restart()

  return
  
 def restart(self,):

  file_restart = '%s/%s.h5' % (self.metadata['restart']['dir'],self.idate.strftime('%Y-%m-%d'))
  if (os.path.exists(file_restart) == False 
      or self.metadata['restart']['flag'] == False):
    print("Cold startup")
    return

  #Read in the restart information
  fp = h5py.File(file_restart,'r')
  self.noahmp.smceq[:] = fp['smceq'][:]
  self.noahmp.albold[:] = fp['albold'][:]
  self.noahmp.sneqvo[:] = fp['sneqvo'][:]
  self.noahmp.stc[:] = fp['stc'][:]
  self.noahmp.sh2o[:] = fp['sh2o'][:]
  self.noahmp.smc[:] = fp['smc'][:]
  self.noahmp.tah[:] = fp['tah'][:]
  self.noahmp.eah[:] = fp['eah'][:]
  self.noahmp.fwet[:] = fp['fwet'][:]
  self.noahmp.canliq[:] = fp['canliq'][:]
  self.noahmp.canice[:] = fp['canice'][:]
  self.noahmp.tv[:] = fp['tv'][:]
  self.noahmp.tg[:] = fp['tg'][:]
  self.noahmp.qsfc1d[:] = fp['qsfc1d'][:]
  self.noahmp.qsnow[:] = fp['qsnow'][:]
  self.noahmp.isnow[:] = fp['isnow'][:]
  self.noahmp.zsnso[:] = fp['zsnso'][:]
  self.noahmp.sndpth[:] = fp['sndpth'][:]
  self.noahmp.swe[:] = fp['swe'][:]
  self.noahmp.snice[:] = fp['snice'][:]
  self.noahmp.snliq[:] = fp['snliq'][:]
  self.noahmp.zwt[:] = fp['zwt'][:]
  self.noahmp.wa[:] = fp['wa'][:]
  self.noahmp.wt[:] = fp['wt'][:]
  self.noahmp.wslake[:] = fp['wslake'][:]
  self.noahmp.lfmass[:] = fp['lfmass'][:]
  self.noahmp.rtmass[:] = fp['rtmass'][:]
  self.noahmp.stmass[:] = fp['stmass'][:]
  self.noahmp.wood[:] = fp['wood'][:]
  self.noahmp.stblcp[:] = fp['stblcp'][:]
  self.noahmp.fastcp[:] = fp['fastcp'][:]
  self.noahmp.plai[:] = fp['plai'][:]
  self.noahmp.psai[:] = fp['psai'][:]
  self.noahmp.cm[:] = fp['cm'][:]
  self.noahmp.ch[:] = fp['ch'][:]
  self.noahmp.tauss[:] = fp['tauss'][:]
  self.noahmp.smcwtd[:] = fp['smcwtd'][:]
  self.noahmp.deeprech[:] = fp['deeprech'][:]
  self.noahmp.rech[:] = fp['rech'][:]
  fp.close()

  return
 
 def general_information(self,info):

  #Define the metadata
  self.dt = info['dt']
  self.dtt = self.dt#info['dtt']
  self.nsoil = len(info['dz'])#['nsoil']
  self.ncores = info['ncores']
  self.idate = info['idate']
  self.fdate = info['fdate']
  self.mkl_flag = info['mkl_flag']
  self.dt_timedelta = datetime.timedelta(seconds=self.dt)
  self.input_fp = nc.Dataset(info['input_file'])
  self.dx = self.input_fp.groups['metadata'].dx
  self.nhru = len(self.input_fp.dimensions['hsu'])
  self.surface_flow_flag = info['surface_flow_flag']
  self.subsurface_module = info['subsurface_module']
  self.hwu_flag = info['water_management']['hwu_flag']
  self.pct = self.input_fp.groups['parameters'].variables['area_pct'][:]/100
  self.pct = self.pct/np.sum(self.pct)
  self.metadata = info
  self.m =  self.input_fp.groups['parameters'].variables['m'][:]  #Noemi

  #Create a list of all the dates
  dates = []
  date = self.idate
  while date < self.fdate:
   dates.append(date)
   date = date + self.dt_timedelta
  self.dates = np.array(dates)

  return

 def initialize_noahmp(self,):

  #Initialize noahmp
  from pyNoahMP.NoahMP import model as noahmp
  self.noahmp = noahmp

  #Initialize parameters
  self.noahmp.ncells = self.nhru
  self.noahmp.nsoil = self.nsoil
  self.noahmp.dt = self.dt
  self.noahmp.nsnow = 3
  #self.noahmp.llanduse[:] = 'MODIFIED_IGBP_MODIS_NOAH'
  self.noahmp.llanduse = assign_string(self.noahmp.llanduse.dtype,'MODIFIED_IGBP_MODIS_NOAH')
  #self.noahmp.lsoil[:] = 'CUST'
  self.noahmp.lsoil = assign_string(self.noahmp.lsoil.dtype,'CUST')
  dir = os.path.dirname(os.path.abspath(__file__))
  VEGPARM = '%s/pyNoahMP/data/VEGPARM.TBL' % dir
  GENPARM = '%s/pyNoahMP/data/GENPARM.TBL' % dir
  MPTABLE = '%s/pyNoahMP/data/MPTABLE.TBL' % dir
  #self.noahmp.vegparm_file[0:len(VEGPARM)] = VEGPARM
  #self.noahmp.genparm_file[0:len(GENPARM)] = GENPARM
  #self.noahmp.mptable_file[0:len(MPTABLE)] = MPTABLE
  self.noahmp.vegparm_file = assign_string(self.noahmp.vegparm_file.dtype,VEGPARM)
  self.noahmp.genparm_file = assign_string(self.noahmp.genparm_file.dtype,GENPARM)
  self.noahmp.mptable_file = assign_string(self.noahmp.mptable_file.dtype,MPTABLE)
  #Define the options
  self.noahmp.idveg = 1#4 # dynamic vegetation (1 -> off ; 2 -> on)
  self.noahmp.iopt_crs = 1#2 # canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
  self.noahmp.iopt_btr = 1#2 # soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
  self.noahmp.iopt_run = 5#2 # runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
  self.noahmp.iopt_sfc = 1#2 # surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
  self.noahmp.iopt_frz = 2#1#1 # supercooled liquid water (1-> NY06; 2->Koren99)
  self.noahmp.iopt_inf = 2#2 # frozen soil permeability (1-> NY06; 2->Koren99)
  self.noahmp.iopt_rad = 1 # radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
  self.noahmp.iopt_alb = 2 # snow surface albedo (1->BATS; 2->CLASS)
  self.noahmp.iopt_snf = 3 # rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)]
  self.noahmp.iopt_tbot = 2#1#1 # lower boundary of soil temperature (1->zero-flux; 2->Noah) 
  self.noahmp.iopt_stc = 2#2 # snow/soil temperature time scheme (only layer 1) 1 -> semi-implicit; 2 -> full implicit (original Noah)

  #Allocate memory
  self.noahmp.initialize_general()

  #Set info
  self.noahmp.iz0tlnd = 0
  self.noahmp.sldpth[:] = np.array(self.metadata['dz'])
  self.noahmp.z_ml[:] = np.sum(self.noahmp.sldpth) 
  self.noahmp.zsoil[:] = -np.cumsum(self.noahmp.sldpth[:],axis=1)
  self.noahmp.zsnso[:] = 0.0
  self.noahmp.zsnso[:,3::] = self.noahmp.zsoil[:]
  
  #Set all to land (can include lake in the future...)
  self.noahmp.ist[:] = 1
  #Set soil color
  self.noahmp.isc[:] = 4#4
  self.noahmp.ice[:] = 0
  #Initialize snow info
  self.noahmp.isnow[:] = 0
  #Others
  self.noahmp.dx = self.dx
  self.noahmp.ice[:] = 0
  self.noahmp.foln[:] = 1.0
  self.noahmp.ficeold[:] = 0.0
  self.noahmp.albold[:] = 0.65
  self.noahmp.sneqvo[:] = 0.0
  self.noahmp.ch[:] = 0.0001
  self.noahmp.cm[:] = 0.0001
  self.noahmp.canliq[:] =  0.0
  self.noahmp.canice[:] = 0.0
  self.noahmp.sndpth[:] = 0.0
  self.noahmp.swe[:] = 0.0
  self.noahmp.snice[:] = 0.0
  self.noahmp.snliq[:] = 0.0
  self.noahmp.wa[:] = 0.0
  self.noahmp.wt[:] = self.noahmp.wa[:]
  self.noahmp.zwt[:] = 0.0
  self.noahmp.wslake[:] = 0.0
  self.noahmp.lfmass[:] = 9.0
  self.noahmp.rtmass[:] = 500.0
  self.noahmp.stmass[:] = 3.33
  self.noahmp.wood[:] = 500.0
  self.noahmp.stblcp[:] = 1000.0
  self.noahmp.fastcp[:] = 1000.0
  self.noahmp.plai[:] = 0.5
  self.noahmp.tauss[:] = 0.0
  self.noahmp.smcwtd[:] = self.noahmp.sh2o[:,1]
  self.noahmp.deeprech[:] = 0.0
  self.noahmp.rech[:] = 0.0
  self.noahmp.eah[:] = 1000.0
  self.noahmp.fwet[:] = 0.0
  self.noahmp.psai[:] = 0.1
  self.noahmp.stc[:] = 285.0
  self.noahmp.slopetyp[:] = 3
  self.noahmp.albold[:] = 0.5
  #Define the data
  self.noahmp.vegtyp[:] = self.input_fp.groups['parameters'].variables['land_cover'][:]
  self.noahmp.soiltyp[:] = np.arange(1,self.noahmp.ncells+1)
  self.noahmp.clay_pct  = self.input_fp.groups['parameters'].variables['clay'][:] # Noemi
  self.noahmp.smcmax[:] = self.input_fp.groups['parameters'].variables['MAXSMC'][:]
  self.noahmp.smcref[:] = self.input_fp.groups['parameters'].variables['REFSMC'][:]
  self.noahmp.smcdry[:] = self.input_fp.groups['parameters'].variables['DRYSMC'][:]
  for ilayer in range(self.noahmp.sh2o.shape[1]):
   self.noahmp.sh2o[:,ilayer] = self.input_fp.groups['parameters'].variables['MAXSMC'][:]
   #self.noahmp.sh2o[:,ilayer] = self.input_fp.groups['parameters'].variables['REFSMC'][:] #Noemi, start at REFSMC
  self.noahmp.smc[:] = self.noahmp.sh2o[:]
  self.noahmp.smcwtd[:] = self.noahmp.sh2o[:,0]
  #Initialize the soil parameters
  self.noahmp.bb0[:] = self.input_fp.groups['parameters'].variables['BB'][:]
  self.noahmp.drysmc0[:] = self.input_fp.groups['parameters'].variables['DRYSMC'][:]
  self.noahmp.f110[:] = self.input_fp.groups['parameters'].variables['F11'][:]
  self.noahmp.maxsmc0[:] = self.input_fp.groups['parameters'].variables['MAXSMC'][:]
  self.noahmp.refsmc0[:] = self.input_fp.groups['parameters'].variables['REFSMC'][:]
  self.noahmp.satpsi0[:] = self.input_fp.groups['parameters'].variables['SATPSI'][:]
  self.noahmp.satdk0[:] = self.input_fp.groups['parameters'].variables['SATDK'][:]
  self.noahmp.satdw0[:] = self.input_fp.groups['parameters'].variables['SATDW'][:]
  self.noahmp.wltsmc0[:] = self.input_fp.groups['parameters'].variables['WLTSMC'][:]
  self.noahmp.qtz0[:] = self.input_fp.groups['parameters'].variables['QTZ'][:]

  #Set lat/lon (declination calculation)
  self.noahmp.lat[:] = 0.0174532925*self.input_fp.groups['metadata'].latitude
  self.noahmp.lon[:] = 0.0174532925*self.input_fp.groups['metadata'].longitude

  #Initialize output
  self.noahmp.tg[:] = 285.0
  self.noahmp.tv[:] = 285.0
  for ilayer in range(self.noahmp.nsnow-1,self.noahmp.stc.shape[1]):
   self.noahmp.stc[:,ilayer] = 285.0
  self.noahmp.tah[:] = 285.0
  self.noahmp.t2mv[:] = 285.0
  self.noahmp.t2mb[:] = 285.0
  self.noahmp.runsf[:] = 0.0
  self.noahmp.runsb[:] = 0.0
  #forcing
  self.noahmp.fveg[:] = 1.0
  self.noahmp.fvgmax[:] = 1.0
  self.noahmp.tbot[:] = 285.0

  #Define the parameters
  self.noahmp.initialize_parameters()

  return

 def initialize_subsurface(self,):

  if self.subsurface_module == 'dtopmodel':self.initialize_dtopmodel()
  elif self.subsurface_module == 'richards':self.initialize_richards()

  return

 def initialize_richards(self,):
   
  from pyRichards import richards
  
  #Initialize richards
  self.richards = richards.richards(self.nhru,self.nsoil)

  #Set other parameters
  self.richards.dx = self.dx
  self.richards.nhru = self.nhru
  #print(self.nhru)
  self.richards.m[:] = self.input_fp.groups['parameters'].variables['m'][:]
  #self.richards.dem[:] = self.input_fp.groups['parameters'].variables['dem'][:]
  self.richards.dem[:] = self.input_fp.groups['parameters'].variables['hand'][:]
  self.richards.slope[:] = self.input_fp.groups['parameters'].variables['slope'][:]
  #self.richards.hand[:] = self.input_fp.groups['parameters'].variables['hand'][:]
  self.richards.area[:] = self.input_fp.groups['parameters'].variables['area'][:]
  self.richards.width = sparse.csr_matrix((self.input_fp.groups['wmatrix'].variables['data'][:],
                                  self.input_fp.groups['wmatrix'].variables['indices'][:],
                                  self.input_fp.groups['wmatrix'].variables['indptr'][:]),
                                  shape=(self.nhru,self.nhru),dtype=np.float64)
                                  #dtype=np.float64)[0:self.nhru,0:self.nhru]
  #print(self.richards.width.shape)
  self.richards.I = self.richards.width.copy()
  self.richards.I[self.richards.I != 0] = 1
  
  return

 def initialize_dtopmodel(self,):

  from pyDTopmodel import dynamic_topmodel as dtopmodel

  #Define some metadata
  nhru_outlet = self.input_fp.groups['outlet'].groups['summary'].variables['hru_dst'].size

  #Initialize Dynamic Topmodel
  self.dtopmodel = dtopmodel.Dynamic_Topmodel(self.nhru,nhru_outlet,self.mkl_flag)

  #Set flags
  self.dtopmodel.surface_flow_flag = self.surface_flow_flag
  #self.dtopmodel.subsurface_flow_flag = self.subsurface_flow_flag

  #Set self.dtopmodel parameters
  self.dtopmodel.dt = self.dt #seconds
  self.dtopmodel.area[:] = 0.0 #meters^2
  self.dtopmodel.dx[:] = self.dx #meters
  #print self.input_fp.groups['parameters'].variables['m'][:]
  self.dtopmodel.m[:] = self.input_fp.groups['parameters'].variables['m'][:] #Noemi
  #self.dtopmodel.sdmax[:] = 100.0#self.input_fp.groups['parameters'].variables['sdmax'][:]
  self.dtopmodel.sdmax[:] = self.input_fp.groups['parameters'].variables['sdmax'][:] #Noemi
  #Set cluster information
  self.dtopmodel.pct[:] = self.input_fp.groups['parameters'].variables['area_pct'][:]/100
  self.dtopmodel.area[:] = self.input_fp.groups['parameters'].variables['area'][:]
  af = 1.0 #anisotropy factor
  self.dtopmodel.T0[:] = af*self.input_fp.groups['parameters'].variables['SATDK'][:]*self.dtopmodel.m
  self.dtopmodel.sti[:] = self.input_fp.groups['parameters'].variables['ti'][:]
  self.dtopmodel.beta[:] = self.input_fp.groups['parameters'].variables['slope'][:]
  self.dtopmodel.carea[:] = self.input_fp.groups['parameters'].variables['carea'][:]
  self.dtopmodel.channel[:] = self.input_fp.groups['parameters'].variables['channel'][:]
  self.dtopmodel.dem[:] = self.input_fp.groups['parameters'].variables['dem'][:] 
  self.dtopmodel.mannings[:] = self.input_fp.groups['parameters'].variables['mannings'][:]
  #Set outlet information
  self.dtopmodel.area_outlet[:] = self.dx**2*self.input_fp.groups['outlet'].groups['summary'].variables['counts'][:]
  self.dtopmodel.pct = self.dtopmodel.pct/np.sum(self.dtopmodel.pct)
  ti_mean = np.sum(self.dtopmodel.pct*self.dtopmodel.sti[:])
 
  #Calculate the sti
  lnT0 = np.log(self.dtopmodel.T0)
  lnTe = np.sum(self.dtopmodel.pct*lnT0)
  self.dtopmodel.sti = self.dtopmodel.sti - (lnT0 - lnTe)

  #Set weight matrix
  self.dtopmodel.flow_matrix = sparse.csr_matrix((self.input_fp.groups['flow_matrix'].variables['data'][:],
			          self.input_fp.groups['flow_matrix'].variables['indices'][:],
			          self.input_fp.groups['flow_matrix'].variables['indptr'][:]),
 				  dtype=np.float64)[0:self.nhru,0:self.nhru]
  self.dtopmodel.flow_matrix.setdiag(self.dtopmodel.flow_matrix.diagonal()) #Ensure the zeros are not sparse (for kinematic wave solution).
  self.dtopmodel.flow_matrix_T = sparse.csr_matrix(self.dtopmodel.flow_matrix.T) #transposed
  self.dtopmodel.flow_matrix_T.setdiag(self.dtopmodel.flow_matrix_T.diagonal()) #Ensure the zeros are not sparse  (for kinematic wave solution).

  #Initialize the solver
  if self.mkl_flag: self.dtopmodel.dtt.initialize(self.dtopmodel.flow_matrix_T.indices,self.dtopmodel.flow_matrix_T.indptr)

  #Initialize the soil moisture deficit values
  self.dtopmodel.si[:] = 0.0
 
  return

 def initialize_hwu(self,info):

  from pyHWU.Human_Water_Use import Human_Water_Use as hwu
  self.hwu = hwu(self,info)
  self.hwu.area = self.input_fp.groups['parameters'].variables['area'][:]

  if self.hwu.hwu_flag == True:
   print("Initializing Human Water Management")
   self.hwu.initialize_allocation(self,)

  return

 def run(self,info):

  #Run the model
  date = self.idate
  tic = time.time()
  self.noahmp.dzwt[:] = 0.0

  if self.subsurface_module == 'dtopmodel': self.dtopmodel.ex[:] = 0.0
  while date < self.fdate:

   #Update input data
   self.update_input(date)

   #Save the original precip
   precip = np.copy(self.noahmp.prcp)

   #Calculate initial NOAH water balance
   self.initialize_water_balance()

   #Update model
   self.update(date)
   
   #Return precip to original value
   self.noahmp.prcp[:] = precip[:] 

   #Calculate final water balance
   self.finalize_water_balance()

   #Update the water balance error
   self.calculate_water_balance_error()

   #Update the energy balance error
   self.calculate_energy_balance_error()
 
   #Update time and date
   self.date = date
   info['date'] = date

   #Update output
   self.update_output(date)

   #Update time step
   date = date + self.dt_timedelta
   self.itime = self.itime + 1

   #Output some statistics
   if (date.hour == 0) and (date.day == 1):
    print(date.strftime("%Y-%m-%d"),'%10.4f'%(time.time()-tic),'et:%10.4f'%self.et,'prcp:%10.4f'%self.prcp,'q:%10.4f'%self.q,'WB ERR:%10.6f' % self.errwat,'ENG ERR:%10.6f' % self.erreng)

  return

 def update_input(self,date):

  self.noahmp.itime = self.itime
  dt = self.dt
  if self.subsurface_module == 'dtopmodel':self.dtopmodel.itime = self.itime
  self.noahmp.nowdate = assign_string(self.noahmp.nowdate.dtype,date.strftime('%Y-%m-%d_%H:%M:%S'))
  #self.noahmp.nowdate[:] = assign_string(len(self.noahmp.nowdate),date.strftime('%Y-%m-%d_%H:%M:%S'))

  #print(self.noahmp.nowdate[:])
  #print(self.noahmp.nowdate.dtype)
  #print(date.strftime('%Y-%m-%d_%H:%M:%S').dtype)
  #tmp = np.char.array(date.strftime('%Y-%m-%d_%H:%M:%S'))
  #exit()
  #self.noahmp.nowdate[:] = tmp[:]#date.strftime('%Y-%m-%d_%H:%M:%S')
  #print(len(self.noahmp.nowdate))

  self.noahmp.julian = (date - datetime.datetime(date.year,1,1,0)).days
  self.noahmp.yearlen = (datetime.datetime(date.year+1,1,1,0) - datetime.datetime(date.year,1,1,1,0)).days + 1

  #Update meteorology
  meteorology = self.input_fp.groups['meteorology']
  if date == self.idate:
   #determine the first time step for the meteorology
   var = meteorology.variables['time']
   ndates = var[:]
   #convert current date to num
   ndate = nc.date2num(date,units=var.units,calendar=var.calendar)
   #identify the position
   self.minitial_itime = np.where(ndates == ndate)[0][0]
  i = self.itime + self.minitial_itime
  self.noahmp.lwdn[:] = meteorology.variables['lwdown'][i,:] #W/m2
  self.noahmp.swdn[:] = meteorology.variables['swdown'][i,:] #W/m2
  self.noahmp.psfc[:] = meteorology.variables['psurf'][i,:] #Pa
  self.noahmp.p_ml[:] = meteorology.variables['psurf'][i,:] #Pa
  self.noahmp.u_ml[:] = (meteorology.variables['wind'][i,:]**2/2)**0.5 #m/s
  self.noahmp.v_ml[:] = (meteorology.variables['wind'][i,:]**2/2)**0.5 #m/s
  self.noahmp.t_ml[:] = meteorology.variables['tair'][i,:] #K
  self.noahmp.q_ml[:] = meteorology.variables['spfh'][i,:] #Kg/Kg
  self.noahmp.qsfc1d[:] = meteorology.variables['spfh'][i,:] #Kg/Kg
  self.noahmp.prcp[:] = meteorology.variables['precip'][i,:] #mm/s

  #Set the partial pressure of CO2 and O2
  self.noahmp.co2air[:] = 355.E-6*self.noahmp.psfc[:]# ! Partial pressure of CO2 (Pa) ! From NOAH-MP-WRF
  self.noahmp.o2air[:] = 0.209*self.noahmp.psfc[:]# ! Partial pressure of O2 (Pa)  ! From NOAH-MP-WRF

  # Update water demands
  if self.hwu.hwu_flag == True:
   if (date.hour*3600)%self.hwu.dta == 0:
    water_use = self.input_fp.groups['water_use']
    if self.hwu.hwu_indust_flag == True:
     self.hwu.demand_indust[:]  = water_use.variables['industrial'][i,:] #m/s
     self.hwu.deficit_indust[:] = np.copy(self.hwu.demand_indust[:])
    if self.hwu.hwu_domest_flag == True:
     self.hwu.demand_domest[:]  = water_use.variables['domestic'][i,:] #m/s
     self.hwu.deficit_domest[:] = np.copy(self.hwu.demand_domest[:])
    if self.hwu.hwu_lstock_flag == True:
     self.hwu.demand_lstock[:]  = water_use.variables['livestock'][i,:] #m/s
     self.hwu.deficit_lstock[:] = np.copy(self.hwu.demand_lstock[:])


  return

 def update(self,date):

  # Apply irrigation
  self.hwu.Human_Water_Irrigation(self,date)
  
  # Update subsurface
  self.update_subsurface()

  # Update NOAH
  self.noahmp.run_model(self.ncores)

  # Calculate water demands and supplies, and allocate volumes
  self.hwu.Calc_Human_Water_Demand_Supply(self,date)

  # Abstract Surface Water and Groundwater
  self.hwu.Water_Supply_Abstraction(self,date)



  return

 def update_subsurface(self,):

  self.noahmp.dzwt[:] = 0.0

  if self.subsurface_module == 'dtopmodel':

   #Calculate the updated soil moisture deficit
   si0 = np.copy(self.noahmp.si0)
   si1 = np.copy(self.noahmp.si1)

   #Calculate the change in deficit
   self.dtopmodel.si[:] = si1[:]
   self.dtopmodel.dsi[:] = np.copy(si1 - si0)
   self.dtopmodel.r[:] = -self.dtopmodel.dsi[:]/self.dtopmodel.dt

   #Add the surface runoff
   self.dtopmodel.qsurf[:] = self.noahmp.runsf[:]/1000.0

   #Update dynamic topmodel
   self.dtopmodel.update(self.ncores)

   #Calculate the change in deficit
   #self.dtopmodel.sideep[:] = self.dtopmodel.sideep.astype(np.float32)
   #self.dtopmodel.si[:] = self.dtopmodel.si.astype(np.float32)
   #dsi = np.copy(si1 - self.dtopmodel.si)

   #Update the soil moisture values
   #print self.dtopmodel.dt*self.dtopmodel.r
   #self.noahmp.dzwt[:] = dsi+self.dtopmodel.dt*self.dtopmodel.ex-self.dtopmodel.dt*self.dtopmodel.r
   #self.noahmp.dzwt[:] = 0
   self.noahmp.hdiv[:] = 0
   ls = []
   for ihru in range(self.nhru):
    cs = np.cumsum(self.noahmp.sldpth[ihru,:])
    m = cs>np.abs(self.noahmp.zwt[ihru])
    if np.sum(m) > 0:
     idx = np.where(m)[0][0]
     #weight the extraction
     fs = self.noahmp.sldpth[ihru,idx:]/np.sum(self.noahmp.sldpth[ihru,idx:])
     self.noahmp.hdiv[ihru,idx:] = 1000*fs*(self.dtopmodel.qout[ihru]-self.dtopmodel.qin[ihru])
    else:
     self.noahmp.hdiv[ihru,:] = 0.0

  elif self.subsurface_module == 'richards':

   #Assign noahmp variables to subsurface module
   self.richards.theta[:] = self.noahmp.smc[:]
   self.richards.thetar[:] = self.noahmp.drysmc0[:]
   self.richards.thetas[:] = self.noahmp.maxsmc0[:]
   self.richards.b[:] = self.noahmp.bb0[:]
   self.richards.satpsi[:] = self.noahmp.satpsi0[:]
   self.richards.ksat[:] = self.noahmp.satdk0[:]
   self.richards.dz[:] = self.noahmp.sldpth[:]

   #Update subsurface module
   self.richards.update()

   #Assign subsurface module variables to noahmp
   self.noahmp.hdiv[:] = self.richards.hdiv[:]
   
  return

 def initialize_water_balance(self,):
 
  smw = np.sum(1000*self.noahmp.sldpth*self.noahmp.smc,axis=1)
  self.beg_wb = np.copy(self.noahmp.canliq + self.noahmp.canice + self.noahmp.swe + self.noahmp.wa + smw)
  self.dzwt0 = np.copy(self.noahmp.dzwt)

  return 

 def finalize_water_balance(self,):

  NOAH = self.noahmp
  smw = np.sum(1000*NOAH.sldpth*NOAH.smc,axis=1) 
  self.end_wb = np.copy(NOAH.canliq + NOAH.canice + NOAH.swe + NOAH.wa + smw)
  
  return

 def calculate_water_balance_error(self,):

  NOAH = self.noahmp
  HWU = self.hwu
  dt = self.dt
  if self.subsurface_module == 'dtopmodel':
   #tmp = np.copy(self.end_wb - self.beg_wb - NOAH.dt*(NOAH.prcp-NOAH.ecan-
   #      NOAH.etran-NOAH.esoil-NOAH.runsf-NOAH.runsb) - 1000*self.dzwt0)
   tmp = np.copy(self.end_wb - self.beg_wb - NOAH.dt*(NOAH.prcp-NOAH.ecan-
         NOAH.etran-NOAH.esoil-NOAH.runsf-NOAH.runsb-np.sum(NOAH.hdiv,axis=1)))
  elif self.subsurface_module == 'richards':
   tmp = np.copy(self.end_wb - self.beg_wb - NOAH.dt*(NOAH.prcp-NOAH.ecan-
         NOAH.etran-NOAH.esoil-NOAH.runsf-NOAH.runsb-np.sum(NOAH.hdiv,axis=1)))
  else:
   tmp = np.copy(self.end_wb - self.beg_wb - NOAH.dt*(NOAH.prcp-NOAH.ecan-
         NOAH.etran-NOAH.esoil-NOAH.runsf-NOAH.runsb))
  self.errwat += np.sum(self.pct*tmp)
  self.q = self.q + dt*np.sum(self.pct*NOAH.runsb) + dt*np.sum(self.pct*NOAH.runsf)
  self.et = self.et + dt*np.sum(self.pct*(NOAH.ecan + NOAH.etran + NOAH.esoil))
  self.etran += dt*np.sum(self.pct*NOAH.etran)
  self.ecan += dt*np.sum(self.pct*NOAH.ecan)
  self.esoil += dt*np.sum(self.pct*NOAH.esoil)
  self.prcp = self.prcp + dt*np.sum(self.pct*NOAH.prcp)

  return

 def calculate_energy_balance_error(self,):

  NOAH = self.noahmp
  tmp = np.copy(NOAH.sav+NOAH.sag-NOAH.fira-NOAH.fsh-NOAH.fcev-NOAH.fgev-NOAH.fctr-NOAH.ssoil)
  self.erreng += np.sum(self.pct*tmp)

  return

 def update_output(self,date):

  NOAH = self.noahmp
  HWU = self.hwu
  HB = self
  itime = self.itime

  #Create the netcdf file
  if date == self.idate: self.create_netcdf_file()

  #General info
  grp = self.output_fp.groups['metadata']
  dates = grp.variables['date']
  dates[itime] = nc.date2num(date,units=dates.units,calendar=dates.calendar)
  
  #Update the variables
  grp = self.output_fp.groups['data']

  
  tmp = {}
  #NoahMP
  tmp['smc'] = np.copy(NOAH.smc) #m3/m3
  tmp['g'] = np.copy(NOAH.ssoil) #W/m2
  tmp['sh'] = np.copy(NOAH.fsh) #W/m2
  tmp['lh'] = np.copy(NOAH.fcev + NOAH.fgev + NOAH.fctr) #W/m2
  tmp['qbase'] = NOAH.dt*np.copy(NOAH.runsb) #mm
  tmp['qsurface'] = NOAH.dt*np.copy(NOAH.runsf) #mm
  tmp['runoff'] = NOAH.dt*(np.copy(NOAH.runsf)+np.copy(NOAH.runsb)) #mm 
  tmp['prcp'] = NOAH.dt*np.copy(NOAH.prcp) #mm
  tmp['trad'] = np.copy(NOAH.trad) #K
  tmp['stc'] = np.copy(NOAH.stc[:,3]) #K
  tmp['tv'] = np.copy(NOAH.tv) #K
  tmp['salb'] = np.copy(NOAH.salb) 
  tmp['wtd'] = np.copy(NOAH.zwt) #m
  #tmp['totsmc'] = np.sum(NOAH.sldpth*NOAH.smc,axis=1)/np.sum(NOAH.sldpth[0]) #m3/m3
  tmp['hdiv'] = np.copy(NOAH.hdiv)

  # root zone
  cs = np.cumsum(NOAH.sldpth[0,:])
  mask = cs <= 0.5
  pct = NOAH.sldpth[0,mask]/np.sum(NOAH.sldpth[0,mask])
  tmp['smc_root'] = np.sum(pct*NOAH.smc[:,mask],axis=1) #m3/m3
  # top soil layer
  tmp['smc1'] = np.copy(NOAH.smc[:,0])
  # total soil moisture / soil water storage -- only until depth to bedrock
  dl = [ np.argmin(np.abs(cs-self.m[i])) for i in range(self.nhru)]
  tmp['totsmc'] = [ np.sum(NOAH.sldpth[i,:dl[i]]*NOAH.smc[i,:dl[i]])/np.sum(NOAH.sldpth[i,:dl[i]])  for i in range(self.nhru)]
  
 
  #Dynamic TOPMODEL
  if self.subsurface_module == 'dtopmodel':
   TOPMODEL = self.dtopmodel
   tmp['swd'] = np.copy(10**3*TOPMODEL.si) #mm
   tmp['qout_subsurface'] = np.copy(TOPMODEL.qout) #m2/s
   tmp['qout_surface'] = np.copy(TOPMODEL.qout_surface) #m2/s
   tmp['sstorage'] = np.copy(TOPMODEL.storage_surface)

  #General
  tmp['errwat'] = np.copy(HB.errwat)

  # Water Management
  if self.hwu.hwu_agric_flag == True:
   tmp['demand_agric'] = np.copy(HWU.demand_agric)*NOAH.dt       #m
   tmp['deficit_agric'] = np.copy(HWU.deficit_agric)*NOAH.dt     #m
   tmp['irrig_agric'] = np.copy(HWU.irrigation)*(NOAH.dt/1000.0)   #m

  if self.hwu.hwu_flag == True:
   if self.hwu.hwu_indust_flag == True:
     tmp['demand_indust'] = np.copy(HWU.demand_indust)*NOAH.dt #m
     tmp['deficit_indust'] = np.copy(HWU.deficit_indust)*NOAH.dt #m
     tmp['alloc_indust'] = np.copy(HWU.demand_indust-HWU.deficit_indust)*NOAH.dt #m
   if self.hwu.hwu_domest_flag == True:
     tmp['demand_domest'] = np.copy(HWU.demand_domest)*NOAH.dt #m
     tmp['deficit_domest'] = np.copy(HWU.deficit_domest)*NOAH.dt #m
     tmp['alloc_domest'] = np.copy(HWU.demand_domest-HWU.deficit_domest)*NOAH.dt #m
   if self.hwu.hwu_lstock_flag == True:
     tmp['demand_lstock'] = np.copy(HWU.demand_lstock)*NOAH.dt #m
     tmp['deficit_lstock'] = np.copy(HWU.deficit_lstock)*NOAH.dt #m
     tmp['alloc_lstock'] = np.copy(HWU.demand_lstock-HWU.deficit_lstock)*NOAH.dt #m
   if self.hwu.hwu_sf_flag == True:
     tmp['alloc_sf'] = np.copy(HWU.alloc_sf) #m
   if self.hwu.hwu_gw_flag == True:                                                                                              tmp['alloc_gw'] = np.copy(HWU.alloc_gw) #m

  #Output the variables
  for var in self.metadata['output']['vars']:
   grp.variables[var][itime,:] = tmp[var][:]


  return

 def create_netcdf_file(self,):

  #Create the output directory if necessary
  os.system('mkdir -p %s' % self.metadata['output']['dir'])

  #Extract the pointers to both input and output files
  ofile = '%s/%s.nc' % (self.metadata['output']['dir'],self.idate.strftime('%Y-%m-%d'))
  self.output_fp = nc.Dataset(ofile,'w',format='NETCDF4')
  fp_out = self.output_fp
  fp_in = self.input_fp

  #Define the metadata
  metadata = {
             'g':{'description':'Ground heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             'sh':{'description':'Sensible heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             'lh':{'description':'Latent heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             #'lwnet':{'description':'Net longwave radiation','units':'W/m2','dims':('time','hru',),'precision':4},
             #'swnet':{'description':'Absorbed shortwave radiation','units':'W/m2','dims':('time','hru',),'precision':4},
             'trad':{'description':'Land surface skin temperature','units':'K','dims':('time','hru',),'precision':2},
             'stc':{'description': 'Snow/soil temperature','units':'K','dims':('time','hru',),'precision':2},
             'tv':{'description': 'Canopy temperature','units':'K','dims':('time','hru',),'precision':2},
             'salb':{'description':'Land surface albedo','units':' ','dims':('time','hru',),'precision':4},

             'qbase':{'description':'Excess runoff','units':'mm/s','dims':('time','hru',),'precision':4},
             'qsurface':{'description':'Surface runoff','units':'mm/s','dims':('time','hru',),'precision':4},
             'runoff':{'description':'Total runoff','units':'mm/s','dims':('time','hru',),'precision':4},
             'prcp':{'description':'Precipitation','units':'mm/s','dims':('time','hru',),'precision':4},
             'smc':{'description':'Soil water content','units':'m3/m3','dims':('time','hru','soil'),'precision':3},
             'swd':{'description':'Soil water deficit','units':'mm','dims':('time','hru',),'precision':4},
             'sstorage':{'description':'Surface storage','units':'mm','dims':('time','hru',),'precision':4},
             #'sndpth':{'description':'Snow depth','units':'mm','dims':('time','hru',)},
             #'swe':{'description':'Snow water equivalent','units':'mm','dims':('time','hru',)},

             'qout_subsurface':{'description':'Subsurface flux','units':'m2/s','dims':('time','hru',),'precision':4},
             'qout_surface':{'description':'Surface flux','units':'m2/s','dims':('time','hru',),'precision':4},
             'wtd':{'description':'WTD','units':'m','dims':('time','hru',),'precision':2},
             'errwat':{'description':'errwat','units':'mm','dims':('time','hru',),'precision':4},
             'totsmc':{'description':'totsmc','units':'m3/m3','dims':('time','hru',),'precision':2},
             'hdiv':{'description':'hdiv','units':'mm/s','dims':('time','hru','soil'),'precision':4},

             'smc1':{'description':'Soil water content at the root zone','units':'m3/m3','dims':('time','hru',),'precision':3},
             'smc_root':{'description':'Soil water content at the root zone','units':'m3/m3','dims':('time','hru',),'precision':3},  

             # Water Management
             'demand_agric':{'description':'Irrigation water demand','units':'m','dims':('time','hru',),'precision':4},
             'deficit_agric':{'description':'Irrigation water deficit','units':'m','dims':('time','hru',),'precision':4},
             'irrig_agric':{'description':'Irrigated water volume','units':'m','dims':('time','hru',),'precision':4},
             'demand_indust':{'description':'Industrial water demand','units':'m','dims':('time','hru',)},'precision':4,
             'deficit_indust':{'description':'Industrial water deficit','units':'m','dims':('time','hru',),'precision':4},
             'alloc_indust':{'description':'Industrial water allocated','units':'m','dims':('time','hru',),'precision':4},
             'demand_domest':{'description':'Domestic demand','units':'m','dims':('time','hru',),'precision':4},
             'deficit_domest':{'description':'Domestic deficit','units':'m','dims':('time','hru',),'precision':4},
             'alloc_domest':{'description':'Domestic water allocated','units':'m','dims':('time','hru',),'precision':4},
             'demand_lstock':{'description':'Livestock demand','units':'m','dims':('time','hru',),'precision':4},
             'deficit_lstock':{'description':'Livestock deficit','units':'m','dims':('time','hru',),'precision':4},
             'alloc_lstock':{'description':'Livestock water allocated','units':'m','dims':('time','hru',),'precision':4},
             'alloc_sf':{'description':'Surface water allocated','units':'m','dims':('time','hru',),'precision':4},
             'alloc_gw':{'description':'Groundwater water allocated','units':'m','dims':('time','hru',),'precision':4}
             }

  #Create the dimensions
  print('Creating the dimensions')
  ntime = 24*3600*((self.fdate - self.idate).days)/self.dt
  nhru = len(fp_in.dimensions['hsu'])
  fp_out.createDimension('hru',nhru)
  fp_out.createDimension('time',ntime)
  fp_out.createDimension('soil',self.nsoil)

  #Create the output
  print('Creating the data group')
  grp = fp_out.createGroup('data')
  for var in self.metadata['output']['vars']:
   ncvar = grp.createVariable(var,'f4',metadata[var]['dims'],least_significant_digit=metadata[var]['precision'])#,zlib=True)
   ncvar.description = metadata[var]['description']
   ncvar.units = metadata[var]['units']

  #Create the metadata
  print('Creating the metadata group')
  grp = fp_out.createGroup('metadata')
  #Time
  grp.createVariable('time','f8',('time',))
  dates = grp.createVariable('date','f8',('time',))
  dates.units = 'hours since 1900-01-01'
  dates.calendar = 'standard'

  #HRU percentage coverage
  print('Setting the HRU percentage coverage')
  pcts = grp.createVariable('pct','f4',('hru',))
  pcts[:] = fp_in.groups['parameters'].variables['area_pct'][:]
  pcts.description = 'hru percentage coverage'
  pcts.units = '%'

  #HRU area
  print('Setting the HRU areal coverage')
  area = grp.createVariable('area','f4',('hru',))
  area[:] = fp_in.groups['parameters'].variables['area'][:]
  area.units = 'meters squared'
  print('Defining the HRU ids')
  hru = grp.createVariable('hru','i4',('hru',))
  hrus =[]
  for value in range(nhru):hrus.append(value)
  hru[:] = np.array(hrus)
  hru.description = 'hru ids'

  return

 def finalize(self,):

  #Create the restart directory if necessary
  os.system('mkdir -p %s' % self.metadata['restart']['dir'])

  #Save the restart file
  file_restart = '%s/%s.h5' % (self.metadata['restart']['dir'],self.fdate.strftime('%Y-%m-%d'))
  fp = h5py.File(file_restart,'w')
  fp['smceq'] = self.noahmp.smceq[:]
  fp['albold'] = self.noahmp.albold[:]
  fp['sneqvo'] = self.noahmp.sneqvo[:]
  fp['stc'] = self.noahmp.stc[:]
  fp['sh2o'] = self.noahmp.sh2o[:]
  fp['smc'] = self.noahmp.smc[:]
  fp['tah'] = self.noahmp.tah[:]
  fp['eah'] = self.noahmp.eah[:]
  fp['fwet'] = self.noahmp.fwet[:]
  fp['canliq'] = self.noahmp.canliq[:]
  fp['canice'] = self.noahmp.canice[:]
  fp['tv'] = self.noahmp.tv[:]
  fp['tg'] = self.noahmp.tg[:]
  fp['qsfc1d'] = self.noahmp.qsfc1d[:]
  fp['qsnow'] = self.noahmp.qsnow[:]
  fp['isnow'] = self.noahmp.isnow[:]
  fp['zsnso'] = self.noahmp.zsnso[:]
  fp['sndpth'] = self.noahmp.sndpth[:]
  fp['swe'] = self.noahmp.swe[:]
  fp['snice'] = self.noahmp.snice[:]
  fp['snliq'] = self.noahmp.snliq[:]
  fp['zwt'] = self.noahmp.zwt[:]
  fp['wa'] = self.noahmp.wa[:]
  fp['wt'] = self.noahmp.wt[:]
  fp['wslake'] = self.noahmp.wslake[:]
  fp['lfmass'] = self.noahmp.lfmass[:]
  fp['rtmass'] = self.noahmp.rtmass[:]
  fp['stmass'] = self.noahmp.stmass[:]
  fp['wood'] = self.noahmp.wood[:]
  fp['stblcp'] = self.noahmp.stblcp[:]
  fp['fastcp'] = self.noahmp.fastcp[:]
  fp['plai'] = self.noahmp.plai[:]
  fp['psai'] = self.noahmp.psai[:]
  fp['cm'] = self.noahmp.cm[:]
  fp['ch'] = self.noahmp.ch[:]
  fp['tauss'] = self.noahmp.tauss[:]
  fp['smcwtd'] = self.noahmp.smcwtd[:]
  fp['deeprech'] = self.noahmp.deeprech[:]
  fp['rech'] = self.noahmp.rech[:]
  fp.close()
   
  #Close the LSM
  self.noahmp.finalize()
  del self.noahmp

  #Close the subsurface module
  if self.subsurface_module == 'dtopmodel':
   if self.mkl_flag: self.dtopmodel.dtt.finalize()
   del self.dtopmodel

  #Close the files
  self.input_fp.close()
  self.output_fp.close()

  return
