import os
import numpy as np
import netCDF4 as nc
import h5py
import datetime
import time
import sys
import scipy.sparse as sparse
import pickle

def assign_string(nelem,pstring):
 #def assign_string(dtype,pstring):

 #print(str(dtype))
 #nelem = int(str(dtype)[2:])
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
  print("Initializing Noah-MP",flush=True)
  self.initialize_noahmp()
  self.tsno=np.zeros((self.noahmp.ncells,self.noahmp.nsnow),order='F').astype(np.float32) #Initial condition TSNO Laura
  #Initialize subsurface module
  print("Initializing subsurface module",flush=True)
  self.initialize_subsurface()

  #Initialize human water use module
  #print("Initializing Human Water Management",flush=True)
  #self.initialize_hwu(info)

  #Initialize routing module
  print("Initializing the routing module",flush=True)
  self.initialize_routing()
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
    print("Cold startup",flush=True)
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
  #routing
  if self.routing_module == 'kinematic':
   self.routing.Q0[:] = fp['Q0'][:]
   self.routing.u0[:] = fp['u0'][:]
   self.routing.A0[:] = fp['A0'][:]
   self.routing.Q1[:] = fp['Q1'][:]
   self.routing.A1[:] = fp['A1'][:]
   self.routing.dA[:] = fp['dA'][:]
   self.routing.bcs[:] = fp['bcs'][:]
   self.routing.qss[:] = fp['qss'][:]
   self.routing.hband_inundation[:] = fp['hband_inundation'][:]
   self.routing.reach2hband_inundation[:] = fp['reach2hband_inundation'][:]
  fp.close()

  return
 
 def general_information(self,info):

  #Parallel information
  self.MPI  = info['MPI']

  #Define the metadata
  self.cid = info['cid']
  self.dt = info['dt']
  self.dtt = self.dt#info['dtt']
  self.nsoil = len(info['dz'])#['nsoil']
  self.ncores = info['ncores']
  self.idate = info['idate']
  self.fdate = info['fdate']
  self.cdir = info['cdir']
  self.Qobs_file = info['Qobs_file']
  self.dt_timedelta = datetime.timedelta(seconds=self.dt)
  self.input_fp = nc.Dataset(info['input_file'])
  self.dx = self.input_fp.groups['metadata'].dx
  self.nhru = len(self.input_fp.dimensions['hru'])
  self.subsurface_module = info['subsurface_module']
  self.routing_module = info['routing_module']['type']
  self.routing_surface_coupling = info['routing_module']['surface_coupling']
  #self.hwu_flag = info['water_management']['hwu_flag']
  self.area = self.input_fp.groups['parameters'].variables['area'][:]
  self.pct = self.input_fp.groups['parameters'].variables['area_pct'][:]/100
  self.pct = self.pct/np.sum(self.pct)
  self.metadata = info
  self.ncsbasins=info['hmc_parameters']['number_of_characteristic_subbasins'] #laura
  self.m = self.input_fp.groups['parameters'].variables['m'][:]  #Noemi
  self.m[:] = 10.0 #m
  self.input_fp_meteo_time = self.input_fp.groups['meteorology'].variables['time']

  #Create a list of all the dates
  dates = []
  date = self.idate
  while date < self.fdate:
   dates.append(date)
   date = date + self.dt_timedelta
  self.dates = np.array(dates)

  #Read in all the meteorology for the given time period
  #determine the first time step for the meteorology
  var = self.input_fp['meteorology'].variables['time']
  ndates = var[:]
  #convert current date to num
  ndate = nc.date2num(self.idate,units=var.units,calendar=var.calendar)
  minitial_time = np.where(ndates == ndate)[0][0]
  #define period to extract
  idate = nc.date2num(self.idate,units=var.units,calendar=var.calendar)
  fdate = nc.date2num(self.fdate,units=var.units,calendar=var.calendar)
  idx  = np.where((ndates >= idate) & (ndates < fdate))[0]
  self.meteodata = {}
  for var in ['lwdown','swdown','psurf','wind','tair','spfh','precip']:
   self.meteodata[var] = self.input_fp['meteorology'][var][idx,:]
  #Define subtimestep
  dt = self.dt
  dt_meteo = (ndates[1] - ndates[0])*3600.0
  self.nt_timestep = int(dt_meteo/dt)

  return

 def initialize_noahmp(self,):

  #Initialize noahmp
  from model.pyNoahMP import NoahMP
  #self.noahmp = pyNoahMP.NoahMP
  self.noahmp = NoahMP

  #Initialize parameters
  self.noahmp.ncells = self.nhru
  self.noahmp.nsoil = self.nsoil
  self.noahmp.nsnow = 3
  self.noahmp.dt = self.dt

  #Allocate memory
  #1d,integer
  vars = ['isnow','isurban','slopetyp','soiltyp','vegtyp','ice','isc','ist','root_depth']
  for var in vars:
   exec('self.noahmp.%s = np.zeros(self.nhru).astype(np.int32)' % var)
  #1d,real
  vars = ['z_ml','lwdn','swdn','p_ml','psfc','prcp','t_ml','q_ml','u_ml','v_ml','fsh','ssoil','salb',\
          'fsno','swe','sndpth','emissi','qsfc1d','tv','tg','canice','canliq','eah','tah','cm','ch',\
          'fwet','sneqvo','albold','qsnow','wslake','zwt','dzwt','wa','wt','smcwtd','deeprech','rech',\
          'lfmass','rtmass','stmass','wood','stblcp','fastcp','plai','psai','tauss','t2mv','t2mb','q2mv',\
          'q2mb','trad','nee','gpp','npp','fvegmp','runsf','runsb','ecan','etran','esoil','fsa','fira',\
          'apar','psn','sav','sag','rssun','rssha','bgap','wgap','tgv','tgb','chv','chb','irc','irg','shc',\
          'shg','evg','ghv','irb','shb','evb','ghb','tr','evc','chleaf','chuc','chv2','chb2','cosz','lat',\
          'lon','fveg','fvgmax','fpice','fcev','fgev','fctr','qsnbot','ponding','ponding1','ponding2','fsr',\
          'co2pp','o2pp','foln','tbot','smcmax','smcdry','smcref','errwat','si0','si1','zwt0','minzwt','co2air',\
          'o2air','bb0','drysmc0','f110','maxsmc0','refsmc0','satpsi0','satdk0','satdw0','wltsmc0','qtz0',\
          'mozb','fvb','mozv','fvv','zpd','zpdg','tauxb','tauyb','tauxv','tauyv','sfcheadrt']
  for var in vars:
   exec('self.noahmp.%s = np.zeros(self.nhru).astype(np.float32)' % var)
  #2d,real
  #vars = ['stc','sh2o','smc','smceq','zsnso','snice','snliq','ficeold','zsoil','sldpth','hdiv']
  vars = ['sh2o','smc','smceq','zsoil','sldpth','hdiv']
  for var in vars:
   exec('self.noahmp.%s = np.zeros((self.nhru,self.nsoil),order=\'F\').astype(np.float32)' % var)
  self.noahmp.zsnso = np.zeros((self.nhru,self.noahmp.nsoil+self.noahmp.nsnow),order='F').astype(np.float32)
  self.noahmp.stc = np.zeros((self.nhru,self.noahmp.nsoil+self.noahmp.nsnow),order='F').astype(np.float32)
  self.noahmp.snice = np.zeros((self.nhru,self.noahmp.nsnow),order='F').astype(np.float32)
  self.noahmp.snliq = np.zeros((self.nhru,self.noahmp.nsnow),order='F').astype(np.float32)
  self.noahmp.ficeold = np.zeros((self.nhru,self.noahmp.nsnow),order='F').astype(np.float32)

  #Define data
  self.noahmp.llanduse = 'MODIFIED_IGBP_MODIS_NOAH'
  #assign_string(self.noahmp.llanduse.dtype,'MODIFIED_IGBP_MODIS_NOAH')
  self.noahmp.lsoil = 'CUST'
  #self.noahmp['lsoil'] = assign_string(self.noahmp.lsoil.dtype,'CUST')
  fdir = os.path.dirname(os.path.abspath(__file__))
  VEGPARM = '%s/pyNoahMP/data/VEGPARM.TBL' % fdir
  GENPARM = '%s/pyNoahMP/data/GENPARM.TBL' % fdir
  MPTABLE = '%s/pyNoahMP/data/MPTABLE.TBL' % fdir
  self.noahmp.vegparm_file = VEGPARM
  self.noahmp.genparm_file = GENPARM
  self.noahmp.mptable_file = MPTABLE
  #self.noahmp['vegparm_file'] = assign_string(self.noahmp.vegparm_file.dtype,VEGPARM)
  #self.noahmp['genparm_file'] = assign_string(self.noahmp.genparm_file.dtype,GENPARM)
  #self.noahmp['mptable_file'] = assign_string(self.noahmp.mptable_file.dtype,MPTABLE)
  #Define the options
  self.noahmp.idveg = 2 # dynamic vegetation (1 -> off ; 2 -> on)
  self.noahmp.iopt_crs = 2#1 #canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
  self.noahmp.iopt_btr = 1#1 # soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
  #Runoff 5 is really messed up
  self.noahmp.iopt_run = 2 # runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
  self.noahmp.iopt_sfc = 1#1#1 # surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
  self.noahmp.iopt_frz = 2#1#2 # supercooled liquid water (1-> NY06; 2->Koren99)
  self.noahmp.iopt_inf = 2#1#2 # frozen soil permeability (1-> NY06; 2->Koren99)
  self.noahmp.iopt_rad = 1#3#2 radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
  self.noahmp.iopt_alb = 1#2 snow surface albedo (1->BATS; 2->CLASS)
  self.noahmp.iopt_snf = 3#1#3 rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)]
  #self.noahmp.iopt_tbot = 1#2#1 # lower boundary of soil temperature (1->zero-flux; 2->Noah) 
  self.noahmp.iopt_tbot = 1#1 # lower boundary of soil temperature (1->zero-flux; 2->Noah) 
  self.noahmp.iopt_stc = 2#1#1#2 snow/soil temperature time scheme (only layer 1) 1 -> semi-implicit; 2 -> full implicit (original Noah)
  self.noahmp.iz0tlnd = -16497279 #0
  self.noahmp.sldpth[:] = np.array(self.metadata['dz'])
  
  self.noahmp.z_ml[:] = 10.0
  self.noahmp.zsoil = -np.cumsum(self.noahmp.sldpth[:],axis=1)
  self.noahmp.zsnso[:] = 0.0
  self.noahmp.zsnso[:,self.noahmp.nsnow:] = self.noahmp.zsoil[:]
  #Set all to land (can include lake in the future...)
  lc = self.input_fp.groups['parameters'].variables['land_cover'][:]
  ist = np.copy(lc).astype(np.int32)
  ist[:] = 1
  ist[lc == 17] = 2
  self.noahmp.ist[:] = ist[:]#1
  self.noahmp.isurban[:] = 13
  #Set soil color
  self.noahmp.isc[:] = 4
  self.noahmp.ice[:] = 0
  #Initialize snow info
  self.noahmp.isnow[:] = 0
  #Others
  self.noahmp.dx = self.dx
  self.noahmp.ice[:] = 0
  self.noahmp.foln[:] = 1.0
  self.noahmp.ficeold[:] = 0.0
  self.noahmp.albold[:] = 0.189999998 #0.65
  self.noahmp.sneqvo[:] = 2.09569975E-04 #0.0
  self.noahmp.ch[:] = 0.108068228
  self.noahmp.cm[:] = 1.81936752E-02
  self.noahmp.canliq[:] =  3.93530267E-04
  self.noahmp.canice[:] = 0.0
  self.noahmp.sndpth[:] = 1.06005312E-03
  self.noahmp.swe[:] = 2.09569975E-04
  self.noahmp.snice[:] = 0.0
  self.noahmp.snliq[:] = 0.0
  self.noahmp.wa[:] = 4900
  self.noahmp.wt[:] = self.noahmp.wa[:]
  self.noahmp.zwt[:] = (25.0 + 2.0) - self.noahmp.wa/1000.0/0.2#1.55999994
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
  self.noahmp.eah[:] = 373.016357 #It doesn't really matter. It's recomputed afterwards. Laura
  self.noahmp.fwet[:] = 0.0
  self.noahmp.psai[:] = 0.1
  stc_array=np.array([0,0,0,266.0995,274.0445,276.8954,279.9152])
  smc_array=np.array([0.2981597,0.2940254,0.2713114,0.3070948])
  self.noahmp.stc[:,0:7]=stc_array
  self.noahmp.smc[:,0:4]=smc_array
  if self.noahmp.nsoil>4:
   self.noahmp.stc[:,7:]=279.9152
   self.noahmp.smc[:,4:]=0.3070948
  self.noahmp.sh2o=self.noahmp.smc
  self.noahmp.slopetyp[:] = 3
  #self.noahmp.albold[:] = 0.5
  #Define the data
  self.noahmp.vegtyp[:] = self.input_fp.groups['parameters'].variables['land_cover'][:]
  hand = self.input_fp.groups['parameters'].variables['hand'][:]
  self.hrus = self.input_fp.groups['parameters'].variables['hru'][:].astype(np.int32)
  self.hbands = self.input_fp.groups['parameters'].variables['hband'][:].astype(np.int32)
  self.nhband = np.unique(self.hbands).size
  #print(self.nhband,flush=True)
  #mask = hand == 0
  #self.noahmp.vegtyp[mask] = 17#16
  self.noahmp.soiltyp[:] = np.arange(1,self.noahmp.ncells+1)
  self.noahmp.clay_pct  = self.input_fp.groups['parameters'].variables['clay'][:] # Noemi
  self.noahmp.smcmax[:] = self.input_fp.groups['parameters'].variables['MAXSMC'][:]
  self.noahmp.smcref[:] = self.input_fp.groups['parameters'].variables['REFSMC'][:]
  self.noahmp.smcdry[:] = self.input_fp.groups['parameters'].variables['DRYSMC'][:]
  #for ilayer in range(self.noahmp.sh2o.shape[1]):
   #self.noahmp.sh2o[:,ilayer] = self.input_fp.groups['parameters'].variables['MAXSMC'][:]
   #self.noahmp.sh2o[:,ilayer] = self.input_fp.groups['parameters'].variables['REFSMC'][:] #Noemi, start at REFSMC
  #self.noahmp.sh2o[:,:]=np.array([[0.298159689,0.294025391,0.271311402,0.307094812],[0.298159689,0.294025391,0.271311402,0.307094812]]) #Laura need to generalize for more layers
  #self.noahmp.smc[:] = self.noahmp.sh2o[:]
  #self.noahmp.smc[:,0]=0.298159689 #Laura
  self.noahmp.smcwtd[:] = self.noahmp.sh2o[:,0]
  #Initialize the soil parameters
  self.noahmp.bb0[:] = self.input_fp.groups['parameters'].variables['BB'][:]
  self.noahmp.drysmc0[:] = self.input_fp.groups['parameters'].variables['DRYSMC'][:]
  self.noahmp.f110[:] = 0.0#self.input_fp.groups['parameters'].variables['F11'][:]
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
  self.noahmp.tg[:] = 263.690887
  self.noahmp.tv[:] = 263.690887
  #for ilayer in range(self.noahmp.nsnow,self.noahmp.stc.shape[1]):
   #self.noahmp.stc[:,ilayer] = 285.0
  self.noahmp.tah[:] = 285.0
  self.noahmp.t2mv[:] = 0 #285
  self.noahmp.t2mb[:] = 0 #285
  self.noahmp.runsf[:] = 0.0
  self.noahmp.runsb[:] = 0.0
  #forcing
  self.noahmp.fveg[:] = 0.98
  self.noahmp.fvgmax[:] = 0.98
  self.noahmp.tbot[:] = 285.0

  #Define the parameters
  noah = self.noahmp
  self.noahmp.initialize_parameters(noah.llanduse,noah.lsoil,noah.vegparm_file,\
         noah.genparm_file,noah.iopt_crs,noah.iopt_btr,noah.iopt_run,noah.iopt_sfc,\
         noah.iopt_frz,noah.iopt_inf,noah.iopt_rad,noah.iopt_alb,noah.iopt_snf,\
         noah.iopt_tbot,noah.iopt_stc,noah.idveg,noah.mptable_file,noah.bb0,\
         noah.drysmc0,noah.f110,noah.maxsmc0,\
         noah.refsmc0,noah.satpsi0,noah.satdk0,noah.satdw0,noah.wltsmc0,noah.qtz0)#,\
         #noah.nsnow,noah.nsoil,noah.zsoil,noah.swe,noah.tg,noah.sndpth,noah.zsnso,\
         #noah.snice,noah.snliq,noah.isnow)

  return

 def initialize_routing(self,):

  #Determine what rank has which cid
  self.comm = self.MPI.COMM_WORLD
  self.size = self.comm.Get_size()
  self.rank = self.comm.Get_rank()
  if self.rank != 0:
   dest = 0
   db_ex = {self.cid:self.rank}
   self.comm.send(db_ex,dest=dest,tag=11)
  elif self.rank == 0:
   db = {}
   db[self.cid] = 0
   for i in range(1,self.size):
    db_ex = self.comm.recv(source=i,tag=11)
    db[list(db_ex.keys())[0]] = db_ex[list(db_ex.keys())[0]]
  #Wait until completed
  self.comm.Barrier()
  #Send the list to all the cores now
  if self.rank == 0:
   for i in range(1,self.size):
    self.comm.send(db,dest=i,tag=11)
  if self.rank != 0:
   db = self.comm.recv(source=0,tag=11)
  #Memorize links
  self.cid_rank_mapping = db

  if self.routing_module == 'kinematic':self.initialize_kinematic()

  return

 def initialize_kinematic(self,):

  from model.pyRouting import routing

  #Initialize kinematic wave routing
  self.routing = routing.kinematic(self.MPI,self.cid,self.cid_rank_mapping,self.dt,
                                   self.nhband,self.nhru,self.cdir,self.Qobs_file)
  self.routing.calculate_inundation_height_per_hband = routing.calculate_inundation_height_per_hband

  #Initialize hru to reach IRF
  ntmax = 100
  self.routing.IRF = {}
  uhs = self.routing.uhs[:]
  uhs = uhs/np.sum(uhs,axis=1)[:,np.newaxis]
  self.routing.IRF['uh'] = uhs[:]
  #Setup qfuture
  self.routing.IRF['qfuture'] = np.zeros((self.nhband,ntmax-1))

  return

 def update_routing(self,):

  if self.routing_module == 'kinematic':

   #Clean up runoff
   self.noahmp.runsf[self.noahmp.runsf < 0] = 0.0
   self.noahmp.runsb[self.noahmp.runsb < 0] = 0.0

   #Compute total runoff for the HRU
   runoff = self.noahmp.runsf+self.noahmp.runsb

   #Back out the actual runoff for this time step
   runoff_true = runoff - self.routing.hru_inundation*1000/self.dt #mm/s
   runoff_true[runoff_true < 0] = 0.0
   r = np.copy(runoff_true)
   r[:] = 1.0
   r[runoff > 0] = runoff_true[runoff > 0]/runoff[runoff > 0]
   self.noahmp.runsf[:] = r*self.noahmp.runsf[:]
   self.noahmp.runsb[:] = r*self.noahmp.runsb[:]
   self.routing.hru_runoff_inundation[:] = runoff - self.noahmp.runsf - self.noahmp.runsb
   self.routing.hru_runoff_inundation[self.routing.hru_runoff_inundation<0] = 0.0

   #Calculate runoff per hband (This needs to be optimized; Numba?)
   runoff_hband = np.zeros(self.nhband)
   area_hband = np.zeros(self.nhband)
   for hru in self.hrus:
    runoff_hband[self.hbands[hru]] += self.area[hru]*runoff[hru]
    area_hband[self.hbands[hru]] += self.area[hru]
   runoff_hband = runoff_hband/area_hband

   if self.routing_surface_coupling == True:
    #Zero out hband_inundation1
    self.routing.hband_inundation1[:] = 0.0

    #Over inundated hbands set runoff to be inundation and zero out corresponding runoff HERE
    m = self.routing.hband_inundation > 0
    self.routing.hband_inundation1[m] = self.dt*runoff_hband[m]/1000.0 #PROBLEMATIC
    runoff_hband[m] = 0.0

    #Determine the change per hband per reach
    scaling = self.routing.hband_inundation1[m]/(self.routing.fct_infiltrate*self.routing.hband_inundation[m])
    self.routing.reach2hband_inundation[:,m] = (1-self.routing.fct_infiltrate)*self.routing.reach2hband_inundation[:,m] + scaling*self.routing.fct_infiltrate*self.routing.reach2hband_inundation[:,m]
    #self.routing.reach2hband_inundation[:,m] = scaling*self.routing.reach2hband_inundation[:,m]

    #Determine the change per hband per reach (option 2)
    '''for ihband in range(self.routing.reach2hband_inundation.shape[1]):
     #Calculate total volume of water that was removed from the height band
     m = self.routing.reach2hband[:,ihband] > 0
     clength_hband = self.routing.c_length[m]
     Vloss = np.sum(clength_hband*(self.routing.hband_inundation1[ihband] - self.routing.hband_inundation[ihband]))
     #Split up the total volume of water among the reaches
     for ireach in range(self.routing.reach2hband_inundation.shape[0]):
      if m[ireach] == False:continue
      self.routing.reach2hband_inundation[ireach,ihband] = self.routing.reach2hband_inundation[ireach,ihband] - self.routing.c_length[ireach]/np.sum(clength_hband)*Vloss
      #Vloss += 
     print(ihband,np.sum(m),Vloss)'''

    #dif = self.routing.hband_inundation1[m] - self.routing.fct_infiltrate*self.routing.hband_inundation[m]
    #self.routing.reach2hband_inundation[:,m] = 0.5*self.routing.reach2hband_inundation[:,m]
    #self.routing.reach2hband_inundation[:,m] = self.routing.reach2hband_inundation[:,m] + dif
    #self.routing.reach2hband_inundation[(self.routing.reach2hband > 0) & (self.routing.reach2hband_inundation < 0.0)] = 0.0
    #tmp = self.routing.reach2hband_inundation[range(self.routing.reach2hband_inundation.shape[0]),self.routing.hband_channel]
    #tmp[tmp < 0.1] = 0.1
    #self.routing.reach2hband_inundation[range(self.routing.reach2hband_inundation.shape[0]),self.routing.hband_channel] = tmp[:]
    #print(np.unique(self.routing.reach2hband_inundation))
    

   #Recompute cross sectional area
   A = np.sum(self.routing.reach2hband*self.routing.reach2hband_inundation,axis=1)/self.routing.c_length
   Adb = self.routing.hdb['Ac'] + self.routing.hdb['Af']
   (self.routing.hband_inundation[:],self.routing.reach2hband_inundation[:]) = self.routing.calculate_inundation_height_per_hband(Adb,A,self.routing.hdb['W'],self.routing.hdb['M'],self.routing.hdb['hand'],self.routing.hdb['hband'].astype(np.int64),self.routing.reach2hband_inundation,self.routing.reach2hband)

   #Apply convolution to the runoff of each hband
   qr = runoff_hband[:,np.newaxis]*self.routing.IRF['uh']
   #QC to ensure we are conserving mass
   m = np.sum(qr,axis=1) > 0
   if np.sum(m) > 0:
    qr[m,:] = runoff_hband[m,np.newaxis]*qr[m,:]/np.sum(qr[m,:],axis=1)[:,np.newaxis]
   #Add corresponding qfuture to this time step
   crunoff = self.routing.IRF['qfuture'][:,0]
   #Update qfuture
   self.routing.IRF['qfuture'][:,0:-1] = self.routing.IRF['qfuture'][:,1:]
   self.routing.IRF['qfuture'][:,-1] = 0.0
   #Add current time step
   crunoff = crunoff + qr[:,0]
   #Add the rest to qfuture
   self.routing.IRF['qfuture'] = self.routing.IRF['qfuture'] + qr[:,1:]

   #Aggregate the hband runoff at the reaches
   self.routing.qss[:] += self.routing.reach2hband.dot(crunoff/1000.0)/self.routing.c_length #m2/s
 
   if self.routing_surface_coupling == True:
    #Add changes between Ac1 and Ac0 to qss term
    self.routing.qss[:] += (A - self.routing.A0[:])/self.dt
    #self.routing.A0[:] = A[:]

   #Update routing module
   self.routing.itime = self.itime
   self.routing.update(self.dt)

   #Update the hru inundation values
   for hru in self.hrus:
    #Only allow fct of the inundated height to infiltrate
    self.routing.hru_inundation[hru] = self.routing.hband_inundation[self.hbands[hru]]
   
  return

 def initialize_subsurface(self,):

  if self.subsurface_module == 'richards':self.initialize_richards(self.ncsbasins)
  #if self.subsurface_module == 'richards':self.initialize_richards()
  return

 def initialize_richards(self,ncsbasins):
 #def initialize_richards(self):

  from model.pyRichards import richards
  #Initialize richards characteristic subbasins based, laura
  self.richards=richards.richards(self.nhru,self.nhband,self.nsoil)
  self.richards.dx = {} #self.dx
  self.richards.nhru = self.nhru
  self.richards.nhband = int(self.nhband) #laura
  self.richards.ncsbasins=int(ncsbasins) #number of characteristic subbasins
  self.richards.dem[:] = self.input_fp.groups['parameters'].variables['hand'][:]
  hru_area = self.input_fp.groups['parameters'].variables['area'][:] #laura
  basin_aux=0
  w_dict={}

  for i in range(1,(int(ncsbasins)+1)):
   text='wmatrix_Basin%s'%int(i)
   text2='Basin%s' %int(i)
   print(' ',text2)
   a=sparse.csr_matrix((self.input_fp.groups[text].variables['data'][:],
                                  self.input_fp.groups[text].variables['indices'][:],
                                  self.input_fp.groups[text].variables['indptr'][:]),
                                  shape=(int((np.max(self.input_fp.groups[text].variables['indices']))+1),int((np.max(self.input_fp.groups[text].variables['indices']))+1)),dtype=np.float64)
   w_dict[text2]=np.array(a.todense())
   self.richards.area=np.empty(self.richards.nhband)  

   aux=0
   for h_band in (np.unique(self.hbands)): #laura sum up and average up all the HRUs areas that belong to the same hband
    m=self.hbands==h_band
    self.richards.area[aux]=np.sum(hru_area[m])
    self.richards.m[aux]=np.sum(((self.input_fp.groups['parameters'].variables['m'][m])*(hru_area[m]))/(self.richards.area[aux]))
    self.richards.demhband[aux]=np.sum(((self.input_fp.groups['parameters'].variables['hand'][m])*(hru_area[m]))/(self.richards.area[aux]))
    self.richards.slope[aux]=np.sum(((self.input_fp.groups['parameters'].variables['slope'][m])*(hru_area[m]))/(self.richards.area[aux]))
    aux=aux+1

   with np.errstate(invalid='ignore', divide='ignore'):
    tmp = (self.richards.area[basin_aux:int((np.max(self.input_fp.groups[text].variables['indices']))+1)+basin_aux])/w_dict['Basin%s' %int(i)]
   dx = (tmp + tmp.T)/2
   self.richards.dx[text2] = dx
   basin_aux=basin_aux+w_dict[text2].shape[0]
  self.richards.w=w_dict

  #Initialize richards, only height band based, laura
  #hru_area = self.input_fp.groups['parameters'].variables['area'][:] #laura
  #self.richards=richards.richards(self.nhru,self.nhband,self.nsoil) #laura
  #Set other parameters
  #self.richards.dx = self.dx
  #self.richards.nhru = self.nhru
  #self.richards.nhband = self.nhband #laura
  #self.richards.dem[:] = self.input_fp.groups['parameters'].variables['hand'][:]
  #self.richards.width = sparse.csr_matrix((self.input_fp.groups['wmatrix'].variables['data'][:],
  #                                self.input_fp.groups['wmatrix'].variables['indices'][:],
  #                                self.input_fp.groups['wmatrix'].variables['indptr'][:]),
  #                                shape=(self.nhband,self.nhband),dtype=np.float64) #laura replaced nhru with nhband
  #self.richards.I = self.richards.width.copy()
  #self.richards.I[self.richards.I != 0] = 1
  #self.richards.w = np.array(self.richards.width.todense())
  #self.richards.area=np.empty([self.richards.w.shape[0]])
  #aux=0
  #for h_band in np.unique(self.hbands): #laura sum up and average up all the HRUs areas that belong to the same hband
   #m=self.hbands==h_band
   #self.richards.area[aux]=np.sum(hru_area[m])
   #self.richards.m[aux]=np.sum(((self.input_fp.groups['parameters'].variables['m'][m])*(hru_area[m]))/(self.richards.area[aux]))
   #self.richards.demhband[aux]=np.sum(((self.input_fp.groups['parameters'].variables['hand'][m])*(hru_area[m]))/(self.richards.area[aux]))
   #self.richards.slope[aux]=np.sum(((self.input_fp.groups['parameters'].variables['slope'][m])*(hru_area[m]))/(self.richards.area[aux]))
   #aux=aux+1
  #with np.errstate(invalid='ignore', divide='ignore'):
   #tmp = self.richards.area/self.richards.w
  #dx = (tmp + tmp.T)/2
  #self.richards.dx = dx
  return

 def initialize_hwu(self,info):

  #from pyHWU.Human_Water_Use import Human_Water_Use as hwu
  #self.hwu = hwu(self,info)
  #self.hwu.area = self.input_fp.groups['parameters'].variables['area'][:]

  #if self.hwu.hwu_flag == True:
  # print("Initializing Human Water Management")
  # self.hwu.initialize_allocation(self,)

  return

 def run(self,info):

  #Run the model
  date = self.idate
  tic = time.time()
  self.noahmp.dzwt[:] = 0.0
  i=0

  while date < self.fdate:
   i=i+1 
   #Update input data
   tic0 = time.time()
   self.noahmp.stc[:,0:self.noahmp.nsnow]=self.tsno #Laura
   self.update_input(date)
   if i==1: #Laura. Initial conditions for tah and eah as defined in NOAH
    self.noahmp.tah=self.noahmp.t_ml
    self.noahmp.eah=(self.noahmp.p_ml*self.noahmp.q_ml)/(0.622+self.noahmp.q_ml)
   #Save the original precip
   precip = np.copy(self.noahmp.prcp)

   #Calculate initial NOAH water balance
   self.initialize_water_balance()

   #Update model
   tic0 = time.time()
   self.update(date)

   #print('update model',time.time() - tic0,flush=True)
   #exit()
   
   #Return precip to original value
   #self.noahmp.prcp[:] = precip[:] 

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
   tic0 = time.time()
   self.update_output(date)
   #print('update output',time.time() - tic0,flush=True)

   #Update time step
   date = date + self.dt_timedelta
   self.itime = self.itime + 1
   #print(date,flush=True)
   #exit()

   #Output some statistics
   if (date.hour == 0) and (date.day == 1):
    #print('cid: %d' % info['cid'],date.strftime("%Y-%m-%d"),'%10.4f'%(time.time()-tic),'et:%10.4f'%self.et,'prcp:%10.4f'%self.prcp,'q:%10.4f'%self.q,'WB ERR:%10.6f' % self.errwat,'ENG ERR:%10.6f' % self.erreng,flush=True)
    print(date.strftime("%Y-%m-%d"),'%10.4f'%(time.time()-tic),'et:%10.4f'%self.et,'prcp:%10.4f'%self.prcp,'q:%10.4f'%self.q,'WB ERR:%10.6f' % self.errwat,'ENG ERR:%10.6f' % self.erreng,flush=True)

  return

 def update_input(self,date):

  self.noahmp.itime = self.itime
  dt = self.dt
  self.noahmp.nowdate = assign_string(19,date.strftime('%Y-%m-%d_%H:%M:%S'))
  #self.noahmp.nowdate = date.strftime('%Y-%m-%d_%H:%M:%S')
  self.noahmp.julian = (date - datetime.datetime(date.year,1,1,0)).days
  self.noahmp.yearlen = (datetime.datetime(date.year+1,1,1,0) - datetime.datetime(date.year,1,1,1,0)).days + 1

  #Update meteorology
  #meteorology = self.input_fp.groups['meteorology']
  '''if date == self.idate:
   #determine the first time step for the meteorology
   var = meteorology.variables['time']
   ndates = var[:]
   #convert current date to num
   ndate = nc.date2num(date,units=var.units,calendar=var.calendar)
   #identify the position
   self.minitial_itime = np.where(ndates == ndate)[0][0]
  i = self.itime + self.minitial_itime'''

  '''
  #Update meteorology
  meteorology = self.input_fp.groups['meteorology']
  #convert current date to num
  ndate = nc.date2num(date,units=self.input_fp_meteo_time.units,calendar=self.input_fp_meteo_time.calendar)
  #identify the position
  i = np.argmin(np.abs(self.input_fp_meteo_time[:]-ndate))

  #Downscale meteorology
  self.noahmp.lwdn[:] = meteorology.variables['lwdown'][i,:] #W/m2
  self.noahmp.swdn[:] = meteorology.variables['swdown'][i,:] #W/m2
  self.noahmp.psfc[:] = meteorology.variables['psurf'][i,:] #Pa
  self.noahmp.p_ml[:] = self.noahmp.psfc[:] #Pa
  self.noahmp.u_ml[:] = (meteorology.variables['wind'][i,:]**2/2)**0.5 #m/s
  self.noahmp.v_ml[:] = self.noahmp.u_ml[:] #m/s
  self.noahmp.t_ml[:] = meteorology.variables['tair'][i,:] #K
  self.noahmp.q_ml[:] = meteorology.variables['spfh'][i,:] #Kg/Kg
  self.noahmp.qsfc1d[:] = self.noahmp.q_ml[:] #Kg/Kg
  self.noahmp.prcp[:] = meteorology.variables['precip'][i,:] #mm/s
  '''
  
  #Update meteorology
  i = int(np.floor(self.itime/self.nt_timestep))
  self.noahmp.lwdn[:] = self.meteodata['lwdown'][i,:] #W/m2
  self.noahmp.swdn[:] = self.meteodata['swdown'][i,:] #W/m2
  self.noahmp.psfc[:] = self.meteodata['psurf'][i,:] #Pa
  self.noahmp.p_ml[:] = self.noahmp.psfc[:] #Pa
  self.noahmp.u_ml[:] = (self.meteodata['wind'][i,:]**2/2)**0.5 #m/s
  self.noahmp.v_ml[:] = self.noahmp.u_ml[:] #m/s
  self.noahmp.t_ml[:] = self.meteodata['tair'][i,:] #K
  self.noahmp.q_ml[:] = self.meteodata['spfh'][i,:] #Kg/Kg
  self.noahmp.qsfc1d[:] = self.noahmp.q_ml[:] #Kg/Kg
  self.noahmp.prcp[:] = self.meteodata['precip'][i,:] #mm/s

  #Set the partial pressure of CO2 and O2
  self.noahmp.co2air[:] = 355.E-6*self.noahmp.psfc[:]# ! Partial pressure of CO2 (Pa) ! From NOAH-MP-WRF
  self.noahmp.o2air[:] = 0.209*self.noahmp.psfc[:]# ! Partial pressure of O2 (Pa)  ! From NOAH-MP-WRF

  for ilayer in range(self.noahmp.nsnow):
   self.noahmp.ficeold[:,ilayer]=self.noahmp.snice[:,ilayer]/(self.noahmp.snice[:,ilayer]+self.noahmp.snliq[:,ilayer])

  # Update water demands
  '''if self.hwu.hwu_flag == True:
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
     self.i    hwu.deficit_lstock[:] = np.copy(self.hwu.demand_lstock[:])'''


  return

 def update(self,date):

  if self.routing_surface_coupling == True:
   self.noahmp.sfcheadrt[:] = self.routing.fct_infiltrate*self.routing.hru_inundation[:]*1000 #mm
  else:
   self.noahmp.sfcheadrt[:] = 0.0
  
  # Update subsurface
  self.update_subsurface()

  # Update NOAH
  noah = self.noahmp
  self.noahmp.run_model(self.ncores,noah.yearlen,noah.idveg,noah.iopt_crs,noah.iopt_btr,\
                        noah.iopt_run,noah.iopt_sfc,noah.iopt_frz,noah.iopt_inf,noah.iopt_rad,\
                        noah.iopt_alb,noah.iopt_snf,noah.iopt_tbot,noah.iopt_stc,noah.itime,\
                        noah.iz0tlnd,noah.dt,noah.dx,noah.julian,noah.isnow,noah.isurban,\
                        noah.slopetyp,noah.soiltyp,noah.vegtyp,noah.ice,noah.isc,noah.ist,\
                        noah.root_depth,noah.nowdate,noah.z_ml,noah.lwdn,noah.swdn,noah.p_ml,\
                        noah.psfc,noah.prcp,noah.t_ml,noah.q_ml,noah.u_ml,noah.v_ml,noah.fsh,\
                        noah.ssoil,noah.salb,noah.fsno,noah.swe,noah.sndpth,noah.emissi,\
                        noah.qsfc1d,noah.tv,noah.tg,noah.canice,noah.canliq,noah.eah,\
                        noah.tah,noah.cm,noah.ch,noah.fwet,noah.sneqvo,noah.albold,noah.qsnow,\
                        noah.wslake,noah.zwt,noah.dzwt,noah.wa,noah.wt,noah.smcwtd,noah.deeprech,\
                        noah.rech,noah.lfmass,noah.rtmass,noah.stmass,noah.wood,noah.stblcp,\
                        noah.fastcp,noah.plai,noah.psai,noah.tauss,noah.t2mv,noah.t2mb,noah.q2mv,\
                        noah.q2mb,noah.trad,noah.nee,noah.gpp,noah.npp,noah.fvegmp,noah.runsf,\
                        noah.runsb,noah.ecan,noah.etran,noah.esoil,noah.fsa,noah.fira,noah.apar,\
                        noah.psn,noah.sav,noah.sag,noah.rssun,noah.rssha,noah.bgap,noah.wgap,\
                        noah.tgv,noah.tgb,noah.chv,noah.chb,noah.irc,noah.irg,noah.shc,noah.shg,\
                        noah.evg,noah.ghv,noah.irb,noah.shb,noah.evb,noah.ghb,noah.tr,noah.evc,\
                        noah.chleaf,noah.chuc,noah.chv2,noah.chb2,noah.cosz,noah.lat,noah.lon,\
                        noah.fveg,noah.fvgmax,noah.fpice,noah.fcev,noah.fgev,noah.fctr,noah.qsnbot,\
                        noah.ponding,noah.ponding1,noah.ponding2,noah.fsr,noah.co2pp,noah.o2pp,\
                        noah.foln,noah.tbot,noah.smcmax,noah.smcdry,noah.smcref,noah.errwat,noah.si0,\
                        noah.si1,noah.zwt0,noah.minzwt,noah.co2air,noah.o2air,noah.bb0,noah.drysmc0,\
                        noah.f110,noah.maxsmc0,noah.refsmc0,noah.satpsi0,noah.satdk0,noah.satdw0,\
                        noah.wltsmc0,noah.qtz0,noah.mozb,noah.fvb,noah.mozv,noah.fvv,\
                        noah.zpd,noah.zpdg,noah.tauxb,noah.tauyb,noah.tauxv,noah.tauyv,\
                        noah.stc,noah.sh2o,noah.smc,noah.smceq,noah.zsnso,\
                        noah.snice,noah.snliq,noah.ficeold,noah.zsoil,noah.sldpth,noah.hdiv,
                        noah.sfcheadrt)
  self.tsno=self.noahmp.stc[:,0:self.noahmp.nsnow] #Laura
  # Calculate water demands and supplies, and allocate volumes
  #self.hwu.Calc_Human_Water_Demand_Supply(self,date)

  # Abstract Surface Water and Groundwater
  #self.hwu.Water_Supply_Abstraction(self,date)

  # Update routing
  self.update_routing()

  return

 def update_subsurface(self,):

  self.noahmp.dzwt[:] = 0.0

  if self.subsurface_module == 'richards':
   #Assign noahmp variables to subsurface module
   hru_area = self.input_fp.groups['parameters'].variables['area'][:] #laura
   aux=0
   for h_band in np.unique(self.hbands): #laura average up all the HRUs areas that belong to the same hband
    m=self.hbands==h_band
    self.richards.theta[aux,:]=(np.sum(((self.noahmp.smc[m,:])*(hru_area[m])[:,None]),axis=0))/(self.richards.area[aux])
    self.richards.thetar[aux]=(np.sum(((self.noahmp.drysmc0[m])*(hru_area[m])),axis=0))/(self.richards.area[aux])
    self.richards.thetas[aux]=(np.sum(((self.noahmp.maxsmc0[m])*(hru_area[m])),axis=0))/(self.richards.area[aux])
    self.richards.b[aux]=(np.sum(((self.noahmp.bb0[m])*(hru_area[m])),axis=0))/(self.richards.area[aux])
    self.richards.satpsi[aux]=(np.sum(((self.noahmp.satpsi0[m])*(hru_area[m])),axis=0))/(self.richards.area[aux])
    self.richards.ksat[aux]=(np.sum(((self.noahmp.satdk0[m])*(hru_area[m])),axis=0))/(self.richards.area[aux])
    self.richards.dz[aux,:]=(np.sum(((self.noahmp.sldpth[m,:])*(hru_area[m])[:,None]),axis=0))/(self.richards.area[aux])
    aux=aux+1

#   print('theta: ',self.richards.theta)
#   print('thetar: ',self.richards.thetar)
#   print('thetas: ',self.richards.thetas)
#   print('b: ',self.richards.b)
#   print('ksat: ',self.richards.ksat)
#   print('satpsi: ',self.richards.satpsi)
#   print('dz: ', self.richards.dz)
#   import sys
#   sys.exit('PARAR!')
   #self.richards.theta[:] = self.noahmp.smc[:]
   #self.richards.thetar[:] = self.noahmp.drysmc0[:]
   #self.richards.thetas[:] = self.noahmp.maxsmc0[:]
   #self.richards.b[:] = self.noahmp.bb0[:]
   #self.richards.satpsi[:] = self.noahmp.satpsi0[:]
   #self.richards.ksat[:] = self.noahmp.satdk0[:]
   #self.richards.dz[:] = self.noahmp.sldpth[:]

   #Update subsurface module
   #0.Update hand value to account for hru inundation (This is a hack to facilitate a non-flooding stream to influence its surrounding hrus)
   if self.routing_module == 'kinematic': #laura kept the inundation computation at a HRU level
     self.richards.dem1 = self.richards.dem+self.routing.hru_inundation
     m = self.richards.dem1[0:-1] > self.richards.dem1[1:]
     self.richards.dem1[0:-1][m] = self.richards.dem1[1:][m]
   else:
     self.richards.dem1 = self.richards.dem
   aux=0
   for h_band in np.unique(self.hbands): #laura average up all the HRUs areas that belong to the same hband
    m=self.hbands==h_band
    self.richards.dem1hband[aux]=(np.sum(((self.richards.dem1[m])*(hru_area[m])),axis=0))/(self.richards.area[aux])
    aux=aux+1
   #self.richards.update()
   self.richards.update_numba()
   #print(self.noahmp.hdiv.shape)
   #Assign subsurface module variables to noahmp
   aux=0 #laura convert hband level divergence to HRU level. TO DO: limit the dovergence for different land covers 
   for h_band in np.unique(self.hbands):
    m=self.hbands==h_band
    self.noahmp.hdiv[m,:]=self.richards.hdiv[aux,:]
    aux=aux+1
   #self.noahmp.hdiv[:] = self.richards.hdiv[:]
   
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
  #HWU = self.hwu
  dt = self.dt
  if self.subsurface_module == 'richards':
   tmp = np.copy(self.end_wb - self.beg_wb - NOAH.dt*(NOAH.prcp-NOAH.ecan-
         NOAH.etran-NOAH.esoil-NOAH.runsf-NOAH.runsb-np.sum(NOAH.hdiv,axis=1)))
  else:
   tmp = np.copy(self.end_wb - self.beg_wb - NOAH.dt*(NOAH.prcp-NOAH.ecan-
         NOAH.etran-NOAH.esoil-NOAH.runsf-NOAH.runsb))
  if self.routing_module == 'kinematic':
   tmp = tmp - np.copy(NOAH.sfcheadrt) + np.copy(self.routing.hru_runoff_inundation)*dt
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
  #HWU = self.hwu
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
  tmp['lwnet'] = np.copy(NOAH.fira) #W/m2
  tmp['swnet'] = np.copy(NOAH.sav + NOAH.sag) #W/m2
  tmp['shb'] = np.copy(NOAH.shb) #W/m2
  tmp['shc'] = np.copy(NOAH.shc) #W/m2
  tmp['shg'] = np.copy(NOAH.shg) #W/m2
  tmp['evb'] = np.copy(NOAH.evb) #W/m2
  tmp['evc'] = np.copy(NOAH.evc) #W/m2
  tmp['evg'] = np.copy(NOAH.evg) #W/m2
  tmp['swdn'] = np.copy(NOAH.swdn) #W/m2
  tmp['t2mv'] = np.copy(NOAH.t2mv) #W/m2
  tmp['t2mb'] = np.copy(NOAH.t2mb) #W/m2
  tmp['qbase'] = NOAH.dt*np.copy(NOAH.runsb) #mm
  tmp['qsurface'] = NOAH.dt*np.copy(NOAH.runsf) #mm
  tmp['runoff'] = NOAH.dt*(np.copy(NOAH.runsf)+np.copy(NOAH.runsb)) #mm 
  tmp['prcp'] = NOAH.dt*np.copy(NOAH.prcp) #mm
  tmp['trad'] = np.copy(NOAH.trad) #K
  tmp['stc'] = np.copy(NOAH.stc[:,3]) #K
  tmp['tv'] = np.copy(NOAH.tv) #K
  tmp['salb'] = np.copy(NOAH.salb) 
  tmp['wtd'] = np.copy(NOAH.zwt) #m
  tmp['swe'] = np.copy(NOAH.swe) #m
  tmp['totsmc'] = np.sum(NOAH.sldpth*NOAH.smc,axis=1)/np.sum(NOAH.sldpth[0]) #m3/m3
  tmp['hdiv'] = np.copy(NOAH.hdiv)
  tmp['mozb'] = np.copy(NOAH.mozb)
  tmp['mozv'] = np.copy(NOAH.mozv)
  tmp['fvb'] = np.copy(NOAH.fvb)
  tmp['fvv'] = np.copy(NOAH.fvv)
  tmp['fveg'] = np.copy(NOAH.fvegmp)
  tmp['zpd'] = np.copy(NOAH.zpd)
  tmp['zpdg'] = np.copy(NOAH.zpdg)
  tmp['tauxb'] = np.copy(NOAH.tauxb)
  tmp['tauyb'] = np.copy(NOAH.tauyb)
  tmp['tauxv'] = np.copy(NOAH.tauxv)
  tmp['tauyv'] = np.copy(NOAH.tauyv)
  tmp['cm'] = np.copy(NOAH.cm)
  tmp['ch'] = np.copy(NOAH.ch)
  tmp['tg'] = np.copy(NOAH.tg)
  #routing
  if self.routing_module == 'kinematic':
   tmp['sfcheadrt'] = np.copy(NOAH.sfcheadrt)
   tmp['inundation'] = np.copy(self.routing.hru_inundation)


  # root zone
  cs = np.cumsum(NOAH.sldpth[0,:])
  mask = cs <= 0.5
  pct = NOAH.sldpth[0,mask]/np.sum(NOAH.sldpth[0,mask])
  tmp['smc_root'] = np.sum(pct*NOAH.smc[:,mask],axis=1) #m3/m3
  # top soil layer
  tmp['smc1'] = np.copy(NOAH.smc[:,0])
  #print(tmp['smc1'])
  # total soil moisture / soil water storage -- only until depth to bedrock
  #dl = [ np.argmin(np.abs(cs-self.m[i])) for i in range(self.nhru)]
  #tmp['totsmc'] = [ np.sum(NOAH.sldpth[i,:dl[i]]*NOAH.smc[i,:dl[i]])/np.sum(NOAH.sldpth[i,:dl[i]])  for i in range(self.nhru)]
  #print(tmp['totsmc'])
  #exit()
  
  #General
  tmp['errwat'] = np.copy(HB.errwat)

  # Water Management
  '''if self.hwu.hwu_agric_flag == True:
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
   if self.hwu.hwu_gw_flag == True:                                                                                              tmp['alloc_gw'] = np.copy(HWU.alloc_gw) #m'''

  #Output the variables
  sep = 100
  if itime == 0:
   self.output = {}
   for var in self.metadata['output']['vars']:
    shp = grp.variables[var].shape
    self.output[var] = np.zeros((sep,shp[1]))
  #Fill the data (MISSING!)
  val = itime % sep
  for var in self.metadata['output']['vars']:
    self.output[var][val,:] = tmp[var]
  # self.output[itime] =
  if (itime+1) % sep == 0:
   for var in self.metadata['output']['vars']:
    grp.variables[var][itime-sep+1:itime+1,:] = self.output[var][:]
  if (itime+1) == self.ntime:
   val = int(sep*np.ceil((itime - sep)/sep))
   for var in self.metadata['output']['vars']:
    grp.variables[var][val:itime+1,:] = self.output[var][0:itime-val+1,:]

  #Output routing variables
  if self.routing_module == 'kinematic':
   grp = self.output_fp.groups['data_routing']
   tmp = {}
   tmp['A'] = self.routing.A1[:]
   tmp['Q'] = self.routing.Q1[:]
   tmp['Qc'] = self.routing.Qc[:]
   tmp['Qf'] = self.routing.Qf[:]
   tmp['reach_inundation'] = self.routing.reach2hband_inundation[:]
   sep = 100
   if itime == 0:
    self.output_routing = {}
    for var in self.metadata['output']['routing_vars']:
     shp = grp.variables[var].shape
     self.output_routing[var] = np.zeros((sep,shp[1]))
   #Fill the data (MISSING!)
   val = itime % sep
   for var in self.metadata['output']['routing_vars']:
     self.output_routing[var][val,:] = tmp[var]
   # self.output[itime] =
   if (itime+1) % sep == 0:
    for var in self.metadata['output']['routing_vars']:
     grp.variables[var][itime-sep+1:itime+1,:] = self.output_routing[var][:]
   if (itime+1) == self.ntime:
    val = int(sep*np.ceil((itime - sep)/sep))
    for var in self.metadata['output']['routing_vars']:
     grp.variables[var][val:itime+1,:] = self.output_routing[var][0:itime-val+1,:]

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
             'lwnet':{'description':'Net longwave radiation','units':'W/m2','dims':('time','hru',),'precision':4},
             'swnet':{'description':'Absorbed shortwave radiation','units':'W/m2','dims':('time','hru',),'precision':4},
             "t2mv":{'description':'Vegetated air temperature','units':'W/m2','dims':('time','hru',),'precision':4},
	     "t2mb":{'description':'Bare air temperature','units':'W/m2','dims':('time','hru',),'precision':4},
	     "fveg":{'description':'Vegetated fraction','units':'','dims':('time','hru',),'precision':4},
	     "mozb":{'description':'Bare stability parameter','units':'','dims':('time','hru',),'precision':4},
	     "mozv":{'description':'Vegetated stability parameter','units':'','dims':('time','hru',),'precision':4},
	     "zpd":{'description':'Vegetated zero plane displacement','units':'','dims':('time','hru',),'precision':4},
	     "zpdg":{'description':'Ground zero plane displacement','units':'','dims':('time','hru',),'precision':4},
             "cm":{'description':'Momentum drag coefficient','units':'','dims':('time','hru',),'precision':4},
             "ch":{'description':'Sensible heat exchange coefficient','units':'','dims':('time','hru',),'precision':4},
	     "tauxb":{'description':'Bare wind stress (e-w)','units':'','dims':('time','hru',),'precision':4},
	     "tauyb":{'description':'Bare wind stress (n-s)','units':'','dims':('time','hru',),'precision':4},
	     "tauxv":{'description':'Bare wind stress (e-w)','units':'','dims':('time','hru',),'precision':4},
	     "tauyv":{'description':'Bare wind stress (n-s)','units':'','dims':('time','hru',),'precision':4},
	     "fvb":{'description':'Bare friction velocity','units':'m/s','dims':('time','hru',),'precision':4},
	     "fvv":{'description':'Vegetated friction velocity','units':'m/s','dims':('time','hru',),'precision':4},
             "shb":{'description':'Bare sensible heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             "shc":{'description':'Vegetated canopy sensible heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             "shg":{'description':'Vegetated ground sensible heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             "evb":{'description':'Bare latent heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             "evc":{'description':'Vegetated canopy latent heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             "evg":{'description':'Vegetated ground latent heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             "swdn":{'description':'Shortwave down','units':'W/m2','dims':('time','hru',),'precision':4},
             'trad':{'description':'Land surface skin temperature','units':'K','dims':('time','hru',),'precision':2},
             'stc':{'description': 'Snow/soil temperature','units':'K','dims':('time','hru',),'precision':2},
             'tv':{'description': 'Canopy temperature','units':'K','dims':('time','hru',),'precision':2},
             'tg':{'description': 'Ground temperature','units':'K','dims':('time','hru',),'precision':2},
             'salb':{'description':'Land surface albedo','units':' ','dims':('time','hru',),'precision':4},

             'qbase':{'description':'Excess runoff','units':'mm/s','dims':('time','hru',),'precision':4},
             'qsurface':{'description':'Surface runoff','units':'mm/s','dims':('time','hru',),'precision':4},
             'runoff':{'description':'Total runoff','units':'mm/s','dims':('time','hru',),'precision':4},
             'prcp':{'description':'Precipitation','units':'mm/s','dims':('time','hru',),'precision':4},
             'smc':{'description':'Soil water content','units':'m3/m3','dims':('time','hru','soil'),'precision':3},
             'swd':{'description':'Soil water deficit','units':'mm','dims':('time','hru',),'precision':4},
             'sstorage':{'description':'Surface storage','units':'mm','dims':('time','hru',),'precision':4},
             #'sndpth':{'description':'Snow depth','units':'mm','dims':('time','hru',)},
             'swe':{'description':'Snow water equivalent','units':'mm','dims':('time','hru',),'precision':4},

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
             'alloc_gw':{'description':'Groundwater water allocated','units':'m','dims':('time','hru',),'precision':4},
             'Q':{'description':'Discharge','units':'m3/s','dims':('time','channel',),'precision':4},
             'Qc':{'description':'Discharge (channel)','units':'m3/s','dims':('time','channel',),'precision':4},
             'Qf':{'description':'Discharge (floodplain)','units':'m3/s','dims':('time','channel',),'precision':4},
             'reach_inundation':{'description':'Inundation height (reach level)','units':'m','dims':('time','channel','hband'),'precision':3},
             'A':{'description':'Cross section','units':'m2','dims':('time','channel',),'precision':4},
             'inundation':{'description':'Inundation height','units':'m','dims':('time','hru',),'precision':4},
             }

  #Create the dimensions
  print('Creating the dimensions',flush=True)
  ntime = 24*3600*((self.fdate - self.idate).days)/self.dt
  self.ntime = ntime
  nhru = len(fp_in.dimensions['hru'])
  fp_out.createDimension('hru',nhru)
  fp_out.createDimension('time',ntime)
  fp_out.createDimension('soil',self.nsoil)
  fp_out.createDimension('hband',self.nhband)
  if self.routing_module == 'kinematic':
   fp_out.createDimension('channel',self.routing.nchannel)

  #Create the output
  print('Creating the data group',flush=True)
  grp = fp_out.createGroup('data')
  for var in self.metadata['output']['vars']:
   ncvar = grp.createVariable(var,'f4',metadata[var]['dims'],least_significant_digit=metadata[var]['precision'])#,zlib=True)
   ncvar.description = metadata[var]['description']
   ncvar.units = metadata[var]['units']

  #Create the routing output
  if self.routing_module == 'kinematic':
   print('Creating the routing group',flush=True)
   grp = fp_out.createGroup('data_routing')
   for var in self.metadata['output']['routing_vars']:
    ncvar = grp.createVariable(var,'f4',metadata[var]['dims'],least_significant_digit=metadata[var]['precision'])#,zlib=True)
    ncvar.description = metadata[var]['description']
    ncvar.units = metadata[var]['units']


  #Create the metadata
  print('Creating the metadata group',flush=True)
  grp = fp_out.createGroup('metadata')
  #Time
  grp.createVariable('time','f8',('time',))
  dates = grp.createVariable('date','f8',('time',))
  dates.units = 'hours since 1900-01-01'
  dates.calendar = 'standard'

  #HRU percentage coverage
  print('Setting the HRU percentage coverage',flush=True)
  pcts = grp.createVariable('pct','f4',('hru',))
  pcts[:] = fp_in.groups['parameters'].variables['area_pct'][:]
  pcts.description = 'hru percentage coverage'
  pcts.units = '%'

  #HRU area
  print('Setting the HRU areal coverage',flush=True)
  area = grp.createVariable('area','f4',('hru',))
  area[:] = fp_in.groups['parameters'].variables['area'][:]
  area.units = 'meters squared'
  print('Defining the HRU ids',flush=True)
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
  #routing
  if self.routing_module == 'kinematic':	
   fp['Q0'] = self.routing.Q0[:]
   fp['u0'] = self.routing.u0[:]
   fp['A0'] = self.routing.A0[:]
   fp['Q1'] = self.routing.Q1[:]
   fp['A1'] = self.routing.A1[:]
   fp['dA'] = self.routing.dA[:]
   fp['bcs'] = self.routing.bcs[:]
   fp['qss'] = self.routing.qss[:]
   #fp['hru_inundation'] = self.routing.hru_inundation[:]
   fp['hband_inundation'] = self.routing.hband_inundation[:]
   fp['reach2hband_inundation'] = self.routing.reach2hband_inundation[:]
  fp.close()
   
  #Close the LSM
  del self.noahmp
 
  #Close the routing module
  if self.routing_module == 'kinematic':
   del self.routing

  #Close the files
  self.input_fp.close()
  self.output_fp.close()

  return
