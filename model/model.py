from mpi4py import MPI
import json
import fiona
import glob
import numpy as np
import os
import sys

def Read_Metadata_File(file):

 metadata = json.load(open(file))

 return metadata

def Run_HydroBlocks(metadata,edir,cid,rdir):

 import datetime
 from dateutil.relativedelta import relativedelta
 #import sys
 #print("%s/model" % hb)
 #sys.path.append("%s/model" % hb)
 import model.HydroBlocks as HydroBlocks

 #Read in the metadata file
 info = metadata
 info['ncores'] = 1
 info['cid'] = cid
 info['cdir'] = '%s/%s' % (edir,cid)
 info['Qobs_file'] = '%s/data/obs/obs.pck' % rdir
 info['routing_file'] = '%s/%s/octopy.pck' % (edir,cid)
 info['input_file'] = '%s/%s/input_file_routing.nc' % (edir,cid)
 info['output'] = {"dir":"%s/output_data/%s" % (edir,cid),
     "vars":info['output']['vars'],
     "routing_vars":info['output']['routing_vars']}
 info['restart'] = {"flag":info['restart']['flag'],
                    "dir":"%s/restart_data/%d" % (edir,cid)}
 os.system('mkdir -p %s' % (info['restart']['dir']))

 #Define idate and fdate
 idate = datetime.datetime(metadata['startdate']['year'],metadata['startdate']['month'],metadata['startdate']['day'],0)
 fdate = datetime.datetime(metadata['enddate']['year'],metadata['enddate']['month'],metadata['enddate']['day'],0) + datetime.timedelta(days=1)
 
 #Run the segments for the model
 restart_frequency = metadata['segment']['restart_frequency']
 sidate = idate
 sfdate = idate
 while sidate < fdate:
  #sfdate = sidate + relativedelta(years=metadata['segment']['years_per_segment'])
  if restart_frequency == 'daily':
   sfdate = sidate + relativedelta(days=1)
  if restart_frequency == 'monthly':
   sfdate = sidate + relativedelta(months=1)
  if restart_frequency == 'yearly':
   sfdate = sidate + relativedelta(years=1)
  if sfdate > fdate: sfdate = fdate
  #Set the parameters
  info['idate'] = sidate
  info['fdate'] = sfdate 
  info['MPI'] = MPI
  #Run the model
  #Initialize
  HB = HydroBlocks.initialize(info)
  #Run the model
  HB.run(info)
  #Finalize
  HB.finalize()
  #Update initial time step
  sidate = sfdate

 #Delete HB and library
 del HB
 del HydroBlocks

 return

def run(comm,metadata_file):

 #Get some general info
 metadata = Read_Metadata_File(metadata_file)
 rdir = metadata['rdir']
 edir = '%s/experiments/simulations/%s' % (rdir,metadata['experiment'])
 size = int(comm.Get_size())
 rank = int(comm.Get_rank())

 #Get cid configuration
 dfile = '%s/data/shp/domain.shp' % rdir
 fp = fiona.open(dfile,'r')
 cids = np.array(range(1,len(list(fp))+1))
 fp.close()

 #Iterate per catchment
 for cid in cids[rank::size]:
  #print('cid:',cid,flush=True)
  #Create output directory
  os.system('mkdir -p %s/output_data/%d' % (edir,cid))
  #os.chdir('%s' % (edir,))
  #mdfile = '%s/%s/metadata.json' % (edir,cid)
  Run_HydroBlocks(metadata,edir,cid,rdir)
