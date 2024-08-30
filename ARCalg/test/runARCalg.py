#
# gaudi steering file that runs ARCalg
#
# to execute:
# k4run runARCalg.py

from Gaudi.Configuration import INFO,DEBUG
from Configurables import EventDataSvc
from k4FWCore import ApplicationMgr, IOSvc

svc = IOSvc("IOSvc")
svc.input = [ "arcsim_K+_30GeV.root"]
svc.output = "dummyout.root"

from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = ['/afs/cern.ch/user/s/spezzulo/private/ARC_detector/compact/arc_v0.xml']

from Configurables import ARCalg
ARCalg = ARCalg("ARCalg")
ARCalg.ARC_simhits=["ARC_HITS"]
ARCalg.MCParticles =["MCParticles"]
ARCalg.ARC_name="ARC_DETECTORNAME"

ARCalg.OutputLevel=INFO

mgr = ApplicationMgr(
    TopAlg=[ARCalg],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[geoservice,EventDataSvc("EventDataSvc")],
    OutputLevel=INFO,
)
