import ROOT

from xvfbwrapper import Xvfb

# with Xvfb() as xvfb:
# launch stuff inside virtual display here.
# It starts/stops around this code block.
vdisplay = Xvfb()
#vdisplay.start()
_file0 = ROOT.TFile("kk.root");
mygeo_man = _file0.Get("default");
myworld = mygeo_man.FindVolumeFast("world_volume");
myworld.Draw("ogl");
v = ROOT.gPad.GetViewer3D();
v.SavePicture("a_geometry.jpg");
#vdisplay.stop()
