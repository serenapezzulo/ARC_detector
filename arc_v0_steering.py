# ddsim steering file steering.py
# to run it: ddsim --steeringFile arc_v0_steering.py -N 10 --outputFile myout.root
from DDSim.DD4hepSimulation import DD4hepSimulation
SIM = DD4hepSimulation()

## The compact XML file, or multiple compact files, if the last one is the closer.
SIM.compactFile = ['compact/arc_v0.xml']

# Register optical physics and Cerenkov process to the default physics
def setupCerenkov(kernel):
        from DDG4 import PhysicsList

        seq = kernel.physicsList()
        cerenkov = PhysicsList(kernel, "Geant4CerenkovPhysics/CerenkovPhys")
        cerenkov.MaxNumPhotonsPerStep = 10
        cerenkov.MaxBetaChangePerStep = 10.0
        cerenkov.TrackSecondariesFirst = False
        cerenkov.VerboseLevel = 0
        cerenkov.enableUI()
        seq.adopt(cerenkov)
        ph = PhysicsList(kernel, "Geant4OpticalPhotonPhysics/OpticalGammaPhys")
        ph.addParticleConstructor("G4OpticalPhoton")
        ph.VerboseLevel = 0
        ph.BoundaryInvokeSD = True
        ph.enableUI()
        seq.adopt(ph)
        return None
SIM.physics.setupUserPhysics(setupCerenkov)

# Associate the Geant4OpticalTrackerAction to these detectors
# this action register total energy of the opt photon as a single hit
# and kills the optical photon, so no time is wasted tracking them
SIM.action.mapActions["ARC_DETECTORNAME"] = "Geant4OpticalTrackerAction"


# Register hit with low energy, compatible with zero
SIM.filter.tracker = "edep0"

# Define filter, so detector is only sensitive to optical photons
# In real life charged particles deposit energy as well!!
SIM.filter.filters["opticalphotons"] = dict(
        name="ParticleSelectFilter/OpticalPhotonSelector",
        parameter={"particle": "opticalphoton"},
        )
SIM.filter.mapDetFilter["ARC_DETECTORNAME"] = "opticalphotons"


# Particle gun settings: pions with fixed energy, random direction
SIM.numberOfEvents = 1000
SIM.enableGun = True
SIM.gun.energy = "50*GeV"
SIM.gun.particle = "pi+"
SIM.gun.distribution = "uniform"
SIM.gun.multiplicity = 1
SIM.gun.position = "0 0 0"


# without visualization
SIM.runType = "batch"

#with visualization
# SIM.runType = "qt"
# SIM.macroFile = "vis.mac"

