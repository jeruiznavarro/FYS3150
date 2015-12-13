from mdconfig import *
program = MD(dt=10, name="molecular-dynamics-fys3150") # dt in femtoseconds, name is the name of the executable you have
program.runNewSystem(numberOfUnitCells = 10, FCCLatticeConstant = 5.26) # lattice constant in angstroms
program.runThermostat(temperature = 300, timesteps = 100)
program.runThermalize(timesteps = 100)
program.runThermalize(timesteps = 10000)
program.saveState(path="states/myFavoriteState")
