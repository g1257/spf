#!/usr/bin/perl -w

use strict;

while(<STDIN>) {
	s/^OPTIONS +/EngineOptions=/;
	s/^VERSION +/Version=/;
	s/^FILENAME +/OutputFilename=/;
	s/^CARRIERS +/Carriers=/;
	s/^MU +/ChemicalPotential=/;
	s/^BETA +/Beta=/;
	s/^MCTHERMALIZATIONS +/MonteCarloThermalizations=/;
	s/^MCEFFECTIVE +/MonteCarloEffectiveIteractions=/;
	s/^MCUNMEASURED +/MonteCarloUnmeasuredPlusOne=/;
	s/^MCWINDOW/MonteCarloWindows/;
	s/^MCSTARTTYPE +/MonteCarloStartTypeOrFile=/;
	s/^MCSTARTLEVEL +/MonteCarloStartLevel=/;
	s/^HISTOGRAMSTEPS +/HistogramSteps=/;
	s/^BOUNDARYCONDITIONS +/BoundaryConditions=/;
	s/^SEED +/RandomSeed=/;
	s/^LINSIZE +/LatticeLength=/;
	s/^hoppings/Hoppings/;
	s/^PARAMETERSJ +/CouplingJ=/;
	s/^POTENTIALV/PotentialV/;
	s/^JAFNN +/JAFNN=/;
	s/^JAFNNN +/JAFNNN=/;
	s/^MagneticField +/MagneticField=/;
	s/^PARAMETERSJ_AF +/PARAMETERSJ_AF/;
	s/^PARAMETERSJ_AF_NN +/PARAMETERSJ_AF_NN=/;
	s/^SPIN_ORBIT_COUPLING +/SPIN_ORBIT_COUPLING=/;
	print;
}







