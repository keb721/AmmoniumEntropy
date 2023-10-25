# AmmoniumEntropy
Codes for the Blow, Whale, Quigley and Sosso Faraday Discussion "Understanding the impact of ammonium ion substitutions on heterogeneous ice nucleation"

"Makefile" produces "make_configs" and "make_gro"
"make_configs" generates hydrogen bond configurations when an input of water positions (3 atom, layering in $z$ direction - e.g. physical_water ice XI in genice) is used - "make_ice_base.py" may be required to refine this geometry
"make_gro" generates GROMACS input files using these hashes and the TIP4P/Ice geometry
"energy_spiderweb.py" produces spiderweb plots as seem in Fig.4 
"mismatch.f90" is used to identify large rotations of minimised water molecules indicating mechanical instabilities
