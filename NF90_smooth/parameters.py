den_dry_mem = 1.24*(1.e-7*1.e-7*1.e-7) # g/nm^3      1.24 #g/cm^3

C_brine_in = 1.753e-24
Rej = 0.70

Vol = 1.8e-5  #m^3/mol
dP = 670000 #Pa
R_gas = 8.314 # m^3-Pa/K-mol
Temp = 293.15 # K

D_ww = 1.36387090680776E13#2.00e18 #nm^2/s
K_pw = 0.125


den_wet_mem = 1.65*(1.e-7*1.e-7*1.e-7)  #g/nm^3
mass_poly = den_dry_mem / den_wet_mem # unitless