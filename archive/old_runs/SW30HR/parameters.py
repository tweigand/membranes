den_dry_mem = 1.31*(1.e-7*1.e-7*1.e-7) # g/nm^3      1.24 #g/cm^3

C_brine_in = 1.753e-24
Rej = 0.70

Vol = 1.8e-5  #m^3/mol
dP = 670000 #Pa
R_gas = 8.314 # m^3-Pa/K-mol
Temp = 293.15 # K

D_ww = 2.00e18 #nm^2/s
K_pw = 0.19 


den_wet_mem = den_dry_mem + K_pw*1.e-21
mass_poly = den_dry_mem / den_wet_mem # unitless



def den_w(b):
  return 1.e-21 #0.2631*b*b + 0.6947*b + 0.9985

def d_den_w(b):
  return 0.#2.0*0.2631*b + 0.6947

def den_s(w,b):
  return den_dry_mem/(1.-w-b)

def d_den_s(w,b):
  return den_dry_mem/((b+w-1.)*(b+w-1.))