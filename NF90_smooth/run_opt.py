import test_membrane
import numpy as np
from scipy.optimize import minimize
x = test_membrane.Mem_Proteus()
res = minimize(x.proteus_eval,9.08e8, method='nelder-mead',options={'xtol': 1.e5,'disp': True})
res.x