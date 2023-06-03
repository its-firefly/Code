from utils import *
from scipy.optimize import minimize, Bounds

t_exp_vals = np.loadtxt('experimentalData.txt')
t_exp = t_exp_vals[0:1800]

#PCM_k, PCM_c, PCM_rho, PCM_h, co2_h, coeff
xl = np.array([10.00, 1.1e3, 10, 0.1, 0.,0.]) #lower bounds
xu = np.array([20.00, 5.5e3, 100, 2., 3.,20.0]) #upper bounds

x0 = np.array([10.2, 2.18e3, 100, 1.1, 2.2, 2.2]) #initial guess

bounds = Bounds(xl,xu)

obj = lambda x: objfunc(x,t_exp)

res = minimize(obj,x0,bounds=bounds,options={'disp':True}) 
print(res)


"""

Note that the optimized results that are obtained from the above process do not represent actual properties of PCM, instead, they are to be treated as assumptions for the next step of the analysis where the gas layer is replaced with other gases like krypton, xenon and air instead of CO_2

""" 



