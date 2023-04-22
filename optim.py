from utils import *
from scipy.optimize import minimize, Bounds

t_exp = np.loadtxt('data.txt')

#PCM_k, PCM_c, PCM_h, co2_h, pcm_h
xl = np.array([5.00, 0.01e3, 0.1, 0.,0.])
xu = np.array([20.00, 1.0e3, 2., 3.,5.0])
# xl = np.array([10.52, 0.05e3, 0.1, 1.])
# xu = np.array([30.52, 0.5e3, 2., 3.])
# xl = np.array([5.26, 0.5e3, 0.1, 1.])
# xu = np.array([10.52, 1.87e3, 2., 3.])
# xl = np.array([0.53, 1.87e3, 0.1, 1.])
# xu = np.array([5.26, 3.40e3, 2., 3.])
# n_var = 4 
# n_obj = 1
# n_constr = 0
x0 = np.array([10.2, 1.08e3, 1.1, 2.2, 2.2])

bounds = Bounds(xl,xu)

obj = lambda x: objfunc(x,t_exp)
# algorithm = GA(
#     pop_size=100,
#     eliminate_duplicates=True)

# res = minimize(obj, algorithm, seed=1, verbose=False)

# print("Best solution found: \nX = %s\nF = %s" % (res.X, res.F))
res = minimize(obj,x0,bounds=bounds,options={'disp':True}) 
print(res)