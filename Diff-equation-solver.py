import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

def odeFunc(t, R):

	r2 = (R[[0]]**2 + R[[1]]**2 + R[[2]]**2)**(0.5)
        
    # Constants
	mu = 3.98716708E5 # Gravitational parameters [km^3/s^2]
	Re = 6378.137     # Radius of the earth
	J = 0.0010826267
	gamma = (5*R[2]**2)/(r2**2) # Partial coefficients used to calculate ap
	lambd = (-3*J*mu*Re**2)/(2*r2**5)

    # Pertubation due to non-spherical earth
	ap = np.reshape(np.array([[lambd*R[0]*(1-gamma)],[lambd*R[1]*(1-gamma)], [lambd*R[2]*(3-gamma)]]), (3,1))

	r_ddot = (-mu*R[0:3])/(r2**3) + ap # Acceleration

	ddy = [np.zeros((6,1))]
	
	ddy = np.reshape(np.array([[R[3]], [R[4]], [R[5]], [r_ddot[0]], [r_ddot[1]], [r_ddot[2]]]), (6,1)) # Returned values for velocity and acceleration [[dr/dt]; [d^2r/dt^2]]

	return ddy

def rungeKutta(func, tspan, steps, y0, order):

	[m, n] = np.shape(y0)
	dy = np.zeros((m, (steps + 1)))

	t  = np.zeros((1, (steps + 1)))

	dy[:,[0]] = y  = y0
	t[0]  = ti = tspan[0] 
	h = (tspan[1] - tspan[0]) / float(steps)
	
	if order not in [1, 2, 4]:

		print('Error! Order integer must be == 1, 2 or 4!')

		return

	# FIRST ORDER RUNGE-KUTTA
	for i in xrange(1, steps + 1): # ITERATE 
		
		if order == 1:

			k1 = odeFunc(ti, y)

			t[0,i]  = ti = tspan[0] + i * h
			dy[:,[i]] = y  = y + k1 * h
			
	## SECOND ORDER RUNGE-KUTTA
		if order == 2:

			k1 = odeFunc(ti, y) * h
			k2 = odeFunc(ti, y + k1 / 2) * h

			t[0,i]  = ti = tspan[0] + i * h
			dy[:,[i]] = y  = y + k2

	# FOURTH ORDER RUNGE-KUTTA
		if order == 4:

			k1 = odeFunc(ti, y) * h
			k2 = odeFunc(ti + h / 2, y + k1 / 2) * h
			k3 = odeFunc(ti + h / 2, y + k2 / 2) * h
			k4 = odeFunc(ti + h, y + k3) * h

			t[0,i]  = ti = tspan[0] + i * h
			dy[:,[i]] = y  = y + (k1 + 2 * (k2 + k3) + k4) / 6

	return(dy, t)


/2), fargs=(data, lines),
                                   interval=t_steps, blit=False)


plt.show()
