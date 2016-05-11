import numpy as np
import random as rd
import math

pi = 3.14159

f = open('Domain.poly', 'w')
g = open('pores.data', 'w')

xs = np.zeros(500,dtype="float")
ys = np.zeros(500,dtype="float")
rs = np.zeros(500,dtype="float")


ly = 139.1
poro_mem = 0.15
pore_ave = 29.38
pore_stdev = 23.26

d_circ_min =  9.0
r_min = d_circ_min/2.0


rad = 12
lx = 5000.
nx_ele = 5000
dx = lx/nx_ele
ny_ele = 105
dy = ly/ny_ele

V_tot = lx*ly
V_pores = 0.
poro = V_pores/V_tot


amp_mean = 81.0
amp_stdev = 31.0
amp = rd.gauss(amp_mean,amp_stdev)
freq = round(2*np.pi/(8*pore_ave),2)

period = 1./freq

#print period

y_bottom_i = np.zeros(nx_ele+1,dtype="float")
y_bottom = np.zeros(nx_ele+1,dtype="float")
x_bottom = np.zeros(nx_ele+1,dtype="float")
for i in range(nx_ele+1):
	x_bottom[i] = i*dx
	y_bottom_i[i] = amp*np.sin(freq*i*dx)
	#print x_bottom[i],y_bottom_i[i],round(y_bottom_i[i])
	if (y_bottom_i[i] > 0.0 and y_bottom_i[i-1] < 0.0):
		amp = rd.gauss(amp_mean,amp_stdev)
		while (amp > ly - 10.0):
			amp = rd.gauss(amp_mean,amp_stdev)
		#print i,amp
		y_bottom_i[i] = amp*np.sin(freq*i*dx) 
		


# Smooth
"""y_bottom[0] = 0.0
y_bottom[1] = 0.5*y_bottom_i[1] + 0.25*y_bottom_i[2]
for i in range(2,nx_ele-2):
	y_bottom[i] = 0.25*y_bottom_i[i-1] + 0.5*y_bottom_i[i] + 0.5*y_bottom_i[i+1]

y_bottom[nx_ele-2] = 0.25*y_bottom_i[nx_ele-3] + 0.5*y_bottom_i[nx_ele-2]
y_bottom[nx_ele-1] = 0.0"""

y_bottom = y_bottom_i



y_b_min = min(y_bottom)

count_b = 0
while (poro < poro_mem):
	x = rd.uniform(0,lx)
	y = rd.uniform(y_b_min,ly)
	d_c = rd.gauss(pore_ave,pore_stdev)
	r = d_c/2.
	while (r < r_min):
		d_c = rd.gauss(pore_ave,pore_stdev)
		r = d_c/2.
	check = True

	### Check Bottom Boundary ###
	bottom_check  = True
	buffer = 5

	x_find = x_bottom[np.where(  (x_bottom > x - r)    &   (x_bottom < x + r)   )]
	y_find = y_bottom[np.where(  (x_bottom > x - r)    &   (x_bottom < x + r)   )]
	if (x_find.size == 0):
		bottom_check = False
	for i in range(0,x_find.size):
		if (y_find[i] > y):
			bottom_check = False
		d = np.sqrt( pow(x_find[i]-x,2) + pow(y_find[i] - y,2) )
		if (d < (r + buffer) ):
			bottom_check = False
	

	if (count_b == 0 and x-r>0.0 and x+r<lx and y+r<ly and bottom_check == True):
		xs[count_b] = x	
		ys[count_b] = y
		rs[count_b] = r
		count_b = count_b + 1
		V_pores = V_pores + pi*r*r
		poro = V_pores/V_tot
		g.write("%f %f %f \n" %(x,y,r))
	if (count_b > 0):
		for i in range(count_b):
			d = pow((x-xs[i])*(x-xs[i])+(y-ys[i])*(y-ys[i]),0.5)
			if (d<r+rs[i]):
				check = False
		V_pores_check = V_pores + pi*r*r
		poro_check = V_pores_check/V_tot
		#print poro_check,poro_mem
		if ( x-r>0.0 and x+r<lx and y+r<ly and check==True and poro_check < poro_mem + 0.0001 and bottom_check == True):
			xs[count_b] = x	
			ys[count_b] = y
			rs[count_b] = r
			g.write("%f %f %f \n" %(x,y,r))
			count_b = count_b + 1
			V_pores = V_pores + pi*r*r
			poro = V_pores/V_tot

print poro

nb_nodes = nx_ele*2 + ny_ele*2 - 1
n_nodes = nb_nodes + count_b*rad

f.write("%i 2 0 1\n" %(n_nodes))


# For Boundary Nodes #
# 1 = No Flux Side Boundaries
# 2 = Bottom Boundary
# 3 = Top Boundary


# Bottom Boundary 
count = 1
for i in range(nx_ele):
	f.write("%i %f %f 2 \n" %(count,dx*i,y_bottom[i]))
	count = count + 1

# Right Boundary 
for i in range(1,ny_ele):
	f.write("%i %f %f 1 \n" %(count,lx,dy*i))
	count = count + 1

# Top Boundary 
for i in reversed(range(nx_ele+1)):
	f.write("%i %f %f 3 \n" %(count,dx*i,ly))
	count = count + 1

# Left Boundary 
for i in reversed(range(1,ny_ele)):
	f.write("%i %f %f 1 \n" %(count,0,dy*i))
	count = count + 1

count = nb_nodes

for i in range(count_b):
	for j in range(rad):
		count = count + 1
		angle = j*(360./rad)*pi/180.
		xp = xs[i]+rs[i]*math.cos(angle)
		yp = ys[i]+rs[i]*math.sin(angle)
		f.write("%i %f %f 100 \n" %(count,xp,yp))


f.write("%i 1 \n" %(n_nodes))

count = 1 

# Bottom Boundary 
for i in range(nx_ele-1):
	f.write("%i %i %i 2 \n" %(count,count,count+1))
	count = count + 1

# Right Boundary 
for i in range(ny_ele):
	f.write("%i %i %i 1 \n" %(count,count,count+1))
	count = count + 1

# Top Boundary 
for i in range(nx_ele):
	f.write("%i %i %i 3 \n" %(count,count,count+1))
	count = count + 1

# Left Boundary 
for i in range(ny_ele):
	if (count == nb_nodes):
		f.write("%i %i 1 1 \n" %(count,count))
	else:
		f.write("%i %i %i 1 \n" %(count,count,count+1))
	count = count + 1

count = nb_nodes + 1
for i in range(count_b):
	count2 = count
	for j in range(rad-1):	
		f.write("%i %i %i 100 \n" %(count2+j,count2+j,count2+j+1))
	f.write("%i %i %i 100 \n" %(count2+j+1,count2+j+1,count))
	count = count2+j+1
	count = count + 1

f.write("0 \n")

rad = 12

f.write("%i \n"% (count_b))


count = 1
for i in range(count_b):
	xp = xs[i]
	yp = ys[i]
	rp = rs[i]
	f.write("%i %f %f 100 100 \n" %(count,xp,yp))
	count = count + 1


print "Mean of Y_bottom",np.mean(y_bottom_i)
