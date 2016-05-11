import numpy as np
import random as rd
import math

pi = 3.14159

f = open('Domain.poly', 'w')
g = open('pores.data', 'w')

xs = np.zeros(500,dtype="float")
ys = np.zeros(500,dtype="float")
rs = np.zeros(500,dtype="float")


ly = 131.3
poro_mem = 0.15
pore_ave = 30.4
pore_stdev = 7.9


rad = 12
lx = 5000.
nx_ele = 750
dx = lx/nx_ele
ny_ele = 105
dy = ly/ny_ele

V_tot = lx*ly
V_pores = 0.
poro = V_pores/V_tot


count_b = 0
while (poro < poro_mem):
	x = rd.uniform(0,lx)
	y = rd.uniform(0,ly)
	r = rd.gauss(pore_ave,pore_stdev)
	r = r/2.
	print poro
	check = True
	if (count_b == 0 and x-r>0.0 and x+r<lx and y-r>0.0 and y+r<ly):
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
		print poro_check,poro_mem
		if ( x-r>0.0 and x+r<lx and y-r>0.0 and y+r<ly and check==True and poro_check < poro_mem + 0.0001):
			xs[count_b] = x	
			ys[count_b] = y
			rs[count_b] = r
			g.write("%f %f %f \n" %(x,y,r))
			count_b = count_b + 1
			V_pores = V_pores + pi*r*r
			poro = V_pores/V_tot

print poro

nb_nodes = nx_ele*2 + ny_ele*2 
n_nodes = nb_nodes + count_b*rad

f.write("%i 2 0 1\n" %(n_nodes))


# For Boundary Nodes #
# 1 = No Flux Side Boundaries
# 2 = Bottom Boundary
# 3 = Top Boundary


# Bottom Boundary 
count = 1
for i in range(nx_ele+1):
	f.write("%i %f %f 2 \n" %(count,dx*i,0))
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
for i in range(nx_ele):
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