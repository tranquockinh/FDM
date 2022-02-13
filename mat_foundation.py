import numpy as np
from pandas import *
import matplotlib.pyplot as plt
# geometry

L = 10 # m
B = 1.0 # m
h = 0.3 # m
n_elems = 10

n_nodes = n_elems + 1
dx = L / n_elems
x = []

for j in range(n_nodes):
	x.append(j * dx)

# material

E_rc = 60E3 # MN/m2
I_rc = B * h**3 / 12
K_s = 200E3 # MN/m

# loads

P = np.array([[1, 6E03],
			  [5, 6E03],
			  [9, 6E03]])

M = np.array([[1, 0],
			  [3, 0]])

force_vector = np.zeros((len(x), 1))
moment_vector = np.zeros((len(x), 1))

def loading_full(x, force):
	output = np.zeros((len(x), 1))
	for i,j in zip(force[:,0], force[:,1]):
		output[int(i),0] = j
	return output

force_vector = loading_full(x, P)
moment_vector = loading_full(x, M)

# initialization

constant_matrix = np.zeros((n_nodes, n_nodes))
load_vector = np.zeros((n_nodes,1))
m = K_s * dx**2 * B
n = E_rc*I_rc / dx**2
k = K_s * dx * B

def matrix_show(weightMatrix, index ='row {}',columns='column {}',Name = 'Matrix'):

        lambdaLabels = []
        layerLabels = []
        for i in range(len(weightMatrix[:,0])):
            lambdaLabels.append(index.format(i+1))
        for j in range(len(weightMatrix[0,:])):
            layerLabels.append(columns.format(j+1))

        df = DataFrame(weightMatrix,index=Index(lambdaLabels, name=Name),columns=layerLabels)

        return df 

def constant_matrix_assemble(x):

	at_point = 1

	for i in range(1, n_elems):
		
		pakage_1 = np.zeros((n_nodes))
		pakage_2 = np.zeros((n_nodes))
		
		pakage_1[at_point-1] = 1 * n
		pakage_1[at_point] = -2 * n
		pakage_1[at_point+1] = 1 * n
		
		constant_matrix[i,:] = pakage_1

		at_point += 1

		temp = 0
		
		for j in range(i):

			temp = m * (i-j)
			
			constant_matrix[i,j] += temp

	# Boundary
	for j in range(n_nodes):
		
		constant_matrix[0,j] = k
		constant_matrix[-1,j] = j*dx * k

	constant_matrix_show = matrix_show(constant_matrix)
	print (constant_matrix_show)

	return constant_matrix_show


def loads_assemble(x):

	load_vector[0,0] = sum(P[:,1])
	load_vector[-1,0] = sum(M[:,1])

	for i,j in zip(P[:,0], P[:,1]):
		
		load_vector[-1,0] += j * x[int(i)]

	for i in range(1, n_elems):

		temp = 0
		
		for j in range(i):

			temp += force_vector[j,0] * (i-j)
			load_vector[i,0] = temp

	print (load_vector)

	return load_vector 


def solver(constant_matrix, load_vector):

	constant_matrix_show = constant_matrix_assemble(x)

	load_vector = loads_assemble(x)

	displacement = np.matmul(np.linalg.inv(constant_matrix), load_vector)

	subgound_reaction = displacement * K_s

	print (displacement)

	return displacement, subgound_reaction

displacement, subgound_reaction = solver(constant_matrix, load_vector)

def ploting(displacement, subgound_reaction):

	fig, ax = plt.subplots(2,1)
	ax[0].plot(x, displacement)
	ax[0].invert_yaxis()
	ax[0].set_ylabel('displacement, m')
	ax[1].plot(x, subgound_reaction)
	ax[1].set_xlabel('x')
	ax[1].set_ylabel('subgrade pressure, MN/m2')
	plt.show()

	return None

ploting(displacement, subgound_reaction)
