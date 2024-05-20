# This program solves steady state Poisson equation explicitly
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
data = np.loadtxt('D:\PYTHON\LAPLACE_SOLUTION/INPUT_DATA.txt')
print(data)
print('\n Assigning: \n n_rows: number of rows', '\n n_columns: number of columns', \
      '\n step_size=step size in x and y direction')
n_rows = int(data[0])
n_columns = int(data[1])
step_size = float(data[2])
max_iter = int(data[3])

print('\n n_row =', n_rows, '\n n_columns =', n_columns, '\n step_size =', step_size, '\n Max iterations=', max_iter)

print('\n Assigning hi: Initial Hydraulic Head')
hi = np.loadtxt('D:\PYTHON\LAPLACE_SOLUTION/hi.txt', dtype='f', delimiter='\t')
print(type(hi))
print(hi)

print('\n Assigning kx: Hydraulic Conductivity along x axis')
kx = np.loadtxt('D:\PYTHON\LAPLACE_SOLUTION/kx.txt', dtype='f', delimiter='\t')

print('\n Assigning ky: Hydraulic Conductivity along y axis')
ky = np.loadtxt('D:\PYTHON\LAPLACE_SOLUTION/ky.txt', dtype='f', delimiter='\t')

print('\n Assigning w: net abstraction')
w = np.loadtxt('D:\PYTHON\LAPLACE_SOLUTION/w.txt', dtype='f', delimiter='\t')

h=hi

for l in range(1,max_iter+1):
      print(l,h[5,5])
      for i in range(1,n_rows-1):
            for j in range(1,n_columns-1):
                  h[i,j]=(kx[i,j-1]*h[i,j-1]+kx[i,j+1]*h[i,j+1]+ky[i-1,j]*h[i-1,j]+ky[i+1,j]*h[i+1,j] \
                          -w[i,j]*step_size*step_size)/(2*kx[i,j]+2*ky[i,j])


print(h)

np.savetxt('D:\PYTHON\LAPLACE_SOLUTION/FINAL_HEAD.txt',h,fmt='% 3.1f')

print(np.shape(h))

print('Drawing plot')
colorinterpolation = 50
colourMap = plt.cm.jet

# Set meshgrid
X, Y = np.meshgrid(np.arange(0, n_columns), np.arange(0, n_rows))
print(np.shape(X),np.shape(Y))
# Configure the contour
plt.title("Head Distribution (in m)")
plt.contourf(X, Y, h, colorinterpolation, cmap=colourMap)

# Set Colorbar
plt.colorbar()

# Show the result in the plot window
plt.show()
