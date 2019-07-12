## Written by - Sheel Nidhan
## Last modified - 24th June 2019
## Extracts the axial slices from data files output of Eddy6 code
# Setting up the environments
import numpy as np
import struct as st

def readgrid(filename):

    f = open(filename)
    nx = int(f.read(4))

    index,x=np.loadtxt(filename,skiprows=1,unpack=1)

    xe = np.zeros(nx+1)
    xc = np.zeros(nx+1)

    xe[:len(x)]=x[:]; xe[-1]=xe[-2] + (xe[-2]-xe[-3]); #only required for consistency of the allocation

    for ii in range(1, nx+1):

        xc[ii] = (xe[ii] + xe[ii-1]) * 0.5 

    xc[0] = xe[0] - (xc[1]-xe[0]) 

    return  nx, index, x, xe, xc  

def readres(filename):

    f = open(filename,'rb')

    _,i,j,k,jp,_ = st.unpack('i'*6,f.read(24))

    data=np.zeros([i,j,k],dtype=float)

    for kk in range (0,k):

        dummy1 = st.unpack('i',f.read(4))
        data[:,:,kk]=np.reshape(st.unpack('d'*i*j,f.read(8*i*j)),(i,j),order='F')
        dummy2 = st.unpack('i',f.read(4))

        if dummy1!=dummy2:
            print('Error reading',kk)
            break

    _,it,_,_,time,_,_,dt,grav,_=st.unpack('iiiidiiddi',f.read(4*7+3*8))


    return i,j,k,it,time,dt,grav,data

# Read 5p data file
def read5p2(filename):
        f = open(filename,'rb')

        _,i,j,k,jp,_ = st.unpack('i'*6,f.read(24))
        
        k = 220
        i = i - 10
        print('i',i)
        print('j',j)
        
        data=np.zeros([i,j,k],dtype=float)

        for kk in range (0,k):
    
                dummy1 = st.unpack('i',f.read(4))    
                data[:,:,kk]=np.reshape(st.unpack('d'*i*j,f.read(8*i*j)),(i,j),order='F')
                dummy2 = st.unpack('i',f.read(4))
    
                if dummy1 != dummy2:
                        print('Error reading',kk)
                        break
        
        
        _,it,_,_,time,_,_,dt,grav,_=st.unpack('iiiidiiddi',f.read(4*7+3*8))
        
        return i,j,k,it,time,dt,grav,data


# Parameters for the 5 plane data files
kmin = 298
kmax = 302

NKRR = np.int(kmax - kmin + 1 + 5*np.floor((4610-kmax)/100))
print('NKRR', NKRR) 

# Read Grid
nx, index, x, xe, xc = readgrid('./x3_grid.in')
nr, index, r, re, rc = readgrid('./x1_grid.in')
print("Number of grid points in x direction", nx)
print("Number of grid points in r direction", nr)


# Map the grid to actual location of the data files
kmin = 298
kmax = 302

NKRR = np.int(kmax - kmin + 1 + 5*np.floor((4610-kmax)/100))
print('NKRR', NKRR) 

xc_data_file = np.zeros([NKRR],dtype=float)
xc_data_file[0:kmax-kmin+1] = xc[kmin-1:kmax] # Storing the disk location

mk = (kmax-kmin+1)   
nk = 1

for i in range(0,np.int(np.floor((4610-kmax)/100))):
    i_actual = kmax + nk*100 - 1
    xc_data_file[mk]   = xc[i_actual-2]  
    xc_data_file[mk+1] = xc[i_actual-1]  
    xc_data_file[mk+2] = xc[i_actual]  
    xc_data_file[mk+3] = xc[i_actual+1] 
    xc_data_file[mk+4] = xc[i_actual+2]  
    mk = mk + 5
    nk = nk + 1

rc_data_file = rc[0:nr-10+1]



# Extracting array location at a given x/D
v = np.array([5, 10, 15, 60, 65, 70, 100])
idx = np.zeros((7))
size_v = np.shape(v)
for i in range(0, size_v[0]):
    idx[i] =  np.int((np.abs(xc_data_file - v[i])).argmin())

idx = idx.astype(int)

print('Spatial locations of extraction ', xc_data_file[idx]) 

nstart = 1892600
nend   = 2000000
stride = 100
ntheta = 258
tot_size = np.int((nend-nstart)/stride + 1)

path_up = '/p/work/snidhan/spod_re5e4/frinf/DATA_FILES/up_datafiles/'
path_vp = '/p/work/snidhan/spod_re5e4/frinf/DATA_FILES/vp_datafiles/'
path_wp = '/p/work/snidhan/spod_re5e4/frinf/DATA_FILES/wp_datafiles/'
path_up_slice = '/p/work/snidhan/spod_re5e4/frinf/DATA_FILES/slices/up_slice/'
path_vp_slice = '/p/work/snidhan/spod_re5e4/frinf/DATA_FILES/slices/vp_slice/'
path_wp_slice = '/p/work/snidhan/spod_re5e4/frinf/DATA_FILES/slices/wp_slice/'

file_up = 'up_'
file_vp = 'vp_'
file_wp = 'wp_'

count = 0

file = 'time_stamp.txt'
time_filename = './' + file
    
for step in range(nstart,nend+100,100):
    
    print('Writing ', step)
    
    i_pad = str(step).zfill(8)
    
    
    filename = path_up + file_up +  i_pad + '.res'
    print('Reading up file ', filename)
    i,j,k,it,time,dt,grav,data = read5p2(filename)
    
    for i in range(0,size_v[0]):
        filename1 = path_up_slice + file_up + i_pad + '_'+ str(v[i]) + '.res'
        print(filename1)
        data_f = data[:,:,idx[i]].reshape(nr-10+1,ntheta,order='F')
        data_f.T.tofile(filename1)

    filename = path_vp + file_vp +  i_pad + '.res'
    print('Reading vp file ', filename)
    i,j,k,it,time,dt,grav,data = read5p2(filename)
    
    for i in range(0,size_v[0]):
        filename1 = path_vp_slice + file_vp + i_pad + '_'+ str(v[i]) + '.res'
        print(filename1)
        data_f = data[:,:,idx[i]].reshape(nr-10+1,ntheta,order='F')
        data_f.T.tofile(filename1)

    filename = path_wp + file_wp +  i_pad + '.res'
    print('Reading wp file ', filename)
    i,j,k,it,time,dt,grav,data = read5p2(filename)
    
    for i in range(0,size_v[0]):
        filename1 = path_wp_slice + file_wp + i_pad + '_'+ str(v[i]) + '.res'
        print(filename1)
        data_f = data[:,:,idx[i]].reshape(nr-10+1,ntheta,order='F')
        data_f.T.tofile(filename1)
        
    print('Written ', step)
    print('Time ', time)
    with open(time_filename, 'a') as f:
	f.write(str(step) + "    ")
	f.write(str(time))
	f.write('\n')
    f.close()
	
    count = count + 1
