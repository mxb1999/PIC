import numpy as np
import math

from numpy.core.arrayprint import _make_options_dict

sigma = 1.7e-4;
e0 =8.85418782e-12;
me =9.10938356e-31;
pi =3.14159265359;
kb= 1.3806485279e-16;   
kb2= 1.3806485279e-23;   
ec= 1.60217662e-19;
c= 299792458.0;              
estat=4.80320427e-10; 	     
Z = 3.1;                        
Te = 2.0e3*11604.5052;          
Te_eV = 2.0e3;
Ti = 1.0e3*11604.5052;
mi_kg = 10230.0*me
def discretize(position: float,
               minx: float,
               maxx: float,
               nx: int):
    return int((position - minx)/(maxx-minx)*nx)

def push_1d(E: np.ndarray,
         B: np.ndarray, 
         positions: np.ndarray, 
         momenta: np.ndarray, 
         dt: float,
         dims: np.ndarray):
    num_particles = len(positions)
    minx = dims[0]
    maxx = dims[1]
    nx = num_particles
    q = ec*Z
    for i in range(num_particles):
        print(positions[i][0])
        ploc = momenta[i]
        positions[i][0] += ploc[0]*dt
        print(positions[i][0])
        print()
        print(ploc)
        local_E = E[i]
        local_B = B[i]
        epsilon_dt = q/2*local_E*dt
        p_minus = ploc + epsilon_dt
        speed = np.linalg.norm(p_minus/(c*mi_kg))
        gamma = 1/math.sqrt(1-(speed)**2)
        bmag = np.linalg.norm(local_B)
        theta = q*dt/gamma*bmag
        t = math.tan(theta/2)*local_B/bmag
        p_prime = p_minus + np.cross(p_minus, t)
        tmag = math.tan(theta/2)
        p_plus = p_minus + 2/(1+tmag**2)*np.cross(p_prime, t)
        momenta[i] = p_plus + epsilon_dt
        print(momenta[i])
        print('\n')


def initialize_one_ppc(positions, momenta, num_zones, minx, maxx):
    dx = (maxx-minx)/num_zones
    center = dx/2
    for i in range(num_zones):
        lower_x = dx*i
        positions[i][0] = lower_x+center
       # print((np.random.rand(3)*c*mi_kg)/(mi_kg*c))
        momenta[i] = np.random.rand(3)*c/2*mi_kg



def initialize(num_particles: int,
            mesh_shape: np.ndarray,
            dims: np.ndarray):
    E = np.zeros((*mesh_shape,3))
    B = np.zeros((*mesh_shape,3))
    for i in range(mesh_shape[0]):
        E[i] = np.array([1,0,0])
        B[i] = np.array([0,0,1])
    positions = np.zeros((num_particles, 3))
    momenta = np.zeros((num_particles, 3))
    dt = 1e-5
    initialize_one_ppc(positions, momenta, num_particles, dims[0], dims[1])
    push_1d(E,B,positions,momenta,dt,dims)

if __name__ == '__main__':
    initialize(10, np.array((10,)), np.array([-5e-4, 5e-4]))