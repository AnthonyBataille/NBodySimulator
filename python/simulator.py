#!/home/antho/anaconda3/bin/python
# -*- coding: utf-8 -*-


#from numba import jit
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#import cProfile

M_sun = 1.99e30
L_ref = 1.5e11
T_ref = 3600*24*365*10.0
fdpi = 100

def generate_solar_system(N, flength, G, M):
    point = np.empty([2, N], dtype = np.float)
    point[:, 0] = [0.0, 0.0]
    point[:, 1:] = np.random.uniform(size = (2, N - 1), low = - flength / 2, high = flength / 2)
    
    mass = np.ones([N], dtype = np.float)
    mass[0] = M
    
    vel = np.zeros([2, N])
    vel[0, 1:] = -point[1, 1:]
    vel[1, 1:] = point[0, 1:]
    r1 = np.sqrt(G*M/np.sqrt(np.sum(point[:, 1:]**2, axis = 0)))
    vel[:, 1:] = np.multiply(r1, vel[:, 1:]/np.sqrt(np.sum(vel[:, 1:]**2, axis = 0)))
    return point, mass, vel

class simulator():
    def init_data(self, mass_dist, vel_dist, point_dist):
        if mass_dist == 'ones':
            self.mass = np.ones((self.N), dtype = np.float)
        elif mass_dist == 'uniform':
            self.mass = np.uniform(size = (self.N), low = 1.0, high = 100.0)
            
        if vel_dist == 'zeros':
            self.vel = np.zeros((2, self.N), dtype = np.float)
        elif vel_dist == 'normal':
            self.vel = np.random.normal(size = (2, self.N), loc = 0.0, scale = 100.0)
            
        if mass_dist == 'uniform':
            self.vel = np.random.uniform(size = (2, self.N), low = 0.0, high = self.flength / 2)
            
        if point_dist == 'uniform':
            self.point = np.random.uniform(size = (2, self.N), low = -self.flength / 2, high = self.flength / 2)
            
    def __init__(self, preset, N, flength, fsize, mass_dist, vel_dist, point_dist):
        
        self.N = N
        
        self.flength = flength
        self.fsize = fsize
        self.K = (flength / (fdpi * fsize))
        
        self.M_ref = M_sun
        self.L_ref = L_ref
        self.T_ref = T_ref

        self.G = 6.674e-11 * self.M_ref * self.T_ref**2 / (self.L_ref)**3
        
        if(preset == 'default'):
            self.init_data(mass_dist, vel_dist, point_dist)
        elif(preset == 'solar'):
            self.point, self.mass, self.vel = generate_solar_system(N, flength, self.G, 1000.0)
        
        self.accel = self.update_accel(self.point, self.mass, self.N, self.G)
        
        self.fig = plt.figure(figsize = (fsize,fsize), dpi = fdpi)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlim(-flength / 2, flength / 2)
        self.ax.set_ylim(-flength / 2, flength / 2)
        self.sc = None
            
    #@jit(nopython=True)
    def _update_masks(self, N):
        del_elems = np.empty((2 * N), dtype = np.int)
        for i in range(0, 2 * N, 2):
            del_elems[i] = N * i + i / 2
            del_elems[i + 1] = del_elems[i] + N
            
        reorder_elems = np.empty((N * 2 * (N - 1)), dtype = np.int)
        for i in range(N):
            for j in range(N - 1):
                reorder_elems[i * 2 * (N - 1) + j] = (i // 2) * 2 * (N - 1) + (i % 2) * (N - 1) + j
                reorder_elems[i * 2 * (N - 1) + N - 1 + j] = ((N + i) // 2) * 2 * (N - 1) + ((N + i) % 2) * (N - 1) + j
        return (del_elems, reorder_elems)
    
    def update_accel(self, point, mass, N, G):
        elems = self._update_masks(N)
        del_elems = elems[0]
        reorder_elems = elems[1]
        
        p = np.array(np.broadcast_to(point, (N, 2, N)))
        p = np.delete(p.ravel(), del_elems, None).reshape([N, 2, N - 1])
        p2 = np.repeat(point, N - 1, 1).reshape([N, 2, N - 1])
        p2 = (p2.ravel())[reorder_elems].reshape([N, 2, N - 1])
        r = p - p2
        r2 = 1 / np.sum(r**2, axis = 1)
        r = np.multiply(np.multiply(r, np.repeat(np.sqrt(r2), 2, 0).reshape(N, 2, N - 1)), np.repeat(r2, 2, 0).reshape(N, 2, N - 1))
        m = np.delete(np.broadcast_to(mass, (N, N)).ravel(), np.arange(0, N**2, N + 1), None).reshape(N, N - 1, 1)
        #accel = G * np.multiply(mass.repeat(2, 0).reshape(N, 2), np.matmul(r, m).reshape(N, 2)) #force
        accel = G * np.matmul(r, m).reshape(N, 2)
    
        return accel.transpose()
    
    #@jit(nopython=True)
    def _remove_indexes(self, distances, mass, vel, N, L):
        indexes = np.ones((N), dtype = np.bool)
        for i in range(N):
            for j in range(i):
                if distances[i, j] < L[i, j] and indexes[i] != False and indexes[j] != False:
                    sum_mass = mass[i] + mass[j]
                    if mass[i] > mass[j]:
                        indexes[j] = False
                        vel[:, i] = mass[i] / sum_mass * vel[:, i] + mass[j] / sum_mass * vel[:, j] 
                        mass[i] += mass[j]
                    else:
                        indexes[i] = False
                        vel[:, j] = mass[i] / sum_mass * vel[:, i] + mass[j] / sum_mass * vel[:, j] 
                        mass[j] += mass[i]
                    N -= 1
        return indexes, mass, vel, N
    
    def merge_point(self, point, vel, accel, mass, N, K):
        distances_x = np.broadcast_to(point[0, :], [N, N]) - point[0, :].repeat(N, 0).reshape([N, N])
        distances_y = np.broadcast_to(point[1, :], [N, N]) - point[1, :].repeat(N, 0).reshape([N, N])
        L = K * (np.broadcast_to(mass, [N, N]) + mass.repeat(N, 0).reshape([N, N]))
        distances = distances_x**2 + distances_y**2
        ri, mass, vel, N2 = self._remove_indexes(distances, mass, vel, N, L)
        point = point[:, ri]
        mass = mass[ri]
        vel = vel[:, ri]
        accel = accel[:, ri]
        return point, vel, accel, mass, N2
    
    def update_point(self, point, vel, accel, dt):
        point += vel * dt
        vel += accel * dt
        return point, vel
    
    def init_figure(self):
        self.sc = self.ax.scatter(self.point[0, :], self.point[1, :], color = 'blue', s = 2 * self.mass)
        #print('init_figure')
        return self.sc,
    
    def update_figure(self, frame):
        #print(frame)
        self.point, self.vel, self.accel, self.mass, self.N = self.merge_point(self.point, self.vel, self.accel, self.mass, self.N, self.K)
        #cProfile.runctx('self.merge_point(self.point, self.vel, self.accel, self.mass, self.N, self.K)', None, locals())
        self.accel = self.update_accel(self.point, self.mass, self.N, self.G)
        #cProfile.runctx('self.update_accel(self.point, self.mass, self.N, self.G)', None, locals())
        self.point, self.vel = self.update_point(self.point, self.vel, self.accel, 0.001)
        
        #print((self.point[0, 1] - self.point[0, 0])**2 + (self.point[1, 1] - self.point[1, 0])**2)
        self.sc.set_offsets(self.point.transpose())
        self.sc.set_sizes(2 * self.mass)
        #hp = np.argmax(self.mass)
        
        #self.ax.set_xlim(self.point[0, hp] - self.flength / 2, self.point[0, hp] + self.flength / 2)
        #self.ax.set_ylim(self.point[1, hp] - self.flength / 2, self.point[1, hp] + self.flength / 2)
        
        return self.sc,
    
    def run(self):
        anim = animation.FuncAnimation(self.fig, self.update_figure, 100, init_func = self.init_figure, interval = 1000/60, blit = True, repeat = True)
        #raise Exception('Waiting..')
        plt.show()
        
#if __name__ == '__main__':
sim = simulator('solar', 2, 200, 10, 'ones', 'normal', 'uniform')
sim.run()
