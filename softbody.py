
import pygame
import numpy as np

X = 1024
Y = 720
pygame.init()
win = pygame.display.set_mode((X, Y))


delaytime = 10  
timecorrect = 100
white = (255, 255, 255)
black = (0, 0, 0)
#radius of individual physics object
radius = 5

# Physics parameters
equilibrium = 15  # Distance between points
damp = 1        # Damping factor
gravity = np.array([0.0, 1])  # Gravity vector
e = 0.7  # Coefficient of restitution (for collision with boundaries)




def magnitude(vector):
    return np.linalg.norm(vector)

class PhysObj:
    def __init__(self, mass, pos, vel, color):
        self.mass = mass
        self.pos = np.array(pos, dtype=float)
        self.vel = np.array(vel, dtype=float)
        self.color = color

    def add_force(self, accel):
        self.vel += accel * (delaytime / timecorrect)

    def move(self):
        self.pos += self.vel * (delaytime / timecorrect)
        pygame.draw.circle(win, self.color, self.pos.astype(int), 5, 0)

    #Applies collision to the boundaries
    def boundary_check(self):
        

        if self.pos[0] > X - 2*radius:
            self.vel[0] *= -e
            self.pos[0] = X - 2*radius
        if self.pos[0] < 2*radius:
            self.vel[0] *= -e
            self.pos[0] = 2*radius
        if self.pos[1] > Y - 2*radius:
            self.vel[1] *= -e
            self.pos[1] = Y - 2*radius
        if self.pos[1] < 2*radius:
            self.vel[1] *= -e
            self.pos[1] = 2*radius

        if self.pos[0] >= X-2*radius:
            self.add_force(-(1+e)*timecorrect*self.vel[0]*np.array([1.0,0.0])/delaytime)           
            
        elif self.pos[0] <= 2*radius:
            self.add_force((1+e)*timecorrect*self.vel[0]*np.array([-1.0,0.0])/delaytime)
                       
        elif self.pos[1] >= Y-2*radius :
            self.add_force(-(1+e)*timecorrect*self.vel[1]*np.array([0.0,1.0])/delaytime)
                  
        elif self.pos[1] <= 2*radius:
            self.add_force((1+e)*timecorrect*self.vel[1]*np.array([0.0,-1.0])/delaytime)

#A class for a rectangular grid that behaves like a softbody
class Softbody:
    def __init__(self, pos, vel, omega, length, width, color, elast):
        self.points = []
        self.pos = pos
        self.vel = vel
        self.omega = omega
        self.length = length
        self.width = width
        for i in range(length):
            row = []
            for j in range(width):
                row.append(PhysObj(radius, pos + np.array([i * equilibrium, j * equilibrium]), vel, color))
            self.points.append(row)
            
        self.elast = elast
        self.bounds = []
        for i in range(length):
            for j in range(width):
                if i == 0 or j == 0 or i == len(self.points) - 1 or j == len(self.points[0]) - 1:
                    self.bounds.append(self.points[i][j])

    def rotate(self,angle):
        rot_matrix = np.array([[np.cos(angle),np.sin(angle)],[-np.sin(angle),np.cos(angle)]])
        for i in range(self.length):
            for j in range(self.width):
                rpos = self.points[i][j].pos - self.pos
                rpos_new = rot_matrix.dot(rpos)
                self.points[i][j].pos = rpos_new + self.pos

        if self.omega != 0:
            center = (self.points[0][0].pos + self.points[-1][-1].pos)/2

            for i in range(self.length):
                for j in range(self.width):
                    point = self.points[i][j]
                    rad = point.pos - center

                    point.vel += self.omega* np.array([-1*rad[1],rad[0]])


    #This function allows the softbody to appear and as a softbody, without this it wont appear in the screen
    def action(self):

        # Apply gravity and movement
        for row in self.points:
            for p in row:
                p.add_force(gravity)
                p.boundary_check()
                p.move()

        # Apply spring forces only to direct neighbors
        for i in range(self.length):
            for j in range(self.width):
                point = self.points[i][j]
                # Connect to right neighbor
                if i < self.length - 1:
                    self.apply_spring(point, self.points[i + 1][j])

                # Connect to bottom neighbor
                if j < self.width - 1:
                    self.apply_spring(point, self.points[i][j + 1])

                # Diagonal connections for better stability
                if i < self.length - 1 and j < self.width - 1:
                    self.apply_spring(point, self.points[i + 1][j + 1],np.sqrt(2))
                if i < self.length - 1 and j > 0:
                    self.apply_spring(point, self.points[i + 1][j - 1],np.sqrt(2))


    def apply_spring(self, p1, p2, factor = 1):
        #Uses simple spring force (F  = -kx)
        rad = p2.pos - p1.pos
        rad_mag = magnitude(rad)
        rad_norm = rad / rad_mag
        if rad_mag != 0:
            
            spring_force = self.elast * (rad_mag - factor*equilibrium)
            damping_force = -damp * rad_norm.dot(p1.vel - p2.vel)
            force = (spring_force + damping_force) * rad_norm

            p1.add_force(force / p1.mass)
            p2.add_force(-force / p2.mass)
        elif rad_mag < 0.5:
            p1.add_force(2*rad_norm)
            p2.add_force(-2*rad_norm)

#Applies collision between softbodies
def collide(obj1,obj2,colfactor = 1):
    for i in obj1.bounds:
        for j in obj2.bounds:
            rad = i.pos - j.pos
            rad_mag = magnitude(rad)
            rad_norm = rad / rad_mag

            if rad_mag < (equilibrium/5)*radius:
                i.add_force(colfactor*timecorrect*magnitude(i.vel - j.vel)*rad_norm/delaytime)
                j.add_force(colfactor*timecorrect*-1*magnitude(j.vel - i.vel)*rad_norm/delaytime)


#Creating Softbodies, it takes the position of the upper left corner point (x,y), velocity of the softbody (x,y), angular velocity, length, width, color and elasticity value

block1 = Softbody(np.array([200.0, 200.0]), np.array([30.0, 0.0]), 0.2, 10, 10, (250,220,1),90)
block2 = Softbody(np.array([550.0, 200.0]), np.array([0.0, 0.0]), 0, 10, 10, (20,220,1),90)

#To rotate softbody (radians) and to give softbody initial angular velocity, best used before runtime loop
block1.rotate(0.785)

pygame.display.set_caption('Soft Body Simulation')

while True:

       
    for event in pygame.event.get():
            if event.type == pygame.QUIT : 
                pygame.quit() 
                quit()
    win.fill(black)

    #Action is used in the runtime loop, for the softbody to exist
    block1.action()
    block2.action()

    #Applies collision between softbodies
    collide(block1,block2,0.75)
    pygame.display.update()


pygame.quit()
    

    
    
    




















    
