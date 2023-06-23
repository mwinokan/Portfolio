
"""

python3 -m http.server
localhost:8000/test.html

Make benzene a rigid body


Place them on a grid

Do gridded collisions?

NaOH
PO4
Sulfur?
Calcium hydroxide

"""

# debugging
debug = False
run = 2

# sim parameters
substeps = 8
dt = 0.0005
scale = 20
max_steps = 10000

# drawing
draw_indices = True
draw_electrostatics = False
draw_acceleration = True
draw_velocity = True

# forces
electrostatics = False
electrostatic_cutoff = 4
intramolecular_electrostatics = False
electrostatic_strength = 50000

# gravity
center_gravity = False
center_gravity_strength = 200

# linear gravity
linear_gravity = False
linear_gravity_strength = 10000

initial_velocity = 0.05

# # random
# random_acceleration = 1
# random_acceleration_strength = 100000

t = 0.0

import time
START = time.perf_counter()

import random

from js import document, Math, setInterval, clearInterval, setTimeout
import pyodide

import numpy as np

class Molecule:

	num_molecules = 0

	def __init__(self,name):
		self.name = name
		self.num_atoms = 0
		self.index = self.num_molecules
		self.num_molecules += 1
		self.atoms = []
		self.angles = []
		self.bonds = []

	def add_atom(self,atom):
		atom.index = self.num_atoms
		atom.mol_index = self.index
		self.num_atoms += 1
		self.atoms.append(atom)

	def add_angle(self,angle):
		self.angles.append(angle)

	def add_bond(self,i,j,order=1,length=None):
		if debug:
			print(f'Molecule.add_bond({i},{j},{order})')
		atom1 = self.atoms[i]
		atom2 = self.atoms[j]
		bond = Bond(atom1,atom2,order=order,length=length)
		self.bonds.append(bond)

	def print(self):
		print(f"Molecule {self.name}:")
		print("\nAtoms:")
		for atom in self.atoms:
			print(atom.index,atom.symbol)
		print("\nBonds:")
		for bond in self.bonds:
			print(bond.atom1.index,bond.atom2.index,bond.atom1.symbol,bond.atom2.symbol,bond.length)
		print("\nAngles:")
		for angle in self.angles:
			print(angle.atom1.index,angle.atom2.index,angle.atom3.index,angle.atom1.symbol,angle.atom2.symbol,angle.atom3.symbol,angle.angle)

class Bond:
	def __init__(self,atom1,atom2,order=1,length=None): 
		self.atom1 = atom1
		self.atom2 = atom2
		self.length = length or COVALENT_RADII[atom1.symbol] + COVALENT_RADII[atom2.symbol]
		self.atom1.bonds.append(self)
		self.atom2.bonds.append(self)

	def apply_constraint(self):

		v = self.atom1.position - self.atom2.position
		d = np.linalg.norm(v)
		n = v/d

		delta = self.length - d

		self.atom1.position += n * delta * 0.5
		self.atom2.position -= n * delta * 0.5

class Angle:
	def __init__(self,atom1,atom2,atom3,angle,strength=100,hard=True):
		self.atom1 = atom1
		self.atom2 = atom2
		self.atom3 = atom3
		self.angle = angle
		self.strength = strength
		if hard:
			self.apply_constraint = self.hard_constraint
		else:
			self.apply_constraint = self.soft_constraint

	def soft_constraint(self):

		v1 = self.atom1.position - self.atom2.position
		v2 = self.atom3.position - self.atom2.position

		d1 = np.linalg.norm(v1)
		d2 = np.linalg.norm(v2)

		angle = 180*np.arccos(np.dot(v1,v2)/(d1*d2))/np.pi

		delta = angle - self.angle

		if angle < self.angle:
			n1 = v1/d1
			n2 = v2/d2
			self.atom1.acceleration += self.strength*np.array([-n1[1],n1[0]])*delta
			self.atom3.acceleration += self.strength*np.array([n2[1],-n2[0]])*delta
		elif angle > self.angle:
			n1 = v1/d1
			n2 = v2/d2
			self.atom1.acceleration += self.strength*np.array([n1[1],-n1[0]])*delta
			self.atom3.acceleration += self.strength*np.array([-n2[1],n2[0]])*delta

	def hard_constraint(self,_print=False):

		v1 = self.atom1.position - self.atom2.position
		v2 = self.atom3.position - self.atom2.position

		d1 = np.linalg.norm(v1)
		d2 = np.linalg.norm(v2)

		v1_angle = vec_angle(v1)
		v2_angle = vec_angle(v2)

		bisector = v1/d1 + v2/d2
		bisector_d = np.linalg.norm(bisector)
		if bisector_d == 0.0:
			bisector_angle = Math.PI*(v1_angle+90)/180
			bisector = np.array([np.cos(bisector_angle),np.sin(bisector_angle)])
		else:
			bisector = bisector / bisector_d

		if _print:
			print('\n')

		if _print:
			print(self.atom1.index,self.atom2.index,self.atom3.index,self.atom1.symbol,self.atom2.symbol,self.atom3.symbol)

		if _print:
			draw_circle(ctx,self.atom1.position[0]*scale,self.atom1.position[1]*scale,self.atom1.r*scale,'green')
			draw_circle(ctx,self.atom3.position[0]*scale,self.atom3.position[1]*scale,self.atom3.r*scale,'green')

		# check bisector direction
		bisector_x_v1 = np.cross(bisector,v1)
		# bisector_x_v2 = np.cross(bisector,v2)
		if _print:
			print(f'{bisector_x_v1=}')
			# print(f'{bisector_x_v2=}')
		if bisector_x_v1 > 0:
			bisector = -bisector

		bisector_angle = vec_angle(bisector)

		if _print:
			print(f'{bisector_angle=}')
		if _print:
			print(f'{self.angle=}')

		if _print:
			print(f'{bisector=}')
			print(f'{v1_angle=}')
			print(f'{v2_angle=}')

		if bisector_x_v1 > 0:
			new_angle1 = (bisector_angle-(180-self.angle/2))%360.0
			new_angle2 = (bisector_angle+(180-self.angle/2))%360.0
		else:
			new_angle1 = (bisector_angle-self.angle/2)%360.0
			new_angle2 = (bisector_angle+self.angle/2)%360.0

		# if v1_angle < bisector_angle:
		# 	new_angle1 = (bisector_angle-self.angle/2)%360.0
		# 	new_angle2 = (bisector_angle+self.angle/2)%360.0
		# else:
		# 	new_angle2 = (bisector_angle-self.angle/2)%360.0
		# 	new_angle1 = (bisector_angle+self.angle/2)%360.0

		# new_angle1 = (bisector_angle-self.angle/2)%360.0
		# new_angle2 = (bisector_angle+self.angle/2)%360.0

		# if v1_angle > v2_angle and new_angle1 < new_angle2:
		# 	new_angle1,new_angle2 = new_angle2,new_angle1
		
		# atom1
		if _print:
			print(f'{new_angle1=}')
			print(f'{new_angle2=}')

		new_angle1 = Math.PI*new_angle1/180
		dx =  d1 * np.cos(new_angle1)
		dy =  d1 * np.sin(new_angle1)
		self.atom1.position = np.array([self.atom2.position[0]+dx,self.atom2.position[1]+dy])

		if _print:
			print(f'{dx=}')
		if _print:
			print(f'{dy=}')

		# atom2
		new_angle2 = Math.PI*new_angle2/180
		dx =  d2 * np.cos(new_angle2)
		dy =  d2 * np.sin(new_angle2)
		self.atom3.position = np.array([self.atom2.position[0]+dx,self.atom2.position[1]+dy])

		if _print:
			print(f'{dx=}')
		if _print:
			print(f'{dy=}')

		if _print:
			draw_line(ctx,(self.atom2.position[0])*scale,(self.atom2.position[1])*scale,(self.atom2.position[0]+bisector[0])*scale,(self.atom2.position[1]+bisector[1])*scale,c='green')

class Atom:
	# def __init__(self,symbol,x=0,y=0,q=0.0):
	def __init__(self,symbol,x=0,y=0,q=0.0,v=[0,0]):
		self.symbol = symbol
		self.bonds = []
		self.colour = COLOURS[symbol]
		self.m = MASSES[symbol]
		self.r = COVALENT_RADII[symbol]
		self.q = q
		if debug:
			print('Atom.__init__',symbol)
		self.position = np.array([x,y],dtype=float)
		# self.old_pos = np.array([x-v[0],y-v[1]],dtype=float)
		self.old_pos = np.array([x,y],dtype=float)
		self.acceleration = np.array([0.0,0.0],dtype=float)
		# self.acceleration = np.array([v[0],v[1]],dtype=float)
	@property
	def x(self):
		return self.position[0]
	@property
	def y(self):
		return self.position[1]
	def update(self,dt,v=None):
		# print(len(self.position),len(self.old_pos))
		# print(v is not None)
		if v is None:
			v = self.position - self.old_pos
		# print(np.linalg.norm(v))
		self.old_pos = self.position.copy()
		self.position = self.position + v + self.acceleration*dt
		# print(self.acceleration[0])
		if draw_acceleration:
			draw_line(ctx,self.x*scale,self.y*scale,self.x*scale+self.acceleration[0]*scale*dt,self.y*scale+self.acceleration[1]*scale*dt,w=10,c='yellow')
		if draw_velocity:
			draw_line(ctx,self.x*scale,self.y*scale,self.x*scale+v[0]*scale,self.y*scale+v[1]*scale,w=10,c='aqua')
		self.acceleration = np.array([0.0,0.0])

class Simulation:

	def __init__(self,canvas,ctx,scale=40):
		if debug:
			print('Setting up simulation...')
		self.scale = scale
		self.canvas = canvas
		self.ctx = ctx
		self.molecules = []

	def add_molecule(self,molecule):
		self.molecules.append(molecule)

	def draw(self):
		for molecule in self.molecules:
			for atom in molecule.atoms:
				draw_circle(self.ctx,atom.x*self.scale,atom.y*self.scale,atom.r*self.scale,atom.colour)
				if draw_indices:
					draw_text(f'{atom.index}',atom.x*self.scale,atom.y*self.scale)

			# for bond in molecule.bonds:
			# 	draw_line(self.ctx,bond.atom1.x*self.scale,bond.atom1.y*self.scale,bond.atom2.x*self.scale,bond.atom2.y*self.scale,w=2,c='black')

	def solve(self,dt,substeps=4,random_v=None):
		dt = dt/substeps
		for i in range(substeps):
			if linear_gravity: 
				self.apply_gravity([0,linear_gravity_strength])
			# if random_acc:
			# 	sim.random_acceleration(random_acceleration_strength)
			if center_gravity:
				self.apply_center_gravity(center_gravity_strength)
			if electrostatics:
				self.apply_electrostatics(electrostatic_strength)
			self.apply_bond_constraints()
			self.check_collisions(1.0)
			self.apply_angle_constraints()
			self.clip()
			self.update_positions(dt,random_v=random_v)
			random_v = None
			pass


	def update_positions(self,dt,random_v=None):
		for molecule in self.molecules:
			if random_v is not None:
				random_v = random_velocity(random_v)
				if debug:
					print(random_v)
			for atom in molecule.atoms:
				atom.update(dt,v=random_v)

	def apply_gravity(self,g):
		g = np.array(g)
		for molecule in self.molecules:
			for atom in molecule.atoms:
				atom.acceleration += g

	def apply_center_gravity(self,g):
		for molecule in self.molecules:
			for atom in molecule.atoms:
				atom.acceleration += -g*atom.position

	def apply_electrostatics(self,strength):

		if intramolecular_electrostatics:
			# apply intra-molecularly
			for mol1 in self.molecules:
				for i,atom1 in enumerate(mol1.atoms):
					for atom2 in mol1.atoms[i+1:]:
						v = atom1.position-atom2.position
						d = np.linalg.norm(v)
						if d > 0:
							F = -strength*atom1.q*atom2.q / d**2
							if F > 0:
								continue

							atom1.acceleration
							n = F * v / d

							atom1.acceleration -= n / atom1.m
							atom2.acceleration += n / atom2.m

							# if mol1.index == 0:
							# 	print(atom1.symbol,atom2.symbol,atom1.q,atom2.q,F)

		# apply inter-molecularly
		for m,mol1 in enumerate(self.molecules):
			for mol2 in self.molecules[m+1:]:
				for atom1 in mol1.atoms:
					for atom2 in mol2.atoms:
						v = atom1.position-atom2.position
						d = np.linalg.norm(v)
						if d > 0 and d < electrostatic_cutoff:
							F = -strength*atom1.q*atom2.q / d**2
							atom1.acceleration
							n = F * v / d

							atom1.acceleration -= n / atom1.m
							atom2.acceleration += n / atom2.m

							# print(f'q_prod={atom1.q*atom2.q},{d=},{F=}')

							if draw_electrostatics:
								if F > 0:
									# repulsive
									draw_line(ctx,atom1.x*scale,atom1.y*scale,atom2.x*scale,atom2.y*scale,w=F/10,c='red')
								else:
									draw_line(ctx,atom1.x*scale,atom1.y*scale,atom2.x*scale,atom2.y*scale,w=F/10,c='green')

							# print(atom1.symbol,atom2.symbol,atom1.q,atom2.q,F)

	def check_collisions(self,c_r=0.75):
		# apply intra-molecularly
		for mol1 in self.molecules:
			for i,atom1 in enumerate(mol1.atoms):
				for atom2 in mol1.atoms[i+1:]:

					v = atom1.position - atom2.position
					d = np.linalg.norm(v)

					r_sum = atom1.r + atom2.r

					if d == 0:
						atom1.position += np.array([0.1,0])
						atom2.position -= np.array([0.1,0])
						v = np.array([0.2,0])
						d = 0.2

					if d < r_sum:

						n = v / d

						m_sum = atom1.m + atom2.m

						m_ratio_1 = atom1.m / m_sum
						m_ratio_2 = atom2.m / m_sum

						delta = r_sum - d

						atom1.position += n * m_ratio_2 * delta * 0.5 * c_r
						atom2.position -= n * m_ratio_1 * delta * 0.5 * c_r

		# only apply inter-molecularly
		for m,mol1 in enumerate(self.molecules):
			for mol2 in self.molecules[m+1:]:
				for atom1 in mol1.atoms:
					for atom2 in mol2.atoms:

						v = atom1.position - atom2.position
						d = np.linalg.norm(v)

						r_sum = atom1.r + atom2.r

						if d == 0:
							atom1.position += np.array([0.1,0])
							atom2.position -= np.array([0.1,0])
							v = np.array([0.2,0])
							d = 0.2

						if d < r_sum:

							n = v / d

							m_sum = atom1.m + atom2.m

							m_ratio_1 = atom1.m / m_sum
							m_ratio_2 = atom2.m / m_sum

							delta = r_sum - d

							atom1.position += n * m_ratio_2 * delta * 0.5 * c_r
							atom2.position -= n * m_ratio_1 * delta * 0.5 * c_r

	def apply_bond_constraints(self):
		for molecule in self.molecules:
			for bond in molecule.bonds:
				bond.apply_constraint()
				pass

	def apply_angle_constraints(self):
		for molecule in self.molecules:
			for angle in molecule.angles:
				angle.apply_constraint()
				pass

	def clip(self):
		for molecule in self.molecules:
			for atom in molecule.atoms:
				# print(atom.x,atom.y,w,h)
				w = self.width/self.scale/2
				h = self.height/self.scale/2
				if atom.x + atom.r > w:
					atom.position = np.array([w-atom.r,atom.y])
				elif atom.x - atom.r < -w:
					atom.position = np.array([-w+atom.r,atom.y])

				if atom.y + atom.r > h:
					atom.position = np.array([atom.x,h-atom.r])
				elif atom.y - atom.r < -h:
					atom.position = np.array([atom.x,-h+atom.r])

def vec_angle(v):
	v = v / np.linalg.norm(v)

	if v[0] == 0:
		if v[1] > 0:
			v_angle = 90
		else:
			v_angle = -90
	else:
		v_angle = 180*np.arctan(v[1]/v[0])/Math.PI
		v_angle = 180*np.arccos(v[0])/Math.PI
		if v[1] < 0:
			v_angle = 360 - v_angle

	# if v_angle > 359.5:
	# 	v_angle += 0.5

	return v_angle%360.0

def set_running():
	document.getElementById("py-status").innerHTML = ''

def draw_circle(ctx,x,y,r,c):
	ctx.beginPath()
	ctx.arc(x,y,r,0,2*Math.PI)
	ctx.fillStyle = c
	ctx.fill()

def draw_line(ctx,x1,y1,x2,y2,w=2,c='white'):
	ctx.beginPath()
	ctx.moveTo(x1,y1)
	ctx.lineTo(x2,y2)
	ctx.lineWidth = w
	ctx.strokeStyle = c
	ctx.stroke()

def draw_text(t,x,y,s=10,c='black'):
	ctx.fillStyle = c
	ctx.font = f"{s}px Arial"
	ctx.fillText(t, x, y)

def draw_axes(ctx,canvas):
	ctx.lineWidth = 2
	draw_line(ctx,0,0,0,-canvas.height/2)
	draw_line(ctx,0,0,0,canvas.height/2)
	draw_line(ctx,0,0,canvas.width/2,0)
	draw_line(ctx,0,0,-canvas.width/2,0)

def main():

	global interval_id, t, sim, ctx, canvas, dt, substeps
	set_running()
	
	canvas = document.getElementById("headerCanvas")
	ctx = canvas.getContext("2d")

	if debug:
		print(f'{canvas.width=}')
		print(f'{canvas.height=}')

	ctx.translate(canvas.width/2,canvas.height/2)

	sim = Simulation(canvas,ctx,scale=scale)

	sim.width = canvas.width
	sim.height = canvas.height

	# make_water_grid(sim,7,5,4)
	# make_water_grid(sim,7,5,4)

	# sim.add_molecule(make_H2O(x=-3,y=3))
	# sim.add_molecule(make_H2O(x=3,y=-3))
	sim.add_molecule(make_H2O(x=-3,y=-3))
	# # sim.add_molecule(make_H2O(x=3,y=3))
	
	sim.add_molecule(make_NH3(x=3,y=3))
	
	sim.add_molecule(make_CO2(x=5,y=5))
	
	sim.add_molecule(make_O3(x=-5,y=5))
	
	sim.add_molecule(make_NH2(x=5,y=-5))
	
	sim.add_molecule(make_CH4(x=8,y=-5))
	
	# # wobbly
	# sim.add_molecule(make_C2H4(x=-5,y=-5))

	# sim.molecules[0].angles[0].apply_constraint(_print=True)
	
	# sim.add_molecule(make_benzene(x=0,y=0))
	
	# # benzene angle debug:
	# j = 12
	# for i in range(j):
	# 	sim.molecules[0].angles[i].apply_constraint(_print=i==j-1)
	
	sim.draw()
	
	if run:
		clear_screen()
		# sim.solve(dt=dt,substeps=substeps,random_acc=random_acceleration>0)
		sim.solve(dt=dt,substeps=substeps,random_v=initial_velocity)
		sim.draw()


	draw_loop_proxy = pyodide.ffi.create_proxy(draw_loop)

	if run > 1:
		interval_id = setInterval(draw_loop_proxy,50)

	if debug:
		print('finished.')
	
counter = 0

def make_water_grid(sim,c,r,d=4):
	for i in range(c*r):
		x = d*(i%c) - d*c//2
		y = d*(i//c) - d*r//2
		# print(i,x,y)

		sim.add_molecule(make_H2O(x,y))

def make_H2O(x,y):
	if debug:
		print('making water')
	h2o = Molecule('water')
	h2o.add_atom(Atom('O',x,y,q=-0.834))
	h2o.add_atom(Atom('H',x+1,y,q=0.417))
	h2o.add_atom(Atom('H',x-1,y,q=0.417))
	h2o.add_bond(0,1)
	h2o.add_bond(0,2)
	h2o.add_angle(Angle(h2o.atoms[1],h2o.atoms[0],h2o.atoms[2],BOND_ANGLES['H O H']))
	return h2o
	# def make_H2O(x,y):

def make_CO2(x,y):
	if debug:
		print('making carbon dioxide')
	co2 = Molecule('carbon dioxide')
	co2.add_atom(Atom('C',x,y,q=0.6))
	co2.add_atom(Atom('O',x+1,y,q=-0.3))
	co2.add_atom(Atom('O',x-1,y,q=-0.3))
	co2.add_bond(0,1)
	co2.add_bond(0,2)
	co2.add_angle(Angle(co2.atoms[1],co2.atoms[0],co2.atoms[2],BOND_ANGLES['O C O']))
	return co2

def make_CO(x,y):
	if debug:
		print('making carbon monoxide')
	co = Molecule('carbon monoxide')
	co.add_atom(Atom('C',x,y,q=0.6))
	co.add_atom(Atom('O',x+1,y,q=-0.3))
	co.add_bond(0,1)
	return co

# no charges!
def make_NH2(x,y):
	if debug:
		print('making amidogen')
	if electrostatics:
		print('amidogen: inaccurate charges!!')
	nh2 = Molecule('amidogen')
	nh2.add_atom(Atom('N',x,y,q=0.6))
	nh2.add_atom(Atom('H',x+1,y,q=-0.3))
	nh2.add_atom(Atom('H',x-1,y,q=-0.3))
	nh2.add_bond(0,1)
	nh2.add_bond(0,2)
	nh2.add_angle(Angle(nh2.atoms[1],nh2.atoms[0],nh2.atoms[2],BOND_ANGLES['H N H']))
	return nh2

# no charges!
def make_NH3(x,y):
	if debug:
		print('making ammonia')
	if electrostatics:
		print('ammonia: inaccurate charges!!')
	nh3 = Molecule('ammonia')
	nh3.add_atom(Atom('N',x,y,q=0.6))
	nh3.add_atom(Atom('H',x+1,y-0.5,q=-0.3))
	nh3.add_atom(Atom('H',x-1,y,q=-0.3))
	nh3.add_atom(Atom('H',x,y+1,q=-0.3))
	nh3.add_bond(0,1)
	nh3.add_bond(0,2)
	nh3.add_bond(0,3)
	nh3.add_angle(Angle(nh3.atoms[1],nh3.atoms[0],nh3.atoms[2],120))
	nh3.add_angle(Angle(nh3.atoms[3],nh3.atoms[0],nh3.atoms[1],120))
	return nh3

# no charges!
def make_CH4(x,y):
	if debug:
		print('making methane')
	if electrostatics:
		print('methane: inaccurate charges!!')
	ch4 = Molecule('methane')
	ch4.add_atom(Atom('C',x,y,q=0.6))
	ch4.add_atom(Atom('H',x+1,y-1,q=-0.3))
	ch4.add_atom(Atom('H',x-1,y-1,q=-0.3))
	ch4.add_atom(Atom('H',x-1,y+1,q=-0.3))
	ch4.add_atom(Atom('H',x+1,y+1,q=-0.3))
	ch4.add_bond(0,1)
	ch4.add_bond(0,2)
	ch4.add_bond(0,3)
	ch4.add_bond(0,4)
	ch4.add_angle(Angle(ch4.atoms[1],ch4.atoms[0],ch4.atoms[2],90))
	ch4.add_angle(Angle(ch4.atoms[2],ch4.atoms[0],ch4.atoms[3],90))
	ch4.add_angle(Angle(ch4.atoms[3],ch4.atoms[0],ch4.atoms[4],90))
	# ch4.add_angle(Angle(ch4.atoms[3],ch4.atoms[0],ch4.atoms[1],120))
	return ch4

# no charges!
def make_O3(x,y):
	if debug:
		print('making ozone')
	if electrostatics:
		print('ozone: inaccurate charges!!')
	ozone = Molecule('ozone')
	ozone.add_atom(Atom('O',x,y,q=0.6))
	ozone.add_atom(Atom('O',x+1,y,q=-0.3))
	ozone.add_atom(Atom('O',x-1,y,q=-0.3))
	ozone.add_bond(0,1)
	ozone.add_bond(0,2)
	ozone.add_angle(Angle(ozone.atoms[1],ozone.atoms[0],ozone.atoms[2],120))
	return ozone

# wobbly
def make_C2H4(x,y):
	if debug:
		print('making ethylene')
	ethe = Molecule('ethylene')
	ethe.add_atom(Atom('C',-0.5+x,y,q=-0.42))
	ethe.add_atom(Atom('C',-0.5+x+1,y,q=-0.42))
	ethe.add_atom(Atom('H',-0.5+x-0.5,y+0.5,q=0.21))
	ethe.add_atom(Atom('H',-0.5+x-0.5,y-0.5,q=0.21))
	ethe.add_atom(Atom('H',-0.5+x+1.5,y+0.5,q=0.21))
	ethe.add_atom(Atom('H',-0.5+x+1.5,y-0.5,q=0.21))
	ethe.add_bond(0,1)
	ethe.add_bond(0,2)
	ethe.add_bond(0,3)
	ethe.add_bond(1,4)
	ethe.add_bond(1,5)

	ethe.add_bond(2,4,length=3)

	# ethe.print()
	
	# all of these make it unstable!
	# ethe.add_angle(Angle(ethe.atoms[0],ethe.atoms[1],ethe.atoms[4],BOND_ANGLES['C C H']))
	# ethe.add_angle(Angle(ethe.atoms[5],ethe.atoms[1],ethe.atoms[4],BOND_ANGLES['C C H']))
	# ethe.add_angle(Angle(ethe.atoms[1],ethe.atoms[0],ethe.atoms[2],BOND_ANGLES['C C H']))
	# ethe.add_angle(Angle(ethe.atoms[3],ethe.atoms[0],ethe.atoms[1],BOND_ANGLES['C C H']))
	
	ethe.add_angle(Angle(ethe.atoms[2],ethe.atoms[0],ethe.atoms[3],BOND_ANGLES['H C H']))
	ethe.add_angle(Angle(ethe.atoms[4],ethe.atoms[1],ethe.atoms[5],BOND_ANGLES['H C H']))
	return ethe

# very unstable
def make_benzene(x,y):
	if debug:
		print('making benzene')

	r = 1.5
	
	benz = Molecule('benzene')
	for i in range(6):

		x = r*np.cos(2*Math.PI*i/6)
		y = r*np.sin(2*Math.PI*i/6)
		benz.add_atom(Atom('C',x,y,q=-0.115))

		x = (r+1)*np.cos(2*Math.PI*i/6)
		y = (r+1)*np.sin(2*Math.PI*i/6)
		benz.add_atom(Atom('H',x,y,q=0.115))

		benz.add_bond(2*i,2*i+1)

	for i in range(6):
		benz.add_bond(2*i-2,2*i)

	for i in range(6):
		# print(2*i-2,2*i,(2*i+2)%12)

		# if i%2 == 0:
		# C C C
		benz.add_angle(Angle(benz.atoms[2*i-2],benz.atoms[2*i],benz.atoms[(2*i+2)%12],120,500))

		# break
		
	# for i in range(5):
	# 	# H C C
	# 	benz.add_angle(Angle(benz.atoms[2*i+1],benz.atoms[2*i],benz.atoms[(2*i+2)%12],120,100))
		
		# # C C H
		# benz.add_angle(Angle(benz.atoms[2*i-2],benz.atoms[2*i],benz.atoms[(2*i+1)%12],103.3,100))
	
	# benz.print()

	return benz

def clear_screen():
	ctx.clearRect(-canvas.width/2, -canvas.height/2, canvas.width, canvas.height)

def draw_loop():
	global ctx, counter, t

	start = time.perf_counter

	clear_screen()
	ctx.fillStyle = 'white'
	ctx.font = "30px Arial"
	# ctx.fillText(counter, -canvas.width/2+10, 50)
	ctx.fillText(f'{t=:.3f}', -canvas.width/2+10, -canvas.height/2+50)

	draw_axes(ctx,canvas)
	sim.solve(dt=dt,substeps=substeps)
	sim.draw()

	t += dt
	counter += 1
	
	if counter == max_steps:
		clearInterval(interval_id)

def random_velocity(strength):
	theta = random.uniform(0,2*Math.PI)
	return strength*np.array([np.cos(theta),np.sin(theta)])

COVALENT_RADII = {
'X': 0.2, 'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76, 'Sc': 1.7, 'Ti': 1.6, 'V': 1.53, 'Cr': 1.39, 'Mn': 1.39, 'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 'Ga': 1.22, 'Ge': 1.2, 'As': 1.19, 'Se': 1.2, 'Br': 1.2, 'Kr': 1.16, 'Rb': 2.2, 'Sr': 1.95, 'Y': 1.9, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.4, 'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89, 'Tm': 1.9, 'Yb': 1.87, 'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.7, 'W': 1.62, 'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.4, 'At': 1.5, 'Rn': 1.5, 'Fr': 2.6, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06, 'Pa': 2.0, 'U': 1.96, 'Np': 1.9, 'Pu': 1.87, 'Am': 1.8, 'Cm': 1.69, 'Bk': 0.2, 'Cf': 0.2, 'Es': 0.2, 'Fm': 0.2, 'Md': 0.2, 'No': 0.2, 'Lr': 0.2, 'Rf': 0.2, 'Db': 0.2, 'Sg': 0.2, 'Bh': 0.2, 'Hs': 0.2, 'Mt': 0.2, 'Ds': 0.2, 'Rg': 0.2, 'Cn': 0.2, 'Nh': 0.2, 'Fl': 0.2, 'Mc': 0.2, 'Lv': 0.2, 'Ts': 0.2, 'Og': 0.2
}

MASSES = {
'X': 1.0, 'H': 1.008, 'He': 4.002602, 'Li': 6.94, 'Be': 9.0121831, 'B': 10.81, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998403163, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815385, 'Si': 28.085, 'P': 30.973761998, 'S': 32.06, 'Cl': 35.45, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078, 'Sc': 44.955908, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938044, 'Fe': 55.845, 'Co': 58.933194, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.723, 'Ge': 72.63, 'As': 74.921595, 'Se': 78.971, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90584, 'Zr': 91.224, 'Nb': 92.90637, 'Mo': 95.95, 'Tc': 97.90721, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42, 'Ag': 107.8682, 'Cd': 112.414, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.6, 'I': 126.90447, 'Xe': 131.293, 'Cs': 132.90545196, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116, 'Pr': 140.90766, 'Nd': 144.242, 'Pm': 144.91276, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93033, 'Er': 167.259, 'Tm': 168.93422, 'Yb': 173.054, 'Lu': 174.9668, 'Hf': 178.49, 'Ta': 180.94788, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217, 'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.592, 'Tl': 204.38, 'Pb': 207.2, 'Bi': 208.9804, 'Po': 208.98243, 'At': 209.98715, 'Rn': 222.01758, 'Fr': 223.01974, 'Ra': 226.02541, 'Ac': 227.02775, 'Th': 232.0377, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.04817, 'Pu': 244.06421, 'Am': 243.06138, 'Cm': 247.07035, 'Bk': 247.07031, 'Cf': 251.07959, 'Es': 252.083, 'Fm': 257.09511, 'Md': 258.09843, 'No': 259.101, 'Lr': 262.11, 'Rf': 267.122, 'Db': 268.126, 'Sg': 271.134, 'Bh': 270.133, 'Hs': 269.1338, 'Mt': 278.156, 'Ds': 281.165, 'Rg': 281.166, 'Cn': 285.177, 'Nh': 286.182, 'Fl': 289.19, 'Mc': 289.194, 'Lv': 293.204, 'Ts': 293.208, 'Og': 294.214
}

BOND_ANGLES = {
'H N C': 111.0, 'N C C': 110.0, 'C C O': 110.5, 'C C C': 108.0, 'N C H': 109.5, 'C N C': 112.0, 'H C C': 109.5, 'C O C': 109.6, 'C S C': 95.0, 'H N H': 120.0, 'H O C': 108.0, 'H C H': 109.0, 'H S C': 95.0, 'N C N': 120.0, 'O C C': 118.0, 'O C H': 121.7, 'O C N': 122.5, 'O C O': 124.0, 'S C C': 112.5, 'S C H': 111.3, 'S S C': 103.3, 'C C H': 110.1, 'C C N': 122.0, 'H O H': 104.5, 'C C S': 124.0, 'N C S': 116.4, 'N C O': 124.0, 'S C S': 124.0, 'O C S': 125.0, 'C C Cl': 120.0, 'C C Br': 120.0, 'C C I': 120.0, 'N C Br': 120.0, 'C C F': 118.8, 'Cl C Cl': 109.0, 'Br C Br': 110.5, 'F C F': 107.0, 'Cl C H': 108.5, 'Br C H': 107.0, 'C C P': 117.0, 'P C F': 122.0, 'F C H': 108.9, 'P C H': 110.0, 'C N N': 115.0, 'C N H': 113.0, 'C N O': 116.0, 'O N O': 128.0, 'C N S': 111.0, 'N N N': 102.2, 'N N O': 103.0, 'N N H': 119.5, 'C N P': 118.3, 'P N H': 123.6, 'S N H': 113.1, 'C O N': 108.5, 'C O P': 120.0, 'C O S': 108.0, 'P O P': 143.0, 'C O H': 115.0, 'N O H': 101.5, 'P O H': 115.0, 'O P O': 111.6, 'O P S': 131.0, 'S P S': 128.5, 'C P O': 94.0, 'N P O': 110.6, 'C S N': 103.0, 'C S S': 103.3, 'C S H': 95.0, 'C S O': 98.0, 'N S O': 94.2, 'O S O': 109.5, 'N S N': 102.3, 'F Al F': 109.5, 'H C O': 108.9, 'H O P': 115.0, 'O S N': 103.0, 'S N C': 113.0, 'H C N': 113.5, 'C S Fe': 100.6, 'S Fe N': 90.0, 'Fe N C': 128.1, 'H C Fe': 180.0, 'N Fe C': 90.0, 'N Fe N': 90.0, 'O C Fe': 180.0, 'O Fe N': 90.0, 'O O Fe': 180.0, 'H C P': 110.0, 'S Zn S': 111.8, 'C S Zn': 95.0, 'C N Zn': 120.7, 'S Zn N': 108.1, 'N Zn N': 107.8, 'F C C': 118.8, 'H N P': 123.6, 'N C Se': 122.5, 'C O O': 104.0, 'O O H': 98.3, 'O Si O': 117.0, 'Si O Si': 150.5, 'H Si H': 119.0, 'H Si O': 118.0, 'Si O H': 122.5, 'Al O Al': 98.0, 'Al O Si': 117.0, 'H O Al': 93.4, 'O Al H': 93.4, 'O Al O': 90.0, 'O Si H': 118.0, 'O S C': 99.0
}

COLOURS = {'H':'Azure','O':'Crimson','C':'LightSlateGrey','N':'RoyalBlue'}

main()
