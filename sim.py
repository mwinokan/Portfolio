#!/usr/bin/env python3

"""

* Stop collisions adding so much energy

"""

# debugging
debug = False
run = 2

# sim parameters
substeps = 8
# dt = 0.0005
dt = 0.0001
scale = 20
max_steps = 10000
initial_velocity = 0.02
avg_velocity = 0.02
# avg_velocity = None
v_adjust = 0.01
ang_v_decay = 0.995

# drawing
draw_orientations = False
draw_hitcircles = False

draw_indices = False
draw_acceleration = False
draw_velocity = False
draw_axes_lines = False
draw_time = False

t = 0.0

import math
import time
START = time.perf_counter()

import random

try:
	from js import document, Math, setInterval, clearInterval, setTimeout
	import pyodide
except ModuleNotFoundError:
	print('COULD NOT LOAD WEB MODULES')
	exit()

import numpy as np

class Molecule:

	num_molecules = 0

	def __init__(self,name,x,y,theta):

		self.index = self.num_molecules
		self.num_molecules += 1
		
		self.num_atoms = 0
		self.atoms = []

		self.name = name

		self.p = np.array([x,y],dtype=float)
		self.p_old = np.array([x,y],dtype=float)
		self.t = theta
		self.t_old = theta

		self._m = None
		self._m_r = None
		self._r = None
		self._v = None
		
		self.a = np.array([0.0,0.0])
		self.ang_v = None
		self.ang_a = 0.0

	def center_atoms(self):

		# calculate the CoM
		x_center = sum([a.x*a.m for a in self.atoms])/self.m
		y_center = sum([a.y*a.m for a in self.atoms])/self.m
		
		# shift atoms by the CoM
		for a in self.atoms:
			a.x -= x_center
			a.y -= y_center

		# force a recalculation of the molecule radius
		self._r = None
		self._m_r = None

	@property
	def x(self):
		return self.p[0]
	
	@property
	def y(self):
		return self.p[1]

	@property
	def m(self):
		if not self._m:
			self._m = sum(a.m for a in self.atoms)
		return self._m

	@property
	def r(self):
		if not self._r:
			self._r = max(a.r+a.d for a in self.atoms)
		return self._r

	@property
	def m_r(self):
		if not self._m_r:
			self._m_r = sum(a.m*a.d**2 for a in self.atoms)
		return self._m_r

	def add_atom(self,atom):
		atom.index = self.num_atoms
		atom.mol_index = self.index
		self.num_atoms += 1
		self.atoms.append(atom)

	def draw(self):
		if draw_hitcircles:
			draw_circle_stroke(self.x,self.y,self.r,'aqua')
		for atom in self.atoms:
			x = self.x + atom.d * np.cos(atom.t+self.t)
			y = self.y + atom.d * np.sin(atom.t+self.t)
			draw_circle(x,y,atom.r,atom.c)
		if draw_orientations:
			dx = self.r * np.cos(self.t)
			dy = self.r * np.sin(self.t)
			draw_line(self.x,self.y,self.x+dx,self.y+dy,w=5,c='yellow')

	def print(self):
		print(f"Molecule {self.name}:")
		print("\nAtoms:")
		for atom in self.atoms:
			print(atom.index,atom.symbol,f'{atom.d=:.2f} {atom.t=:.2f} {atom.x=:.2f} {atom.y=:.2f}')

	def get_velocity(self,recalc=True):
		if self._v is None or recalc:
			self._v = self.p - self.p_old
		return self._v

	def update(self,dt,v=None,recalc=True,v_scale=1.0):

		# calculate linear velocity
		if v is None:
			v = self.get_velocity(recalc=recalc)
			if v_scale is not None:
				if v_scale > 1.0:
					v *= 1.0 + v_adjust
				elif v_scale < 1.0:
					v *= 1.0 - v_adjust

		# angular velocity
		ang_v = self.t - self.t_old
		ang_v *= ang_v_decay

		# save current position
		self.p_old = self.p.copy()
		self.t_old = float(self.t)

		# update position
		self.p = self.p + v + self.a*dt
		self.t = self.t + ang_v + self.ang_a*dt

		if debug:
			print(self.index,f'p: {self.p_old} --> {self.p}, t: {self.t_old:.1f} --> {self.t:.1f}')

		# zero accelerations
		self.a = np.array([0.0,0.0])
		self.ang_a = 0.0

	def check_clip(self,axis,boundary):

		# axis should be a 2-vector of the clipping direction

		# clipping with positive x (right)
		if axis[0] > 0:
			for atom in self.atoms:
				x = atom.r + self.x + atom.d * np.cos(atom.t+self.t)
				if x > boundary:
					delta = x - boundary
					self.p -= np.array([2.0*delta,0.0])

					# angular
					lever = atom.d * np.array([np.cos(atom.t+self.t),np.sin(atom.t+self.t)])
					delta = np.array([delta,0.0])
					torque = cw_normal(lever)
					self.t += np.dot(delta,torque) / np.linalg.norm(torque)
					break

		# clipping with negative x (left)
		elif axis[0] < 0:
			for atom in self.atoms:
				x = -atom.r + self.x + atom.d * np.cos(atom.t+self.t)
				if x < boundary:
					delta = x - boundary
					self.p -= np.array([2.0*delta,0.0])

					# angular
					lever = atom.d * np.array([np.cos(atom.t+self.t),np.sin(atom.t+self.t)])
					delta = np.array([delta,0.0])
					torque = cw_normal(lever)
					self.t += np.dot(delta,torque) / np.linalg.norm(torque)
					break

		# clipping with positive y (bottom)
		if axis[1] > 0:
			for atom in self.atoms:
				y = atom.r + self.y + atom.d * np.sin(atom.t+self.t)
				if y > boundary:
					delta = y - boundary
					self.p -= np.array([0.0,2.0*delta])
					
					# angular
					lever = atom.d * np.array([np.cos(atom.t+self.t),np.sin(atom.t+self.t)])
					delta = np.array([0.0,delta])
					torque = cw_normal(lever)
					self.t += np.dot(delta,torque) / np.linalg.norm(torque)
					break

		# clipping with negative y (top)
		elif axis[1] < 0:
			for atom in self.atoms:
				y = -atom.r + self.y + atom.d * np.sin(atom.t+self.t)
				if y < boundary:
					delta = y - boundary
					self.p -= np.array([0.0,2.0*delta])
					
					# angular
					lever = atom.d * np.array([np.cos(atom.t+self.t),np.sin(atom.t+self.t)])
					delta = np.array([0.0,delta])
					torque = cw_normal(lever)
					self.t += np.dot(delta,torque) / np.linalg.norm(torque)
					break

class Atom:
	def __init__(self,symbol,xy=[0.0,0.0],rt=None,q=0.0):

		if debug:
			print(f'Atom.__init__({xy=},{rt=})')

		self.symbol = symbol

		self.c = COLOURS[symbol]
		self.m = MASSES[symbol]
		self.r = COVALENT_RADII[symbol]
		self.q = q

		if rt is not None:
			self._d, self.t = rt
			self._x = self.d*np.cos(self.t)
			self._y = self.d*np.sin(self.t)
			self._p = np.array([self.x,self.y])
		else:
			self._x, self._y = xy
			self._p = np.array([self.x,self.y])
			self._d = np.linalg.norm(self.p)
			self._t = vec_angle(np.array(self.p))

	@property
	def x(self):
		return self._x

	@property
	def y(self):
		return self._y

	@property
	def d(self):
		return self._d

	@property
	def t(self):
		return self._t

	@property
	def p(self):
		return self._p

	@x.setter
	def x(self,a):
		self._x = a
		self._p = np.array([a,self.y])
		self._d = np.linalg.norm(self.p)
		self._t = vec_angle(np.array(self.p))

	@y.setter
	def y(self,a):
		self._y = a
		self._p = np.array([self.x,a])
		self._d = np.linalg.norm(self.p)
		self._t = vec_angle(np.array(self.p))

	@d.setter
	def d(self,a):
		self._d = a
		self._x = self.d*np.cos(self.t)
		self._y = self.d*np.sin(self.t)
		self._p = np.array([self.x,self.y])

	@t.setter
	def t(self,a):
		self._t = a
		self._x = self.d*np.cos(self.t)
		self._y = self.d*np.sin(self.t)
		self._p = np.array([self.x,self.y])
		
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
			molecule.draw()
			# for atom in molecule.atoms:
			# 	draw_circle(self.ctx,atom.x*self.scale,atom.y*self.scale,atom.r*self.scale,atom.colour)
			# 	if draw_indices:
			# 		draw_text(f'{atom.index}',atom.x*self.scale,atom.y*self.scale)

	def solve(self,dt,substeps=4,random_v=None):
		dt = dt/substeps
		for i in range(substeps):
			# self.apply_bond_constraints()
			self.check_collisions(0.5)
			# self.apply_angle_constraints()
			self.clip()
			current_avg_velocity = self.avg_velocity()
			self.update_positions(dt,random_v=random_v,recalc=False,current_avg_velocity=current_avg_velocity)
			random_v = None
			pass

	def avg_velocity(self):
		v_sum = 0.0
		count = 0
		for molecule in self.molecules:
			count += 1
			v = molecule.get_velocity(recalc=True)
			try:
				v_sum += np.linalg.norm(v)
			except:
				print(v,type(v))
				raise Exception
		return v_sum / count

	def update_positions(self,dt,random_v=None,recalc=True,current_avg_velocity=None):
		if random_v is not None:
			v_scale = None
		elif current_avg_velocity is not None and avg_velocity is not None:
			v_scale = avg_velocity / current_avg_velocity
		else:
			v_scale = None
		for molecule in self.molecules:
			if random_v is not None:
				v = random_velocity(random_v)
			else:
				v = None
			# 	if debug:
			# 		print(random_v)
			# for atom in molecule.atoms:
			# 	atom.update(dt,v=random_v,recalc=recalc,v_scale=v_scale)
			molecule.update(dt,v=v,recalc=recalc,v_scale=v_scale)
			pass

	def molecule_collision(self,mol1,mol2,c_r=0.75):
		for atom1 in mol1.atoms:
			x1 = mol1.x + atom1.d*np.cos(atom1.t + mol1.t)
			y1 = mol1.y + atom1.d*np.sin(atom1.t + mol1.t)
			
			for atom2 in mol2.atoms:
				x2 = mol2.x + atom2.d*np.cos(atom2.t + mol2.t)
				y2 = mol2.y + atom2.d*np.sin(atom2.t + mol2.t)

				v = np.array([x1-x2,y1-y2])
				d = np.linalg.norm(v)

				r_sum = atom1.r + atom2.r

				if d == 0:
					mol1.p += np.array([0.1,0])
					mol2.p -= np.array([0.1,0])
					v = np.array([0.2,0])
					d = 0.2

				if d < r_sum:

					# linear

					n = v / d

					m_sum = mol1.m + mol2.m

					m_ratio_1 = mol1.m / m_sum
					m_ratio_2 = mol2.m / m_sum

					delta = r_sum - d
					
					mol1.p += n * m_ratio_2 * delta * 0.5 * c_r
					mol2.p -= n * m_ratio_1 * delta * 0.5 * c_r

					# angular
					d_vec = 0.5 * delta * n

					if atom1.d > 0:
						lever_1 = atom1.d * np.array([np.cos(atom1.t + mol1.t),np.sin(atom1.t + mol1.t)])
						torque_1 = cw_normal(lever_1)
						mol1.t += np.dot(d_vec,torque_1) / (np.linalg.norm(torque_1)*mol1.m_r)

					if atom2.d > 0:
						lever_2 = atom2.d * np.array([np.cos(atom2.t + mol2.t),np.sin(atom2.t + mol2.t)])
						torque_2 = cw_normal(lever_2)
						mol2.t += np.dot(-d_vec,torque_2) / (np.linalg.norm(torque_2)**mol2.m_r)
					
					return

	def check_collisions(self,c_r=0.75):

		for m,mol1 in enumerate(self.molecules):
			for mol2 in self.molecules[m+1:]:

				v = mol1.p - mol2.p
				d = np.linalg.norm(v)

				r_sum = mol1.r + mol2.r

				if d == 0:
					mol1.p += np.array([0.1,0])
					mol2.p -= np.array([0.1,0])
					v = np.array([0.2,0])
					d = 0.2

				if d < r_sum:

					# # apply using the atom radii
					self.molecule_collision(mol1,mol2,c_r=c_r)

					# # apply using the molecular raii
					
					# n = v / d

					# m_sum = mol1.m + mol2.m

					# m_ratio_1 = mol1.m / m_sum
					# m_ratio_2 = mol2.m / m_sum

					# delta = r_sum - d

					# mol1.p += n * m_ratio_2 * delta * 0.5 * c_r
					# mol2.p -= n * m_ratio_1 * delta * 0.5 * c_r

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
		for mol in self.molecules:
			
			w = self.width/self.scale/2
			h = self.height/self.scale/2

			if mol.x + mol.r > w:
				mol.check_clip(axis=[1,0],boundary=w)
			elif mol.x - mol.r < -w:
				mol.check_clip(axis=[-1,0],boundary=-w)

			if mol.y + mol.r > h:
				mol.check_clip(axis=[0,1],boundary=h)
			elif mol.y - mol.r < -h:
				mol.check_clip(axis=[0,-1],boundary=-h)

PI = np.arccos(-1)
TwoPI = 2.0*Math.PI
UP = 1.5*PI

def vec_angle(v,d=None):
	d = d or np.linalg.norm(v)

	if d == 0:
		return UP

	v = v / d

	if v[0] == 0:
		if v[1] > 0:
			v_angle = PI/2
		else:
			v_angle = -PI/2
	else:
		v_angle = np.arccos(v[0])
		if v[1] < 0:
			v_angle = TwoPI - v_angle

	return v_angle%(TwoPI)

def set_running():
	document.getElementById("py-status").innerHTML = ''

def draw_circle(x,y,r,c):
	x *= scale
	y *= scale
	r *= scale
	ctx.beginPath()
	ctx.arc(x,y,r,0,TwoPI)
	ctx.fillStyle = c
	ctx.fill()

def draw_circle_stroke(x,y,r,c,w=1):
	x *= scale
	y *= scale
	r *= scale
	ctx.beginPath()
	ctx.lineWidth = w
	ctx.arc(x,y,r,0,TwoPI)
	ctx.strokeStyle = c
	ctx.stroke()

def draw_line(x1,y1,x2,y2,w=2,c='white'):
	x1 *= scale
	y1 *= scale
	x2 *= scale
	y2 *= scale
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

	# sim.add_molecule(make_H2O(0,0))
	# sim.add_molecule(make_H2O(5,0,0))

	# make_water_grid(sim,7,4)
	
	# makers = [
		# make_H2,
		# make_H2O,
		# make_O2,
		# make_NaCl,
		# make_ClI,
		# make_CO2,
		# make_CO,
		# # make_NH2,
		# make_NH3,
		# make_CH4,
		# make_O3,
		# make_C2H4,
		# make_benzene,
	# ]

	makers = [
		make_H2O,
		make_H2O,
		make_C2H6O,
		make_H2O,
		make_H2O,
		make_benzene,
		make_H2O,
		make_H2O,
		make_C2H4,
		make_H2O,
		make_H2O,
		make_benzene,
	]

	make_molecule_grid(sim,makers,5,3,6)
	
	sim.draw()
	
	if run:
		clear_screen()
		# sim.solve(dt=dt,substeps=substeps,random_acc=random_acceleration>0)
		sim.solve(dt=dt,substeps=substeps,random_v=initial_velocity)
		sim.draw()

	draw_loop_proxy = pyodide.ffi.create_proxy(draw_loop)

	if run > 1:
		# interval_id = setInterval(draw_loop_proxy,50)
		interval_id = setInterval(draw_loop_proxy,20)

	if debug:
		print('finished.')
	
counter = 0

def make_molecule_grid(sim,makers,c,r=None,d=4):

	# print((canvas.width/2)/scale)
	# print((canvas.height/2)/scale)

	r = r or math.ceil(len(makers)/c)

	for i,maker in enumerate(makers):
		x = d*(i%c) - d*(c//2)
		y = d*(i//c) - d*(r//2)
		# print(i,x,y)
		sim.add_molecule(maker(x,y))

def make_water_grid(sim,c,r,d=4):
	for i in range(c*r):
		x = d*(i%c) - d*c//2
		y = d*(i//c) - d*r//2
		sim.add_molecule(make_H2O(x,y))

def make_H2(x,y,t=0):
	mol = Molecule('hydrogen',x,y,t)
	mol.add_atom(Atom('H',rt=[COVALENT_RADII['H'],PI/2],q=0.0))
	mol.add_atom(Atom('H',rt=[COVALENT_RADII['H'],-PI/2],q=0.0))
	mol.center_atoms()
	return mol

def make_H2O(x,y,t=0):
	mol = Molecule('water',x,y,t)
	mol.add_atom(Atom('O',q=-0.834))
	mol.add_atom(Atom('H',rt=[COVALENT_RADII['O']+COVALENT_RADII['H'],BOND_ANGLES['H O H']/2],q=0.417))
	mol.add_atom(Atom('H',rt=[COVALENT_RADII['O']+COVALENT_RADII['H'],-BOND_ANGLES['H O H']/2],q=0.417))
	mol.center_atoms()
	return mol

def make_O2(x,y,t=0):
	mol = Molecule('oxygen',x,y,t)
	mol.add_atom(Atom('O',rt=[COVALENT_RADII['O'],PI/2],q=0.0))
	mol.add_atom(Atom('O',rt=[COVALENT_RADII['O'],-PI/2],q=0.0))
	mol.center_atoms()
	return mol

def make_NaCl(x,y,t=0):
	mol = Molecule('salt',x,y,t)
	mol.add_atom(Atom('Na',rt=[COVALENT_RADII['Na'],PI/2],q=None))
	mol.add_atom(Atom('Cl',rt=[COVALENT_RADII['Cl'],-PI/2],q=None))
	mol.center_atoms()
	return mol

def make_ClI(x,y,t=0):
	mol = Molecule('iodine monochloride',x,y,t)
	mol.add_atom(Atom('Cl',rt=[COVALENT_RADII['Cl'],PI/2],q=None))
	mol.add_atom(Atom('I',rt=[COVALENT_RADII['I'],-PI/2],q=None))
	mol.center_atoms()
	return mol

def make_CO2(x,y,t=0):
	mol = Molecule('carbon dioxide',x,y,t)
	mol.add_atom(Atom('C',q=0.6))
	mol.add_atom(Atom('O',rt=[COVALENT_RADII['C']+COVALENT_RADII['O'],BOND_ANGLES['O C O']/2],q=-0.3))
	mol.add_atom(Atom('O',rt=[COVALENT_RADII['C']+COVALENT_RADII['O'],-BOND_ANGLES['O C O']/2],q=-0.3))
	mol.center_atoms()
	return mol

def make_CO(x,y,t=0):
	mol = Molecule('carbon monoxide',x,y,t)
	mol.add_atom(Atom('C',q=0.6))
	mol.add_atom(Atom('O',rt=[bond_length('CO'),PI/2],q=-0.3))
	mol.center_atoms()
	return mol

# no charges!
def make_NH2(x,y,t=0):
	mol = Molecule('amidogen',x,y,t)
	mol.add_atom(Atom('N',q=None))
	mol.add_atom(Atom('H',rt=[bond_length('OH'),BOND_ANGLES['H N H']/2],q=None))
	mol.add_atom(Atom('H',rt=[bond_length('OH'),-BOND_ANGLES['H N H']/2],q=None))
	mol.center_atoms()
	return mol

# no charges!
def make_NH3(x,y,t=0):
	mol = Molecule('ammonia',x,y,t)
	mol.add_atom(Atom('N',q=None))
	mol.add_atom(Atom('H',rt=[bond_length('NH'),0],q=None))
	mol.add_atom(Atom('H',rt=[bond_length('NH'),TwoPI/3],q=None))
	mol.add_atom(Atom('H',rt=[bond_length('NH'),2*TwoPI/3],q=None))
	mol.center_atoms()
	return mol

# no charges!
def make_CH4(x,y,t=0):
	mol = Molecule('methane',x,y,t)
	mol.add_atom(Atom('C',q=None))
	mol.add_atom(Atom('H',rt=[bond_length('HC'),0],q=None))
	mol.add_atom(Atom('H',rt=[bond_length('HC'),PI/2],q=None))
	mol.add_atom(Atom('H',rt=[bond_length('HC'),PI],q=None))
	mol.add_atom(Atom('H',rt=[bond_length('HC'),UP],q=None))
	mol.center_atoms()
	return mol

# no charges!
def make_O3(x,y,t=0):
	mol = Molecule('ozone',x,y,t)
	mol.add_atom(Atom('O',q=None))
	mol.add_atom(Atom('O',rt=[bond_length('OO'),PI/3],q=None))
	mol.add_atom(Atom('O',rt=[bond_length('OO'),-PI/3],q=None))
	mol.center_atoms()
	return mol

# wobbly
def make_C2H4(x,y,t=0):
	mol = Molecule('ethylene',x,y,t)
	mol.add_atom(Atom('C',rt=[COVALENT_RADII['C'],PI/2],q=-0.42))
	mol.add_atom(Atom('C',rt=[COVALENT_RADII['C'],-PI/2],q=-0.42))

	mol.add_atom(Atom('H',xy=[mol.atoms[1].x+bond_length('CH')*np.sin(BOND_ANGLES['C C H']),mol.atoms[1].y+bond_length('CH')*np.cos(BOND_ANGLES['C C H'])],q=0.21))
	mol.add_atom(Atom('H',xy=[mol.atoms[1].x+bond_length('CH')*np.sin(-BOND_ANGLES['C C H']),mol.atoms[1].y+bond_length('CH')*np.cos(-BOND_ANGLES['C C H'])],q=0.21))

	mol.add_atom(Atom('H',xy=[mol.atoms[0].x-bond_length('CH')*np.sin(BOND_ANGLES['C C H']),mol.atoms[0].y-bond_length('CH')*np.cos(BOND_ANGLES['C C H'])],q=0.21))
	mol.add_atom(Atom('H',xy=[mol.atoms[0].x-bond_length('CH')*np.sin(-BOND_ANGLES['C C H']),mol.atoms[0].y-bond_length('CH')*np.cos(-BOND_ANGLES['C C H'])],q=0.21))
	mol.center_atoms()
	return mol

def make_C2H6O(x,y,t=0):
	mol = Molecule('ethanol',x,y,t)
	mol.add_atom(Atom('C',rt=[bond_length('CC'),BOND_ANGLES['C C O']/2],q=None))
	mol.add_atom(Atom('C',q=None))
	mol.add_atom(Atom('O',rt=[bond_length('CO'),-BOND_ANGLES['C C O']/2],q=None))

	mol.add_atom(Atom('H',xy=[mol.atoms[0].x+bond_length('CH')*np.sin(55/180*PI),mol.atoms[0].y+bond_length('CH')*np.cos(55/180*PI)],q=0.21))
	mol.add_atom(Atom('H',xy=[mol.atoms[0].x+bond_length('CH')*np.sin(-55/180*PI),mol.atoms[0].y+bond_length('CH')*np.cos(-55/180*PI)],q=0.21))
	mol.add_atom(Atom('H',xy=[mol.atoms[0].x,mol.atoms[0].y+bond_length('CH')],q=0.21))

	mol.add_atom(Atom('H',xy=[mol.atoms[1].x-bond_length('CH')*np.cos(30/180*PI),mol.atoms[1].y-bond_length('CH')*np.sin(30/180*PI)],q=0.21))
	mol.add_atom(Atom('H',xy=[mol.atoms[1].x-bond_length('CH')*np.cos(-30/180*PI),mol.atoms[1].y-bond_length('CH')*np.sin(-30/180*PI)],q=0.21))
	
	mol.add_atom(Atom('H',xy=[mol.atoms[2].x+bond_length('OH')*np.cos(-30/180*PI),mol.atoms[2].y+bond_length('OH')*np.sin(-30/180*PI)],q=0.21))

	mol.center_atoms()
	return mol

# very unstable
def make_benzene(x,y,t=0):

	mol = Molecule('benzene',x,y,t)
	r = COVALENT_RADII['C']/np.sin(PI/6)

	for i in range(6):
		mol.add_atom(Atom('C',rt=[r,i*PI/3],q=None))
		mol.add_atom(Atom('H',rt=[r+bond_length('CH'),i*PI/3],q=None))

	return mol

def clear_screen():
	ctx.clearRect(-canvas.width/2, -canvas.height/2, canvas.width, canvas.height)

def draw_loop():
	global ctx, counter, t

	start = time.perf_counter

	clear_screen()
	if draw_time:
		ctx.fillStyle = 'white'
		ctx.font = "30px Arial"
		# ctx.fillText(counter, -canvas.width/2+10, 50)
		ctx.fillText(f'{t=:.3f}', -canvas.width/2+10, -canvas.height/2+50)

	if draw_axes_lines:
		draw_axes(ctx,canvas)

	sim.solve(dt=dt,substeps=substeps)
	sim.draw()

	t += dt
	counter += 1
	
	if counter == max_steps:
		clearInterval(interval_id)

def cw_normal(vec):
	return np.array([vec[1],-vec[0]])

def random_velocity(strength):
	theta = random.uniform(0,TwoPI)
	return strength*np.array([np.cos(theta),np.sin(theta)])

def bond_length(symbols):
	return COVALENT_RADII[symbols[0]]+COVALENT_RADII[symbols[1]]

COVALENT_RADII = {
'X': 0.2, 'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76, 'Sc': 1.7, 'Ti': 1.6, 'V': 1.53, 'Cr': 1.39, 'Mn': 1.39, 'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 'Ga': 1.22, 'Ge': 1.2, 'As': 1.19, 'Se': 1.2, 'Br': 1.2, 'Kr': 1.16, 'Rb': 2.2, 'Sr': 1.95, 'Y': 1.9, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.4, 'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89, 'Tm': 1.9, 'Yb': 1.87, 'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.7, 'W': 1.62, 'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.4, 'At': 1.5, 'Rn': 1.5, 'Fr': 2.6, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06, 'Pa': 2.0, 'U': 1.96, 'Np': 1.9, 'Pu': 1.87, 'Am': 1.8, 'Cm': 1.69, 'Bk': 0.2, 'Cf': 0.2, 'Es': 0.2, 'Fm': 0.2, 'Md': 0.2, 'No': 0.2, 'Lr': 0.2, 'Rf': 0.2, 'Db': 0.2, 'Sg': 0.2, 'Bh': 0.2, 'Hs': 0.2, 'Mt': 0.2, 'Ds': 0.2, 'Rg': 0.2, 'Cn': 0.2, 'Nh': 0.2, 'Fl': 0.2, 'Mc': 0.2, 'Lv': 0.2, 'Ts': 0.2, 'Og': 0.2
}

MASSES = {
'X': 1.0, 'H': 1.008, 'He': 4.002602, 'Li': 6.94, 'Be': 9.0121831, 'B': 10.81, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998403163, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815385, 'Si': 28.085, 'P': 30.973761998, 'S': 32.06, 'Cl': 35.45, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078, 'Sc': 44.955908, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938044, 'Fe': 55.845, 'Co': 58.933194, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.723, 'Ge': 72.63, 'As': 74.921595, 'Se': 78.971, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90584, 'Zr': 91.224, 'Nb': 92.90637, 'Mo': 95.95, 'Tc': 97.90721, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42, 'Ag': 107.8682, 'Cd': 112.414, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.6, 'I': 126.90447, 'Xe': 131.293, 'Cs': 132.90545196, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116, 'Pr': 140.90766, 'Nd': 144.242, 'Pm': 144.91276, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93033, 'Er': 167.259, 'Tm': 168.93422, 'Yb': 173.054, 'Lu': 174.9668, 'Hf': 178.49, 'Ta': 180.94788, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217, 'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.592, 'Tl': 204.38, 'Pb': 207.2, 'Bi': 208.9804, 'Po': 208.98243, 'At': 209.98715, 'Rn': 222.01758, 'Fr': 223.01974, 'Ra': 226.02541, 'Ac': 227.02775, 'Th': 232.0377, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.04817, 'Pu': 244.06421, 'Am': 243.06138, 'Cm': 247.07035, 'Bk': 247.07031, 'Cf': 251.07959, 'Es': 252.083, 'Fm': 257.09511, 'Md': 258.09843, 'No': 259.101, 'Lr': 262.11, 'Rf': 267.122, 'Db': 268.126, 'Sg': 271.134, 'Bh': 270.133, 'Hs': 269.1338, 'Mt': 278.156, 'Ds': 281.165, 'Rg': 281.166, 'Cn': 285.177, 'Nh': 286.182, 'Fl': 289.19, 'Mc': 289.194, 'Lv': 293.204, 'Ts': 293.208, 'Og': 294.214
}

BOND_ANGLES = {
'H N C': 111.0, 'N C C': 110.0, 'C C O': 110.5, 'C C C': 108.0, 'N C H': 109.5, 'C N C': 112.0, 'H C C': 109.5, 'C O C': 109.6, 'C S C': 95.0, 'H N H': 120.0, 'H O C': 108.0, 'H C H': 109.0, 'H S C': 95.0, 'N C N': 120.0, 'O C C': 118.0, 'O C H': 121.7, 'O C N': 122.5, 'O C O': 124.0, 'S C C': 112.5, 'S C H': 111.3, 'S S C': 103.3, 'C C H': 110.1, 'C C N': 122.0, 'H O H': 104.5, 'C C S': 124.0, 'N C S': 116.4, 'N C O': 124.0, 'S C S': 124.0, 'O C S': 125.0, 'C C Cl': 120.0, 'C C Br': 120.0, 'C C I': 120.0, 'N C Br': 120.0, 'C C F': 118.8, 'Cl C Cl': 109.0, 'Br C Br': 110.5, 'F C F': 107.0, 'Cl C H': 108.5, 'Br C H': 107.0, 'C C P': 117.0, 'P C F': 122.0, 'F C H': 108.9, 'P C H': 110.0, 'C N N': 115.0, 'C N H': 113.0, 'C N O': 116.0, 'O N O': 128.0, 'C N S': 111.0, 'N N N': 102.2, 'N N O': 103.0, 'N N H': 119.5, 'C N P': 118.3, 'P N H': 123.6, 'S N H': 113.1, 'C O N': 108.5, 'C O P': 120.0, 'C O S': 108.0, 'P O P': 143.0, 'C O H': 115.0, 'N O H': 101.5, 'P O H': 115.0, 'O P O': 111.6, 'O P S': 131.0, 'S P S': 128.5, 'C P O': 94.0, 'N P O': 110.6, 'C S N': 103.0, 'C S S': 103.3, 'C S H': 95.0, 'C S O': 98.0, 'N S O': 94.2, 'O S O': 109.5, 'N S N': 102.3, 'F Al F': 109.5, 'H C O': 108.9, 'H O P': 115.0, 'O S N': 103.0, 'S N C': 113.0, 'H C N': 113.5, 'C S Fe': 100.6, 'S Fe N': 90.0, 'Fe N C': 128.1, 'H C Fe': 180.0, 'N Fe C': 90.0, 'N Fe N': 90.0, 'O C Fe': 180.0, 'O Fe N': 90.0, 'O O Fe': 180.0, 'H C P': 110.0, 'S Zn S': 111.8, 'C S Zn': 95.0, 'C N Zn': 120.7, 'S Zn N': 108.1, 'N Zn N': 107.8, 'F C C': 118.8, 'H N P': 123.6, 'N C Se': 122.5, 'C O O': 104.0, 'O O H': 98.3, 'O Si O': 117.0, 'Si O Si': 150.5, 'H Si H': 119.0, 'H Si O': 118.0, 'Si O H': 122.5, 'Al O Al': 98.0, 'Al O Si': 117.0, 'H O Al': 93.4, 'O Al H': 93.4, 'O Al O': 90.0, 'O Si H': 118.0, 'O S C': 99.0
}

for key,angle in BOND_ANGLES.items():
	BOND_ANGLES[key] = PI*angle/180.0

COLOURS = {
	'H':'Azure',
	'O':'Crimson',
	'C':'LightSlateGrey',
	'N':'RoyalBlue',
	'Cl':'Mediumspringgreen',
	'Na':'Mediumpurple',
	'I':'Purple',
}

main()
