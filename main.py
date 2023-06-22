
"""

python3 -m http.server
localhost:8000/test.html

"""

debug = False

import time
START = time.perf_counter()

if debug:
	print('importing from js... ',end='')	
from js import document, Math, setInterval, clearInterval, setTimeout
import pyodide

if debug:
	print('Done.')	
	print('importing numpy... ',end='')	
import numpy as np
if debug:
	print('Done.')	

class Molecule:

	num_molecules = 0

	def __init__(self):
		self.num_atoms = 0
		self.index = self.num_molecules
		self.num_molecules += 1
		self.atoms = []
		self.angles = []
		self.bonds = []

	def add_atom(self,atom):
		atom.index = self.num_atoms
		self.num_atoms += 1
		self.atoms.append(atom)

	def add_angle(self,angle):
		self.angles.append(angle)

	def add_bond(self,i,j,order=1):
		if debug:
			print(f'Molecule.add_bond({i},{j},{order})')
		atom1 = self.atoms[i]
		atom2 = self.atoms[j]
		bond = Bond(atom1,atom2,order=order)
		self.bonds.append(bond)

class Bond:
	def __init__(self,atom1,atom2,order=1): 
		self.atom1 = atom1
		self.atom2 = atom2
		self.length = COVALENT_RADII[atom1.symbol] + COVALENT_RADII[atom2.symbol]
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
	def __init__(self,atom1,atom2,atom3,angle,strength=100):
		self.atom1 = atom1
		self.atom2 = atom2
		self.atom3 = atom3
		self.angle = angle
		self.strength = strength

	def apply_constraint(self):

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

		# # d = np.linalg.norm(v)
		# # n = v/d

		# delta = self.length - d

		# self.atom1.position += n * delta * 0.5
		# self.atom2.position -= n * delta * 0.5

class Atom:
	def __init__(self,symbol,x=0,y=0,r=1):
		self.symbol = symbol
		self.bonds = []
		self.colour = COLOURS[symbol]
		self.r = COVALENT_RADII[symbol]
		if debug:
			print('Atom.__init__',symbol)
		self.position = np.array([x,y],dtype=float)
		self.old_pos = np.array([x,y],dtype=float)
		print(symbol,x,y)
		self.acceleration = np.array([0.0,0.0],dtype=float)
	@property
	def x(self):
		return self.position[0]
	@property
	def y(self):
		return self.position[1]
	def update(self,dt):
		# print(len(self.position),len(self.old_pos))
		v = self.position - self.old_pos
		self.old_pos = self.position.copy()
		self.position = self.position + v*dt + self.acceleration*dt**2
		self.acceleration = np.array([0.0,0.0])

class Simulation:

	def __init__(self,canvas,ctx,scale=40):
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
			# for bond in molecule.bonds:
			# 	draw_line(self.ctx,bond.atom1.x*self.scale,bond.atom1.y*self.scale,bond.atom2.x*self.scale,bond.atom2.y*self.scale,w=2,c='black')

	def solve(self,dt,substeps=4):
		dt = dt/substeps
		for i in range(substeps):
			# self.apply_gravity([0,10000])
			self.apply_center_gravity(500)
			# self.check_collisions()
			self.apply_angle_constraints()
			self.apply_bond_constraints()
			self.clip()
			self.update_positions(dt)
			pass

	def update_positions(self,dt):
		for molecule in self.molecules:
			for atom in molecule.atoms:
				atom.update(dt)

	def apply_gravity(self,g):
		g = np.array(g)
		for molecule in self.molecules:
			for atom in molecule.atoms:
				atom.acceleration += g

	def apply_center_gravity(self,g):
		g = np.array(g)
		for molecule in self.molecules:
			for atom in molecule.atoms:
				d = np.linalg.norm(atom.position)
				if d > 0:
					atom.acceleration += -g*atom.position/d

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
					atom.position = np.array([w+atom.r,atom.y])
				elif atom.y + atom.r > h:
					atom.position = np.array([atom.x,h-atom.r])
				elif atom.y - atom.r < -h:
					atom.position = np.array([atom.x,h+atom.r])
				pass

def set_running():
	document.getElementById("py-status").innerHTML = ''

def draw_circle(ctx,x,y,r,c):
	ctx.beginPath()
	ctx.arc(x,y,r,0,2*Math.PI)
	ctx.fillStyle = c
	ctx.fill()

def draw_line(ctx,x1,y1,x2,y2,w=2,c='white'):
	ctx.beginPath()
	ctx.lineWidth = w
	ctx.moveTo(x1,y1)
	ctx.lineTo(x2-x1,y2-y1)
	ctx.strokeStyle = c
	ctx.stroke()

def draw_axes(ctx,canvas):
	ctx.lineWidth = 2
	draw_line(ctx,0,0,0,-canvas.height/2)
	draw_line(ctx,0,0,0,canvas.height/2)
	draw_line(ctx,0,0,canvas.width/2,0)
	draw_line(ctx,0,0,-canvas.width/2,0)

def main():

	global interval_id, t, sim, ctx, canvas, dt
	set_running()
	
	canvas = document.getElementById("headerCanvas")
	ctx = canvas.getContext("2d")

	if debug:
		print(f'{canvas.width=}')
		print(f'{canvas.height=}')

	ctx.translate(canvas.width/2,canvas.height/2)

	sim = Simulation(canvas,ctx)

	sim.width = canvas.width
	sim.height = canvas.height

	sim.add_molecule(make_H2O(x=0,y=0))
	sim.add_molecule(make_H2O(x=3,y=3))

	t = 0.0
	dt = 0.01
	
	sim.draw()
	sim.solve(dt)

	draw_loop_proxy = pyodide.ffi.create_proxy(draw_loop)

	interval_id = setInterval(draw_loop_proxy,50)
	
counter = 0

def make_H2O(x,y):
	h2o = Molecule()
	h2o.add_atom(Atom('O',x,y))
	h2o.add_atom(Atom('H',x+1,y))
	h2o.add_atom(Atom('H',x-1,y))
	h2o.add_bond(0,1)
	h2o.add_bond(0,2)
	h2o.add_angle(Angle(h2o.atoms[1],h2o.atoms[0],h2o.atoms[2],100))
	return h2o

def draw_loop():
	global ctx, counter, t

	start = time.perf_counter

	ctx.clearRect(-canvas.width/2, -canvas.height/2, canvas.width, canvas.height);
	ctx.fillStyle = 'white'
	ctx.font = "30px Arial";
	# ctx.fillText(counter, -canvas.width/2+10, 50);
	ctx.fillText(f'{t=:.2f}', -canvas.width/2+10, -canvas.height/2+50);

	draw_axes(ctx,canvas)
	sim.solve(dt)
	sim.draw()

	t += dt
	counter += 1
	
	if counter == 1000:
		clearInterval(interval_id)

COVALENT_RADII = {'X': 0.2, 'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76, 'Sc': 1.7, 'Ti': 1.6, 'V': 1.53, 'Cr': 1.39, 'Mn': 1.39, 'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 'Ga': 1.22, 'Ge': 1.2, 'As': 1.19, 'Se': 1.2, 'Br': 1.2, 'Kr': 1.16, 'Rb': 2.2, 'Sr': 1.95, 'Y': 1.9, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.4, 'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89, 'Tm': 1.9, 'Yb': 1.87, 'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.7, 'W': 1.62, 'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.4, 'At': 1.5, 'Rn': 1.5, 'Fr': 2.6, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06, 'Pa': 2.0, 'U': 1.96, 'Np': 1.9, 'Pu': 1.87, 'Am': 1.8, 'Cm': 1.69, 'Bk': 0.2, 'Cf': 0.2, 'Es': 0.2, 'Fm': 0.2, 'Md': 0.2, 'No': 0.2, 'Lr': 0.2, 'Rf': 0.2, 'Db': 0.2, 'Sg': 0.2, 'Bh': 0.2, 'Hs': 0.2, 'Mt': 0.2, 'Ds': 0.2, 'Rg': 0.2, 'Cn': 0.2, 'Nh': 0.2, 'Fl': 0.2, 'Mc': 0.2, 'Lv': 0.2, 'Ts': 0.2, 'Og': 0.2}

COLOURS = {'H':'Azure','O':'Crimson'}

main()
