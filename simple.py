

"""

python3 -m http.server
localhost:8000/test.html

"""

debug = False

from datetime import datetime
# import time
# START = time.perf_counter()

from js import document, Math, setInterval, clearInterval, setTimeout
# import pyodide
# from pyodide import create_proxy
from pyodide.ffi import create_proxy
# import numpy as np

def set_running():
	document.getElementById("py-status").innerHTML = ''

def draw_circle(ctx,x,y,r,c='white'):
	ctx.beginPath()
	ctx.arc(x,y,r,0,2*Math.PI)
	ctx.fillStyle = c
	ctx.fill()

# def draw_line(ctx,x1,y1,x2,y2,w=2,c='white'):
# 	ctx.beginPath()
# 	ctx.lineWidth = w
# 	ctx.moveTo(x1,y1)
# 	ctx.lineTo(x2-x1,y2-y1)
# 	ctx.strokeStyle = c
# 	ctx.stroke()

# def draw_axes(ctx,canvas):
# 	ctx.lineWidth = 2
# 	draw_line(ctx,0,0,0,-canvas.height/2)
# 	draw_line(ctx,0,0,0,canvas.height/2)
# 	draw_line(ctx,0,0,canvas.width/2,0)
# 	draw_line(ctx,0,0,-canvas.width/2,0)

def main():

	global canvas,ctx

	set_running()
	
	canvas = document.getElementById("headerCanvas")
	ctx = canvas.getContext("2d")

	ctx.translate(canvas.width/2,canvas.height/2)

	print('create proxy')
	# draw_loop_proxy = pyodide.ffi.create_proxy(draw_loop)
	draw_loop_proxy = create_proxy(draw_loop)
	interval_id = setInterval(draw_loop_proxy,50)
	
counter = 0
def draw_loop():
	global counter, ctx, canvas
	counter += 1
	ctx.clearRect(-canvas.width/2, -canvas.height/2, canvas.width, canvas.height);
	ctx.font = "30px Arial";
	ctx.fillText(datetime.now(), 0, 0);
	# ctx.fillText(counter, 0, 50);
	ctx.fillText((2*counter)%(canvas.height/2), 0, 50);
	# print(counter,counter%(canvas.height/2))
	# print(counter,datetime.now())

	### THIS IS BREAKING!!!
	draw_circle(ctx,-200,(2*counter)%(canvas.height/2),50)

COVALENT_RADII = {'X': 0.2, 'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76, 'Sc': 1.7, 'Ti': 1.6, 'V': 1.53, 'Cr': 1.39, 'Mn': 1.39, 'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 'Ga': 1.22, 'Ge': 1.2, 'As': 1.19, 'Se': 1.2, 'Br': 1.2, 'Kr': 1.16, 'Rb': 2.2, 'Sr': 1.95, 'Y': 1.9, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.4, 'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89, 'Tm': 1.9, 'Yb': 1.87, 'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.7, 'W': 1.62, 'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.4, 'At': 1.5, 'Rn': 1.5, 'Fr': 2.6, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06, 'Pa': 2.0, 'U': 1.96, 'Np': 1.9, 'Pu': 1.87, 'Am': 1.8, 'Cm': 1.69, 'Bk': 0.2, 'Cf': 0.2, 'Es': 0.2, 'Fm': 0.2, 'Md': 0.2, 'No': 0.2, 'Lr': 0.2, 'Rf': 0.2, 'Db': 0.2, 'Sg': 0.2, 'Bh': 0.2, 'Hs': 0.2, 'Mt': 0.2, 'Ds': 0.2, 'Rg': 0.2, 'Cn': 0.2, 'Nh': 0.2, 'Fl': 0.2, 'Mc': 0.2, 'Lv': 0.2, 'Ts': 0.2, 'Og': 0.2}

COLOURS = {'H':'Azure','O':'Crimson'}

main()
