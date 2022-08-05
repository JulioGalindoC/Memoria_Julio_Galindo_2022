from numpy import cross

def calcular_volumen(x,y,z,debug=False):
	ni, nj = x.shape
	printif = lambda s: print(s) if debug else None
	
	printif(f"ni = {ni} nj = {nj}")

	V = 0.
	for i in range(ni-1):
		for j in range(nj-1):
			I = i+1
			J = j+1
			x0, y0, z0 = x[i,j], y[i,j], z[i,j]
			x1, y1, z1 = x[I,j], y[I,j], z[I,j]
			x2, y2, z2 = x[I,J], y[I,J], z[I,J]
			x3, y3, z3 = x[i,J], y[i,J], z[i,J]
			printif(f"x0, y0, z0 = {x0}, {y0}, {z0}")
			printif(f"x1, y1, z1 = {x1}, {y1}, {z1}")
			printif(f"x2, y2, z2 = {x2}, {y2}, {z2}")
			printif(f"x3, y3, z3 = {x3}, {y3}, {z3}")

			dx = x1 - x0

			a = [0,y0,z0]
			b = [0,y3,z3]
			dA1 = cross(a,b)/2

			c = [0,y1,z1]
			d = [0,y2,z2]
			dA2 = cross(c,d)/2

			printif(f"dx = {dx}")
			printif(f"dA1 = {dA1}")
			printif(f"dA2 = {dA2}")

			dV = dx * (dA1[0] + dA2[0]) / 2

			printif(f"dV = {dV}")

			V += dV
	return V

if __name__ == "__main__":

	#Probar con un cilindro...
	from numpy import meshgrid, linspace, pi, cos, sin
	r = 1.
	h = 4.
	th = linspace(0, 2*pi, 10)
	x = linspace(0, h, 3)

	TH,X = meshgrid(th,x)

	Y = r*cos(TH)
	Z = r*sin(TH)


	RealVol = pi*r**2*h

	ApproxVol = calcular_volumen(X,Y,Z, debug=True)

	print(f"RealVol={RealVol}")
	print(f"ApproxVol={ApproxVol}")