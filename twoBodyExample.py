from vpython import *

G = 6.67e-11
m1 = 3e30
m2 = 3e30
rdist = 1.5e11

M = m1 + m2
# Let the initial position of the objects be on the x-axis
x1 = - m2 * rdist / M 
x2 = m1 * rdist / M

star1 = sphere(pos = vector(x1, 0, 0), radius = 8e9,
	color = color.cyan, make_trail = True)
star2 = sphere(pos = vector(x2, 0, 0), radius = 8e9, 
	color = color.yellow, make_trail = True)

r = star2.pos - star1.pos

# Center of mass (must be <0, 0, 0>)
R = (m1 * star1.pos + m2 * star2.pos) / m1
center = sphere(pos = R, radius = 2e9,
	color = color.red)

# Velocity for circle motion
v1circle = sqrt(G * m2 * mag(star1.pos) / mag(r) ** 2)
print(v1circle)

# We are working in xy plane
star1.v = vector(0, 0.7 * v1circle, 0)

# Momentum
star1.p = m1 * star1.v
star2.p = -star1.p

# Reduced mass
mu = m1 * m2 / M

# Magnitude of the linear momentum
l = mag(cross(star1.pos, star1.p) +
	cross(star2.pos, star2.p))

t = 0
# Delta time
dt = 1000

# Graph of potential energy
tgraph = graph(xtitle = "r [m]", ytitle = "U [J]")
fU = gcurve(color = color.blue, dot = True)
fUg = gcurve(color = color.green, dot = True)
fUc = gcurve(color = color.purple, dot = True)

while t < 1e8:
	# Equivalent for Sleep(1 / 1000) (1000 operations per second)
	rate(1000)
	r = star2.pos - star1.pos
	F2 = -G * m1 * m2 * norm(r) / mag(r) ** 2
	# Change of potential due to gravitational force (F = ma; p = mv; p + F * dt = m(v + dt * a))
	star2.p += F2 * dt
	star1.p -= F2 * dt
	# Change of position due to change of speed (p = mv; pos += v * dt; pos += p * dt / m)
	star1.pos += star1.p * dt / m1
	star2.pos +=  star2.p * dt / m2

	# Showing the potential energy Ueff = Ugravitational + Ucentrifugal
	# Ug is antiderivative from Fgravitational
	# Formula for Uc is too hard to explain here
	Ug = -G * m1 * m2 / mag(r)
	Uc = l ** 2 / (2 * mu * mag(r) ** 2)
	Uef = Ug + Uc
	
	fU.plot(mag(r), Uef)
	fUg.plot(mag(r), Ug)
	fUc.plot(mag(r), Uc)
	t += dt
