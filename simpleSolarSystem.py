from vpython import *

G = 6.67e-11

class Body :
	def __init__(self, mass, sphere):
		self.mass = mass
		self.sphere = sphere
		self.v = vector(0, 0, 0)
		self.p = 0

class TwoBodySystem :
	def __init__(self, body1, body2, distance, centerRadius, curvature, R = vector(0, 0, 0)):
		global G

		self.R = vector(0, 0, 0)

		self.body1 = body1
		self.body2 = body2
		self.distance = distance

		self.M1 = body1.mass #type(body1) == Body if body1.mass else body1.mu
		self.M2 = body2.mass #type(body2) == Body if body2.mass else body2.mu

		# Mass
		M = self.M1 + self.M2
		self.mu = self.M1 * self.M2 / M

		# Start X
		x1 = - self.M2 * distance / M
		x2 = self.M1 * distance / M

		# Positions
		self.body1.sphere.pos = R + vector(x1, 0, 0);
		self.body2.sphere.pos = R + vector(x2, 0, 0);
		self.pos1 = vector(x1, 0, 0)
		self.pos2 = vector(x2, 0, 0)

		self.r = self.pos2 - self.pos1

		self.center = sphere(pos = R,
			radius = centerRadius, color = color.red, make_trail = True)

		# Start velocity
		v1 = sqrt(G * self.M2 * mag(self.pos1) / mag(self.r) ** 2)
		self.body1.v = vector(0, curvature * v1, 0)

		# Momentum
		self.body1.p = self.M1 * self.body1.v
		self.body2.p = -self.body1.p

	def iteration(self, dt, R = vector(0, 0, 0)):
		global G

		self.r = self.pos2 - self.pos1
		F2 = -G * self.M1 * self.M2 * norm(self.r) / mag(self.r) ** 2
		self.body2.p += F2 * dt
		self.body1.p -= F2 * dt
		
		self.pos1 += self.body1.p * dt / self.M1
		self.pos2 += self.body2.p * dt / self.M2

		self.body1.sphere.pos = R + self.pos1
		self.body2.sphere.pos = R + self.pos2
		self.center.pos = R;
		pass

class SolarPlanet :
	def __init__(self, dist, mass, curv, colr):
		self.mass = mass
		self.dist = dist
		self.curv = curv
		self.colr = colr

AU = 1.49e11
Me = 6.0e24

SolarSystem = [
	SolarPlanet(0.588 * AU, 0.055 * Me, 0.795, color.white * 0.3), # Mercury
	SolarPlanet(0.718 * AU, 0.815 * Me, 0.993, color.orange * 0.5), # Venus
	SolarPlanet(1.017 * AU, 1.000 * Me, 0.984, color.cyan), # Earth
	SolarPlanet(1.666 * AU, 0.107 * Me, 0.907, color.red * 0.3), # Mars

	SolarPlanet(5.457 * AU, 318.0 * Me, 0.951, color.orange * 0.3), # Jupiter
	SolarPlanet(10.07 * AU, 95.00 * Me, 0.943, color.yellow * 0.5), # Saturn
	SolarPlanet(20.06 * AU, 14.00 * Me, 0.953, vector(0.7, 0.7, 1)), # Uranus
	SolarPlanet(30.47 * AU, 17.00 * Me, 0.992, vector(0.5, 0.5, 1)), # Neptune
]

SunM = 2.0e30

TwoBodies = []
for planet in SolarSystem :
	TwoBodies.append(
		TwoBodySystem(Body(SunM, sphere(radius = 3.0e10, color = color.yellow, emissive = True,
										make_trail = True, shininess = 1.0)),
				Body(planet.mass, sphere(radius = 1.5e10 * ((planet.mass / Me) ** 0.2), 
										color = planet.colr, make_trail = True)),
				planet.dist,
				1,
				planet.curv)
	)

MoonM = 7.3e22
MoonD = 405.4e6

#EarthMoon = TwoBodySystem(Body(EarthM, sphere(radius = 5.0e9, color = color.cyan, make_trail = True)),
#					Body(MoonM, sphere(radius = 0.7e8, color = color.white, make_trail = True)),
#					MoonD,
#					1e7,
#					0.945,
#					SunEarth.body2.sphere.pos)

t = 0
dt = 1200

while t < 1e10:
	rate(72 * 30)
	for System in TwoBodies :
		System.iteration(dt)

	# SunEarth.iteration(dt)
	# EarthMoon.iteration(dt, SunEarth.body2.sphere.pos);
	t += dt