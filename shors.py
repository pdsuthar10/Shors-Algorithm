#!/usr/bin/env python


import math
import random
import argparse


def printNone(str):
	pass

def printVerbose(str):
	print(str)

printInfo = printNone

####################################################################################################                                                                                                  
                                        Quantum Components                                         
####################################################################################################

class Mapping:
	def init(self, state, ampt):
		self.state = state
		self.ampt = ampt


class QuantumState:
	def init(self, ampt, reg):
		self.ampt = ampt
		self.reg = reg
		self.entangled = {}

	def entangle(self, fromState, ampt):
		reg = fromState.reg
		entanglement = Mapping(fromState, ampt)
		try:
			self.entangled[reg].append(entanglement)
		except KeyError:
			self.entangled[reg] = [entanglement]

	def entangles(self, reg = None):
		entangles = 0
		if reg is None:
			for states in self.entangled.values():
				entangles += len(states)
		else:
			entangles = len(self.entangled[reg])

		return entangles


class Qubitreg:
	def init(self, numBits):
		self.numBits = numBits
		self.numStates = 1 << numBits
		self.entangled = []
		self.states = [QuantumState(complex(0.0), self) for x in range(self.numStates)]
		self.states[0].ampt = complex(1.0)

	def propagate(self, fromreg = None):
		if fromreg is not None:
			for state in self.states:
				ampt = complex(0.0)

				try:
					entangles = state.entangled[fromreg]
					for entangle in entangles:
						ampt += entangle.state.ampt * entangle.ampt

					state.ampt = ampt
				except KeyError:
					state.ampt = ampt

		for reg in self.entangled:
			if reg is fromreg:
				continue

			reg.propagate(self)

	# Map will convert any mapping to a unitary tensor given each element v
	# returned by the mapping has the property v * v.conjugate() = 1
	#
	def map(self, toreg, mapping, propagate = True):
		self.entangled.append(toreg)
		toreg.entangled.append(self)

		# Create the covariant/contravariant representations
		mapTensorX = {}
		mapTensorY = {}
		for x in range(self.numStates):
			mapTensorX[x] = {}
			codomain = mapping(x)
			for element in codomain:
				y = element.state
				mapTensorX[x][y] = element

				try:
					mapTensorY[y][x] = element
				except KeyError:
					mapTensorY[y] = { x: element }

		# Normalize the mapping:
		def normalize(tensor, p = False):
			lSqrt = math.sqrt
			for vectors in tensor.values():
				sumProb = 0.0
				for element in vectors.values():
					ampt = element.ampt
					sumProb += (ampt * ampt.conjugate()).real

				normalized = lSqrt(sumProb)
				for element in vectors.values():
					element.ampt = element.ampt / normalized

		normalize(mapTensorX)
		normalize(mapTensorY, True)

		# Entangle the regs
		for x, yStates in mapTensorX.items():
			for y, element in yStates.items():
				ampt = element.ampt
				toState = toreg.states[int(y)]
				fromState = self.states[x]
				toState.entangle(fromState, ampt)
				fromState.entangle(toState, ampt.conjugate())

		if propagate:
			toreg.propagate(self)

	def measure(self):
		measure = random.random()
		sumProb = 0.0

		# Pick a state
		finalX = None
		finalState = None
		for x, state in enumerate(self.states):
			ampt = state.ampt
			sumProb += (ampt * ampt.conjugate()).real

			if sumProb > measure:
				finalState = state
				finalX = x
				break

		# If state was found, update the system
		if finalState is not None:
			for state in self.states:
				state.ampt = complex(0.0)

			finalState.ampt = complex(1.0)
			self.propagate()

		return finalX

	def entangles(self, reg = None):
		entangles = 0
		for state in self.states:
			entangles += state.entangles(None)

		return entangles

	def ampts(self):
		ampts = []
		for state in self.states:
			ampts.append(state.ampt)

		return ampts

def printEntangles(reg):
	printInfo("Entagles: " + str(reg.entangles()))

def printampts(reg):
	ampts = reg.ampts()
	for x, ampt in enumerate(ampts):
		printInfo('State #' + str(x) + '\'s ampt: ' + str(ampt))

def hadamard(x, Q):
	codomain = []
	for y in range(Q):
		ampt = complex(pow(-1.0, bitCount(x & y) & 1))
		codomain.append(Mapping(y, ampt))

	return  codomain

# Quantum Modular Exponentiation
def qModExp(a, exp, mod):
	state = modExp(a, exp, mod)
	ampt = complex(1.0)
	return [Mapping(state, ampt)]

# Quantum Fourier Transform
def qft(x, Q):
	fQ = float(Q)
	k = -2.0 * math.pi
	codomain = []

	for y in range(Q):
		theta = (k * float((x * y) % Q)) / fQ
		ampt = complex(math.cos(theta), math.sin(theta))
		codomain.append(Mapping(y, ampt))

	return codomain

def findPeriod(a, N):
	nNumBits = N.bit_length()
	inputNumBits = (2 * nNumBits) - 1
	inputNumBits += 2 if ((1 << inputNumBits) < (N * N)) else 0
	Q = 1 << inputNumBits

	printInfo("Finding the period...")
	printInfo("Q = " + str(Q) + "\ta = " + str(a))
	
	inputreg = Qubitreg(inputNumBits)
	hmdInputreg = Qubitreg(inputNumBits)
	qftInputreg = Qubitreg(inputNumBits)
	outputreg = Qubitreg(inputNumBits)

	printInfo("regs generated")
	printInfo("Performing Hadamard on input reg")

	inputreg.map(hmdInputreg, lambda x: hadamard(x, Q), False)
	#inputreg.hadamard(False)

	printInfo("Hadamard complete")
	printInfo("Mapping input reg to output reg, where f(x) is a^x mod N")

	hmdInputreg.map(outputreg, lambda x: qModExp(a, x, N), False)

	printInfo("Modular exponentiation complete")
	printInfo("Performing quantum Fourier transform on output reg")

	hmdInputreg.map(qftInputreg, lambda x: qft(x, Q), False)
	inputreg.propagate()

	printInfo("Quantum Fourier transform complete")
	printInfo("Performing a measurement on the output reg")

	y = outputreg.measure()

	printInfo("Output reg measured\ty = " + str(y))

	

	printInfo("Performing a measurement on the periodicity reg")

	x = qftInputreg.measure()

	printInfo("QFT reg measured\tx = " + str(x))

	if x is None:
		return None

	printInfo("Finding the period via continued fractions")

	r = cf(x, Q, N)

	printInfo("Candidate period\tr = " + str(r))

	return r

####################################################################################################                                                                                                  
                                        Classical Components                                                                                                                                       
####################################################################################################

BIT_LIMIT = 12

def bitCount(x):
	sumBits = 0
	while x > 0:
		sumBits += x & 1
		x >>= 1

	return sumBits

# Greatest Common Divisor
def gcd(a, b):
	while b != 0:
		tA = a % b
		a = b
		b = tA

	return a

# Extended Euclidean
def extendedGCD(a, b):
	frac = []
	while b != 0:
		frac.append(a // b)
		tA = a % b
		a = b
		b = tA

	return frac

# Continued frac
def cf(y, Q, N):
	frac = extendedGCD(y, Q)
	depth = 2

	def partial(frac, depth):
		c = 0
		r = 1

		for i in reversed(range(depth)):
			tR = frac[i] * r + c
			c = r
			r = tR

		return c

	r = 0
	for d in range(depth, len(frac) + 1):
		tR = partial(frac, d)
		if tR == r or tR >= N:
			return r

		r = tR

	return r

# Modular Exponentiation
def modExp(a, exp, mod):
	fx = 1
	exp = int(exp)
	while exp > 0:
		if (exp & 1) == 1:
			fx = fx * a % mod
		a = (a * a) % mod
		exp = exp >> 1

	return fx

def pick(N):
	a = math.floor((random.random() * (N - 1)) + 0.5)
	return a

def checkCandidates(a, r, N, nh):
	if r is None:
		return None

	# Check multiples
	for k in range(1, nh + 2):
		tR = k * r
		if modExp(a, a, N) == modExp(a, int(a + tR), N):
			return (tR)

	# Check lower nh
	for tR in range(r - nh, r):
		if modExp(a, a, N) == modExp(a, a + tR, N):
			return tR

	# Check upper neigborhood
	for tR in range(r + 1, r + nh + 1):
		if modExp(a, a, N) == modExp(a, a + tR, N):
			return tR

	return None

def shors(N, atts = 1, nh = 0.0, numprd = 1):
	if(N.bit_length() > BIT_LIMIT or N < 3):
		return False

	prd = []
	nh = math.floor(N * nh) + 1

	printInfo("N = " + str(N))
	printInfo("Neighbourhood = " + str(nh))
	printInfo("Number of periods = " + str(numprd))

	for attempt in range(atts):
		printInfo("\nAttempt #" + str(attempt))

		a = pick(N)
		while a < 2:
			a = pick(N)

		d = gcd(a, N)
		if d > 1:
			printInfo("Found factors classically, re-attempt")
			printInfo("Factor: " + str(a))
			continue

		r = findPeriod(a, N)

		printInfo("Checking candidate period, nearby values, and multiples")

		r = checkCandidates(a, r, N, int(nh))

		if r is None:
			printInfo("Period was not found, re-attempt")
			continue

		if (r % 2) > 0:
			printInfo("Period was odd, re-attempt")
			continue

		d = modExp(a, (r // 2), N)
		if r == 0 or d == (N - 1):
			printInfo("Period was trivial, re-attempt")
			continue

		printInfo("Period found\tr = " + str(r))

		prd.append(r)
		if(len(prd) < numprd):
			continue

		printInfo("\nFinding least common multiple of all prd")

		r = 1
		for period in prd:
			d = gcd(prd, r)
			r = (r * prd) // d

		b = modExp(a, (r // 2), N)
		f1 = gcd(N, b + 1)
		f2 = gcd(N, b - 1)

		return [f1, f2]

	return None

####################################################################################################                                                                                                 
                                   Command-line functionality                                                                                                                              
####################################################################################################

def parseArgs():
	parser = argparse.ArgumentParser(description='Simulate Shor\'s algorithm for N.')
	parser.add_argument('-a', '--atts', type=int, default=20, help='Number of quantum attemtps to perform')
	parser.add_argument('-n', '--nh', type=float, default=0.01, help='nh size for checking candidates (as percentage of N)')
	parser.add_argument('-p', '--prd', type=int, default=2, help='Number of period to get before determining least common multiple')
	parser.add_argument('-v', '--verbose', type=bool, default=True, help='Verbose')
	parser.add_argument('N', type=int, help='The integer to factor')
	return parser.parse_args()

def main():
	args = parseArgs()

	global printInfo
	if args.verbose:
		printInfo = printVerbose
	else:
		printInfo = printNone

	factors = shors(args.N, args.atts, args.nh, args.prd)
	if factors is not None:
		print("Factors:\t" + str(factors[0]) + ", " + str(factors[1]))

if __name__ == "__main__":
	main()