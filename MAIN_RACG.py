#!/bin/python

import time
import random
import sys,math

class atom():

	def __init__(self,element,charge,coord):		# KEEP ATOM INFO 
		self.name = element				# NAME (e.g., 'Bi', 'O')
		self.q = charge					# CHARGE
		self.xyz = [coord[0],coord[1],coord[2]]		# POSITION
		self.f = [0.,0.,0.]				# FORCE ACTION ON THE ATOM

class RACG(atom):

	def __init__(self):

		self.bma = 3000.				# DEFAULT CATION - ANION Born-Mayer PAREMETERS
		self.bmr = 0.30					#
		self.lja = -36.					# Lennard-Jones to avoid too sparsed configuration (1/r^6) scaled attractive energy / force; b/t anion - anion
		self.gnorm_tol = 0.035				#
		self.gnorm_div_tol = 15.00			# gnorm diverge tolerance

		self.max_trial = 80000				# MAX TRIAL FOR RANDOM CONFIG GENERATION (Success conditions are, 1. interatomic distance; 2. interatomic force)

		self.box_len = []				# BOX SIZE : a box for placing atoms (randomly)
		self.element = []				# list saves elements
		self.stochio = []				# list saves chemical stochiometry
		self.charge  = []				# list saves physical charges of each elements
		self.size    = None				# var  saves a size of building unit (e.g., (PbO)n -> size = 1, (Bi2O3)n -> size = 5)

		self.inter_atom_dist_tol = 1.5			# default setups : allowed minimum interatomic distance in initial config
		self.inter_atom_force_tol = 10.			# default setups : allowed minimum interatomic forces   in initial config

		# CHECK LIST BEFORE GETTING INTO ACTUAL RANDOM BUILDING
		self.number_of_checklist = 5						# checklist number
		self.check_list = [False for i in range(self.number_of_checklist)]	# checklist board

		self.opti_max = 4096							# max iteration
		self.ls_max_trial = 512							# max trial for line-search sub-problem solver
		self.opti_flag = False							# if use optimisation flag : default = False
		self.to_ev_const = 14.39964390675221758120				# CGS to Atomic Unit Constant




		####################################################################################################################################################
		#	READ INPUT FILE; default = 'RACG.in'
		####################################################################################################################################################
		try:
			with open("RACG.in","r") as f:

				for line in f:
					if line[0] == "#" or line[0] == "\n":	# SKIP IF LINE IS COMMENTED OR EMPTY
						continue
					else:
						spl = line.split()
						stride = len(spl)

						if "BOX_SIZE:" ==  spl[0]:				# GET BOX SIZE
							self.box_len.append( float(spl[1]) )
							self.box_len.append( float(spl[2]) )
							self.box_len.append( float(spl[3]) )
							self.check_list[0] = True
						elif "BASE_ELEMENT:" == spl[0]:				# GET BASE ATOM SPECIES 
							for i in range(stride-1):
								self.element.append(spl[i+1])
							self.check_list[1] = True
						elif "UNIT_STOCHIO:" == spl[0]:				# GET STOCHIOMETRY
							for i in range(stride-1):
								self.stochio.append(int(spl[i+1]))
							self.check_list[2] = True
						elif "ATOM_CHARGE:" == spl[0]:				# GET ATOM CHARGE
							for i in range(stride-1):
								self.charge.append(float(spl[i+1]))
							self.check_list[3] = True
						elif "UNIT_SIZE:" == spl[0]:
							self.size = int(spl[1])
							self.check_list[4] = True

						# SOME DEFAULT PARAMETERS: IF THERE ARE IN CUSTOM MODE

						elif "INTER_ATOMIC_DISTANCE_TOL(Angs):" == spl[0]:	# THERE ARE DEFAULT SETUPS BUT ... USER CAN CUSTOM
							self.inter_atom_dist_tol = float(spl[1])
						elif "INTER_ATOMIC_FORCE_TOL(eV/Angs):" == spl[0]:
							self.inter_atom_force_tol = float(spl[1])
						elif "BM_A(eV):" == spl[0]:
							self.bma = float(spl[1])
						elif "BM_RHO(/Angs):" == spl[0]:
							self.bmr = float(spl[1])
						elif "LJ_(Angs^6):" == spl[0]:
							self.lja = float(spl[1])
							self.lja = -self.lja
						elif "GNORM_TOL(eV/Angs):" == spl[0]:
							self.gnorm_tol = float(spl[1])


				if False in self.check_list:
					print("Error Invoked !!! : essential input parameter(s) is(are) not specified")
					sys.exit()

		except FileNotFoundError:
			print("Error Invoked !!! : input file 'RAGC.in' is not found")
			sys.exit()
	


	####################################################################################################################################################
	#	SOME FUNCTIONS FOR INTERNAL USES : distance, energy, force calculations
	####################################################################################################################################################

	def get_dist(self,A,B):
		dx = A.xyz[0] - B.xyz[0]
		dy = A.xyz[1] - B.xyz[1]
		dz = A.xyz[2] - B.xyz[2]
		return math.sqrt(dx*dx + dy*dy + dz*dz)

	def get_energy(self,A,B):
		r = self.get_dist(A,B)
		ret =  self.to_ev_const*A.q*B.q/r
		if A.q*B.q < 0.:	# if pair charges are not in the same sign > add short_range repulsion term (Born-Mayer)
			rep = self.bma*math.exp(-r/self.bmr)
			ret += rep

		if A.q < 0 and B.q < 0:	# add lj test
			lj = self.lja/r/r/r/r/r/r
			ret += lj

		return ret

	def get_force(self,A,B):
		ret = []				# data contains force
		r = self.get_dist(A,B)
		dx = A.xyz[0] - B.xyz[0]
		dy = A.xyz[1] - B.xyz[1]
		dz = A.xyz[2] - B.xyz[2]
		fx = self.to_ev_const*A.q*B.q*dx/r/r/r	# FORCE ACTING ON "A"
		fy = self.to_ev_const*A.q*B.q*dy/r/r/r
		fz = self.to_ev_const*A.q*B.q*dz/r/r/r
		ret = [fx,fy,fz]

		if A.q*B.q < 0.:	# if pair charge are not in the same sign > add short_range repulsion term (Born-Mayer)
			rep_x = self.bma/self.bmr*math.exp(-r/self.bmr)/r*dx
			rep_y = self.bma/self.bmr*math.exp(-r/self.bmr)/r*dy
			rep_z = self.bma/self.bmr*math.exp(-r/self.bmr)/r*dz
			ret[0] += rep_x
			ret[1] += rep_y
			ret[2] += rep_z

		if A.q < 0 and B.q < 0:	# add lj test

			lj_x = self.lja*6./r/r/r/r/r/r/r/r*dx
			lj_y = self.lja*6./r/r/r/r/r/r/r/r*dy
			lj_z = self.lja*6./r/r/r/r/r/r/r/r*dz
			ret[0] += lj_x
			ret[1] += lj_y
			ret[2] += lj_z

		return ret

	####################################################################################################################################################
	# BEGIN RANDOM CONFIG GEN
	####################################################################################################################################################

	def build_rand_config(self):

		self.atoms = []			# FOR RECORDING
		noa = 0				# Number of atoms
		random.seed(time.time())	# Random seed 

		ret = True			# Initialise retur 'True'

		for i in range(self.size):	# for the size of building block (unit)

			# 1st element; in 'self.stochio[0]'
			for j in range(self.stochio[0]):
				for k in range(self.max_trial):
				
					if k == self.max_trial - 1:	# if max-trial reached ... then return 'False' failed to find an appropriate configuration
						ret = False

					rand_x = random.uniform(0.,self.box_len[0])						#
					rand_y = random.uniform(0.,self.box_len[1])						#
					rand_z = random.uniform(0.,self.box_len[2])						#
					self.atoms.append( atom(self.element[0],self.charge[0],[rand_x,rand_y,rand_z]))		#
					noa += 1										# generate random atom position in the box
					fc = [0.,0.,0.]										# force initialisation

					if noa == 1:	# if only '1' atom in the configuration ... move on to the next atom generation
						break	#

					elif noa > 1:			# if more then 1 atom in the box; then start checking 'interatomic distance' & 'interatomic forces'
						dist_fail = False	# flag if fail in interatomic distance

						for ii in range(noa-1):	# '-1' for excluding self distance
							dist = self.get_dist(self.atoms[ii],self.atoms[noa-1])		# get distance b/t the new atom ([noa-1]) vs all the others

							if dist < self.inter_atom_dist_tol:	# if distance within the tolerance
								noa -= 1			# reject
								self.atoms.pop()		# delete new atom just added
								dist_fail = True		# set fail flag = 'true'
								break				# back to retry
						
							force = self.get_force(self.atoms[noa-1],self.atoms[ii])	# get force acting on the new atom
							fc[0] += force[0]						# by all the other atoms in the box
							fc[1] += force[1]						#
							fc[2] += force[2]						#

						if dist_fail == False:							# if distance check done
							net_force = math.sqrt(fc[0]*fc[0]+fc[1]*fc[1]+fc[2]*fc[2])	# get net force acting on the new atom
							if net_force > self.inter_atom_force_tol:			# if force check failed
								noa -= 1						# reject trial atom
								self.atoms.pop()					# delete new atom just added
							else:				# if force check done	
								break			# escape

			# 2nd element; repeat the same sequence done above
			for j in range(self.stochio[1]):

				for k in range(self.max_trial):

					if k == self.max_trial - 1:
						ret = False

					rand_x = random.uniform(0.,self.box_len[0])
					rand_y = random.uniform(0.,self.box_len[1])
					rand_z = random.uniform(0.,self.box_len[2])
					self.atoms.append( atom(self.element[1],self.charge[1],[rand_x,rand_y,rand_z]))
					noa += 1
					fc = [0.,0.,0.]

					if noa == 1:
						break
					elif noa > 1:
			
						dist_fail = False
						for ii in range(noa-1):
							dist = self.get_dist(self.atoms[ii],self.atoms[noa-1])

							if dist < self.inter_atom_dist_tol:
								noa -= 1			
								self.atoms.pop()		
								dist_fail = True
								break			
						
							force = self.get_force(self.atoms[noa-1],self.atoms[ii])
							fc[0] += force[0]
							fc[1] += force[1]
							fc[2] += force[2]

						if dist_fail == False:

							net_force = math.sqrt(fc[0]*fc[0]+fc[1]*fc[1]+fc[2]*fc[2])
						
							if net_force > self.inter_atom_force_tol:
								noa -= 1				
								self.atoms.pop()			
							else:
								break
		
		self.noa = noa # save number of atoms 

		return ret
	# End of Rand Gen 
	####################################################################################################################################################



	####################################################################################################################################################
	# Some functions for geometric optimisation; energy_update, force_update, gnorm_update, line-search sub-problem solver, backtracking line-search
	####################################################################################################################################################

	def energy_update(self):								# get energy of current config
		self.sys_e = 0.									#
		for i in range(len(self.atoms)):						#
			for j in range(i+1,len(self.atoms)):					#
				pair_energy = self.get_energy(self.atoms[i],self.atoms[j])	#
				self.sys_e += pair_energy					#
		return self.sys_e								# return current energy

	def force_update(self):									# get force of current config
		# force initilisation
		for i in range(len(self.atoms)):
			self.atoms[i].f[0] = 0.
			self.atoms[i].f[1] = 0.
			self.atoms[i].f[2] = 0.
		# calculate force
		for i in range(len(self.atoms)):
			for j in range(i+1,len(self.atoms)):
		
				force = self.get_force(self.atoms[i],self.atoms[j])
			
				self.atoms[i].f[0] += force[0]
				self.atoms[i].f[1] += force[1]
				self.atoms[i].f[2] += force[2]

				self.atoms[j].f[0] -= force[0]
				self.atoms[j].f[1] -= force[1]
				self.atoms[j].f[2] -= force[2]

	# updated force are saved in 'atoms' (sub-class)

	# get gnorm of current config ... will be used to tell if stop geometric optimisation
	def gnorm_update(self):
		self.gnorm = 0.
		for i in range(len(self.atoms)):
			self.gnorm += (self.atoms[i].f[0]*self.atoms[i].f[0] + self.atoms[i].f[1]*self.atoms[i].f[1] +self.atoms[i].f[2]*self.atoms[i].f[2])
		self.gnorm = math.sqrt(self.gnorm)/float((len(self.atoms)*3))
		return self.gnorm
	
	
	def ls_sub(self):				# line-search subproblem solver
		t = 1; beta = 0.84			# some internal parameters

		while(True):				# if initial step size (atoms to move) is too huge, constrain step-size
			too_huge_step = False
			for i in range(len(self.atoms)):
				if math.fabs(t*self.atoms[i].f[0]) > 0.1 or math.fabs(t*self.atoms[i].f[1]) > 0.1 or math.fabs(t*self.atoms[i].f[2]) > 0.1:
					t = t*beta
					too_huge_step = True
					break
			if too_huge_step == False:
				break

		original_energy = self.energy_update()	# get energy of current config
		self.force_update()			# get forces of current config

		tot_norm_sq = 0.			# ls-sub parameter
		for i in range(len(self.atoms)):
			tot_norm_sq += (self.atoms[i].f[0]*self.atoms[i].f[0] + self.atoms[i].f[1]*self.atoms[i].f[1] +self.atoms[i].f[2]*self.atoms[i].f[2])
		
		'''
			SUB-PROBLEM SOLVER MAIN
		'''
		for i in range(self.ls_max_trial):

			for j in range(len(self.atoms)):						# MOVE ATOMS TOWARD WHERE FORCE IS IN ACTION
				self.atoms[j].xyz[0] = self.atoms[j].xyz[0] + t*self.atoms[j].f[0]	#
				self.atoms[j].xyz[1] = self.atoms[j].xyz[1] + t*self.atoms[j].f[1]	#
				self.atoms[j].xyz[2] = self.atoms[j].xyz[2] + t*self.atoms[j].f[2]	#
			trial_energy = self.energy_update()						# get trial energy
			backtrack_std = original_energy - 0.5*t*tot_norm_sq				# get acceptance condition parameter

			if trial_energy < backtrack_std:						# ls-sub solved; update configuration
				return True								# escape

			else:											# if rejected
				for j in range(len(self.atoms)):						# move atoms back to the original config
					self.atoms[j].xyz[0] = self.atoms[j].xyz[0] - t*self.atoms[j].f[0]	#
					self.atoms[j].xyz[1] = self.atoms[j].xyz[1] - t*self.atoms[j].f[1]	#
					self.atoms[j].xyz[2] = self.atoms[j].xyz[2] - t*self.atoms[j].f[2]	#
				t = beta*t									# reduce step-size

		return False	# if max trial reached ... then return False

	'''
		geometric optimisation main
	'''
	def pre_opti(self):

		self.energy_update()	# update energy
		self.force_update()	# update force
		self.gnorm_update()	# update gnorm
	
		if self.gnorm < self.gnorm_tol:	# if gnorm < gnorm_tol
			return True		# escape

		for i in range(self.opti_max):
		
			if i%32 == 0:
				print(" Cycle:\t%6.3d / Gnorm:\t%12.5f" % (i+1,self.gnorm))				# record opti process

			res = self.ls_sub()

			if res == True:				# ls-sub if solved
				self.energy_update()		# update all
				self.force_update()		#
				self.gnorm_update()		#
			
				if self.gnorm < self.gnorm_tol:	# if meets termination condition
					#print(" Cycle:\t%6.3d / Gnorm:\t%12.5f ... Final Config" % (i+1,self.gnorm))			#
					return True		#
				elif self.gnorm > self.gnorm_div_tol: # if gnorm diverges
					return False

			else:					# if failed to solve ls-subproblem
				print(" No lower configuration is found ")
				return False
	
	# END OF GEOMETRIC OPTIMISATION

	'''
		record output
	'''

	def print_random_generator_info(self):

		print("")
		print(" BOX_INFO(x/y/z):\t %12.2f\t%12.2f\t%12.2f" % (self.box_len[0],self.box_len[1],self.box_len[2]))
		print(" ELEMENT_UNIT:\t %4.3s%d%4.3s%d" % (self.element[0],self.stochio[0],self.element[1],self.stochio[1]))
		print(" ELEMENT_Q   :\t %.3f  %.3f"    % (self.charge[0],self.charge[1]))
		print("")
		print(" ALLOWED_INITIAL_INTERATOMIC_DISTANCE(Angs):\t%12.4f" % (self.inter_atom_dist_tol))
		print(" ALLOWED_INITIAL_INTERATOMIC_FORCE(eV/Angs):\t%12.5f" % (self.inter_atom_force_tol))
		print("")
		print(" MAXIMUM_RANDOM_GENERATION_TRIAL:\t%d" % (self.max_trial))
		print(" MAXIMUM_OPTIMISATION_CYCLE:     \t%d" % (self.opti_max))
		print(" MAXIMUM_LS_SUBPROBLEM_TRIALS:   \t%d" % (self.ls_max_trial))
		print("")
		print(" PRE-OPTIMISATION_Born_Mayer_A(Cation-Anion):  \t%.4f" % (self.bma))
		print(" PRE-OPTIMISATION_Born_Mayer_R(Cation-Anion):  \t%.4f" % (self.bmr))
		print(" PRE-OPTIMISATION_Vdw(Anion-Anion):            \t%.4f" % (self.lja))
		print("")

	def print_xyz(self):			

		with open("racg_out.xyz","w") as f:
		
			self.energy_update()
			self.force_update()
			f.write("\t%d\n" % (self.noa))
			f.write("SCF DONE\t %.6f\n" % (self.sys_e))
			for i in range(self.noa):
				f.write("%.3s\t%12.6f\t%12.6f\t%12.6f\n" % (self.atoms[i].name,self.atoms[i].xyz[0],self.atoms[i].xyz[1],self.atoms[i].xyz[2]))
			f.write("\nForce\n")
			self.gnorm_update()
			f.write("gnorm: %.6lf\n" % (self.gnorm))
			for i in range(self.noa):
				f.write("%.3s\t%18.6e\t%18.6e\t%18.6e\t%10.4f\n" % (self.atoms[i].name,self.atoms[i].f[0],self.atoms[i].f[1],self.atoms[i].f[2],self.atoms[i].q))

		
		



if __name__ == '__main__':

	random_config_generator = RACG()
	generation_result = random_config_generator.build_rand_config()
	
	if generation_result == True:
		print("#"*100)
		print(" Random Atom Configuration Generator version 1")
		print(" Author: Woongkyu Jee, woong.jee.16@ucl.ac.uk")
		print(" Last Edited: 17 / 12 / 2020")
		print("")
		print(" Random generation is successfully done ! ")
		print(" Initialising geometric optimiser ...")
		print("")
		print(" General Input Info")
		random_config_generator.print_random_generator_info()
		print("#"*100)
		print("")
		print(" optimisation log\n")
		opti_result = random_config_generator.pre_opti()
		print("")
		if opti_result == True:
			print(" optimisation done, find detailed info in 'out.xyz'\n")
			print("#"*100)
			print(" RESULT_TAG: SUCCESS ")
			print("#"*100)
		else:
			print(" optimisation is finished but structure may not in a gnorm below the tolerance")
			print(" find detaild info in out.xyz\n")
			print("#"*100)
			print(" RESULT_TAG: FAIL ")
			print("#"*100)
		
		random_config_generator.print_xyz()

	else:
		print("#"*100)
		print(" Failed to generate initial random config !")
		print(" Try tweaking parameters below (see input file 'RACG.in')")
		print("#"*100)
		print("")
		print(" 1. 'BOX_SIZE'\n")
		print(" 2. 'INTER_ATOMIC_DISTANCE_TOL(Angs)\n")
		print(" 3. 'INTER_ATOMIC_FORCE_TOL(eV/Angs)\n")
		print("#"*100)
		print(" RESULT_TAG: FAIL ")
		print("#"*100)
