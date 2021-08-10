#!/bin/python

from math import fabs
import sys

_e_tol = 0.0001
arg_len = len(sys.argv)

if arg_len == 2:

	if sys.argv[1] == "-help":
		print("AUTO READ RACG_statistics, just run!")
		print("Format of input file has to be in two colum base, 1st: structure-code 2nd: energy")
		print("The input file must be pre-sorted (in descending order): use 'sort -k2 ...")
		print("Put input file through command line argument")	
		sys.exit(1)


try:
	with open(sys.argv[1],"r") as f:
		data = []
		for line in f:
			spl = line.split()
			spl[1] = float(spl[1])
			data.append(spl)
		data.reverse()
except FileNotFoundError:
	print("input file is not found ..")
	sys.exit(1)

##############################################################################################################

freq = []
code = []
ener = []
tmp_freq = 0

cur_code = data[0][0]		# FIRST CNT
cur_ener = data[0][1]
tmp_freq += 1

code.append(cur_code)		# APPEND CODE / ENER
ener.append(cur_ener)

for i in range(len(data)-1):

	if i == len(data) - 2:
		freq.append(tmp_freq)
		break

	next_ener = data[i+1][1]
	diff = fabs( cur_ener - next_ener )

	if diff < _e_tol :	# if within torlerance (i.e., same structure)
		tmp_freq += 1
	else:
		freq.append(tmp_freq)		# update freq

		cur_code = data[i+1][0]		# set new struct
		cur_ener = data[i+1][1]
		code.append(cur_code)
		ener.append(cur_ener)

		tmp_freq = 1


# RES PRING
print("\t",len(freq))

for i in range(len(freq)):
	print("%12.9s\t%20.6f\t%10d" % (code[i],ener[i],freq[i]))
