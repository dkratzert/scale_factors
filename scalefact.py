#! python
#-*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
#  "THE BEER-WARE LICENSE" (Revision 42):
#  <dkratzert@gmx.de> wrote this file. As long as you retain this notice you
#  can do whatever you want with this stuff. If we meet some day, and you think
#  this stuff is worth it, you can buy me a beer in return Daniel Kratzert.
#  ----------------------------------------------------------------------------
# This script devides a given XD2006 hkl file in different scale factors batches
#
# Todo: ability to define scale resolution batches
# Translations:
# erstes = first element
#
#la = float(0.71073) #lambda (wellenlaenge)

import sys
import string
import os
import logging
# Python version checks:
# This script doesn't work with the old Python version provided with APEXII
print "Python version is %s.%s.%s " %(sys.version_info[0], sys.version_info[1], sys.version_info[2])
if int(sys.version_info[0]) <= 2:
	if int(sys.version_info[1]) < 7:
		print "Your python version is %s.%s.%s Please Use version 2.7.x or greather, but not >= 3.0.x!"\
		%(sys.version_info[0], sys.version_info[1], sys.version_info[2])
		sys.exit(-1)

import math as np
from argparse import ArgumentParser

#
# The argumentparser for the command line options:
#

parser = ArgumentParser(description='This script devides a given XD2006 hkl file in different scale factors batches')
parser.add_argument("-f", dest="hkl_file", 
					metavar='hkl file', help=".hkl input file")
					
parser.add_argument("-p", dest="p4p_file", 
					metavar='p4p/mas file', 
					default=False, 
					help=".p4p or .mas file with cell data")
					
parser.add_argument("-o", dest="out_file", 
					metavar='hkl file', 
					help=".hkl output file")
					
parser.add_argument("-k", dest="scalef", 
					metavar='"start step numscales"', 
					default=False, 
					help="Define own resolution shells in sin(theta/lambda) with a start value, size of step, \
					and quantity of scale factors with space separation in quotation marks (prohibits -a, -s).")
					
parser.add_argument("-K", dest="defscale", 
					metavar='"resolution shells"', 
					default=False,
					help="Define the resolution shells space separated in quotation marks in any order.")
					
parser.add_argument("-a", dest="resola", action='store_true', 
					default=False, 
					help="Write resolution in Angstroem instead of batch-numbers.")
					
parser.add_argument("-s", dest="resols", action='store_true', 
					default=False, 
					help="Write resolution in sin(theta/lambda) instead of batch-numbers.")
					
options = parser.parse_args()

# options.resola ist True wenn -a angegeben ist

#
# some checks if the command line arguments are correct
# 

if options.hkl_file is None:
	print "\nPlease give the input hkl-file as argument!\n"
	parser.print_help()
	sys.exit(-1)
if options.p4p_file is False:
	print "\nPlease give the input p4p-file as argument!\n"
	parser.print_help()
	sys.exit(-1)
if options.out_file is None:
	print "\nPlease give the output hkl-file as argument!\n"
	parser.print_help()
	sys.exit(-1)
#if options.scalef is not False:1
#	if options.resola is True or options.resols is True:
#		print "\nThe options -a and -s prohibit the option -k!"
#		parser.print_help()
#		sys.exit(-1)
if options.resola is False and options.resols is False and options.defscale is False and options.scalef is False:
	print "\nPlease give more arguments!\n"
	parser.print_help()
	sys.exit(-1)

	
	
def find_line(lines, tofind):
   for i, line in enumerate(lines):
      if line.find(tofind) >= 0:
         return i #gives back the line found in tofind

		 

#
# calculates the resolution of the refelctions:
#		 
def triklin(h, k ,l, a, b, c, alpha, beta, gamma, la):
	'''gives sin(theta)/lambda'''
	v = a * b * c * np.sqrt( 1 + 2 * (np.cos(alpha) * np.cos(beta) * np.cos(gamma)) - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 )	
	astar = (1/v) * b * c * np.sin(alpha)
	bstar = (1/v) * c * a * np.sin(beta)
	cstar = (1/v) * a * b * np.sin(gamma)
	
	alphastar = ( np.cos(beta) * np.cos(gamma) - np.cos(alpha) ) / (np.sin(beta) * np.sin(gamma) )
	betastar  = ( np.cos(gamma) * np.cos(alpha) - np.cos(beta) ) / (np.sin(gamma) * np.sin(alpha) )
	gammastar = ( np.cos(alpha) * np.cos(beta) - np.cos(gamma) ) / (np.sin(alpha) * np.sin(beta) )
				
	teil1 = (h**2 * astar**2 +	 k**2 * bstar**2 + l**2 * cstar**2)
	teil2 = 2.0 * k * l * bstar * cstar * alphastar
	teil3 = 2.0 * l * h * cstar * astar * betastar
	teil4 = 2.0 * h * k * astar * bstar * gammastar
	
	sintsq = ( (la**2 / 4.0) * (teil1 + teil2 + teil3 + teil4) )
	
	sinthl = np.sqrt(sintsq) / la
	
	if options.resola is True:
		sinthl = 1.0 / (2.0 * sinthl)  #if you want resolution in d
	return sinthl
		 
		 

def sinthl_zeilen(hkl):
	'''calculates a list with the resolutions of each reflection'''
	print "calculating resolutions..."
	while True: 
		line = []
		line = hkl.readline()
		if not line: 
			break
		line = line.split()
		if line[2] == 'NDAT':  #leave header out of the list
			continue
		h = int(line[0])
		k = int(line[1])
		l = int(line[2])
		
		sinthl_list.append(triklin(h, k ,l, a, b, c, alpha, beta, gamma, la))
	sinthl_list.reverse()


	
	
def scale_zeilen(hkl, steps):
	'''calculates the resolution batches and fills in the refelctions'''
	print "calculating batches..."
	keylist = []
	
	while True: 
		line = []
		line = hkl.readline()
		if not line: 
			break
		line = line.split()
		if line[2] == 'NDAT':
			line.append('\n')
			outfile.write(' '.join(line))
			continue
		sinthlval = sinthl_list.pop()

		if options.resola is True or options.resols is True:
			# fill all informations of h,k,l,I,sig,resolution in out:
			out = '%s %s %s %s %s  %s\n'%(str(line[0]).rjust(4), str(line[1]).rjust(3), str(line[2]).rjust(3), str(line[4]).rjust(7), str(line[5]).rjust(7), str(round(sinthlval, 4)).ljust(6, '0'))
			outfile.write(out)
			continue
		for key in scales.iterkeys():
			#print key
			erstes = float(scales[key])
			if int(key)+1 not in scales:
				zweites = float(99)
			else:
				zweites = float(scales[int(key)+1])
			
			if erstes < sinthlval < zweites:
				# h k l batch Fsq sig Tbar
				keylist.append(key) #make a list with resolution batches
				
				out = '%s %s %s %s %s %s %s\n'%(str(line[0]).rjust(4), str(line[1]).rjust(3), str(line[2]).rjust(3), str(key+1).rjust(2), str(line[4]).rjust(7), str(line[5]).rjust(7), "1.0000".rjust(4))
				outfile.write(out)
		
	for i in range(steps):
		logging.info('Refl. in batch %2s: %s',i+1 , keylist.count(i) )
			
				
def create_scalelist(defscale):
	'''creates a list of scales from command line input'''
	dscale = options.defscale.split()
	dscale = [float(i) for i in dscale]
	if float(0) not in dscale:
		dscale.append(0.0)
	dscale.sort()
	return dscale
	

def drange(start, step, stop):
	# count from start to stop with float(step)
	r = start
	while r < stop:
		yield r
		r += step

	
		
if __name__ == "__main__":	
	'main function'
	logfile = str(options.out_file+'.lst')
	logging.basicConfig(filename=logfile, filemode='w', level=logging.INFO, format='%(message)s') #format of logging
	console = logging.StreamHandler() #define console logging
	logging.getLogger('').addHandler(console) #adds the console logger
	logging.info('Scale factor calculation script by Daniel Kratzert (dkratzert@gmx.de)\n')
	logging.info('Command line arguments given given:')
	logging.info('hkl-file: %s, p4p/mas-file: %s, output-file: %s\n', options.hkl_file, options.p4p_file, options.out_file)
	if options.resola is True:
		logging.info('writing resolution batches in Angstroem instead of batch numbers.')
	if options.resols is True:
		logging.info('writing resolution batches in sin(theta/lambda) instead of batch numbers.')
		
	if options.defscale is not False:
		# this runs if we define scale factor shells freely:
		scales_list = create_scalelist(options.defscale)
		key = []
		key = range(len(scales_list))
		scales = dict(zip(key, scales_list))
		#print 'scale steps given:', scales.values()
		scalesteps = scales.values()
		
		printsteps = scalesteps
		printsteps.append('inf')
		#print printsteps
		logging.info('scale steps given: %s (%s batches)', printsteps, len(printsteps))
	
	if options.scalef is not False:
		#this runs if we define "start step number" scale factors
		so = options.scalef.split()
		
		scalesteps = []
		start = float(so[0])
		step = float(so[1])
		scalesteps.append(start)
		if float(0) not in scalesteps:
			scalesteps.append(0.0)

		numscales = float(so[2])
		#print scalesteps
		scalesteps.sort()
		# sum up all scale factors until limit reached:
		for i in scalesteps:
			if len(scalesteps) < numscales:
				scalesteps.append(round(scalesteps[-1]+step, 2))
		
		scalesteps.sort()
		printsteps = []
		
		for i in range(len(scalesteps)):
			printsteps.append(scalesteps[i])
		
		printsteps.append('inf')
		
		angststeps = []
		for i in printsteps:
			if i > 0:
				angststeps.append(round(1.0 / (2.0 * float(i)), 2))
		angststeps.append(999)
		angststeps.sort()
		#print angststeps
		#for x in printsteps:
		#	print repr(x).rjust(2)
		#print '%s %s %s%s %s' %('using scale steps (in sintl):', printsteps, " (", len(printsteps)-1, "scale factors)")
		#print '%s %s' %('using scale steps (in angstrom):', angststeps)
		logging.info('using scale steps (in sintl):    %s (%s scale fact.)', printsteps, len(printsteps)-1)
		logging.info('using scale steps (in angstrom): %s', angststeps)
		
		key = []
		key = range(len(scalesteps))
		scales = dict(zip(key, scalesteps))
	

	
	sinthl_list = []
	cell_data = []
	p4pfile = open( options.p4p_file, "r" ).readlines(10)
		
	# find the cell informations:
	cellnum = find_line(p4pfile, 'CELL ')
	cell_data = p4pfile[cellnum].split()[0:7]
	del cell_data[0]
	
	
	
	if int(options.p4p_file.lower().find('.mas')) != -1:
		# we have a XD-Masterfile
		#print "we have a XD-Masterfile"
		wavenum = find_line(p4pfile, 'WAVE')
		la = float(p4pfile[wavenum].split()[1])
		#print "\nWavelength is:", la
		logging.info('\nWavelength is: %s', la)
	elif int(options.p4p_file.lower().find('.p4p')) != -1:
		#print "we have a p4p file"
		# we have a p4p file
		sourcenum = find_line(p4pfile, 'SOURCE')
		la = float(p4pfile[sourcenum].split()[2])
		#print "\nWavelength is:", la
		logging.info('\nWavelength is: %s', la)
	else:
		#print "\nCould not find a valid cell and wavelength in %s\n" %(options.p4p_file)
		logging.info('\nCould not find a valid cell and wavelength in %s\n', options.p4p_file)
		sys.exit(-1)
		
		
	# das sind die Zelldaten:
	a = float(cell_data[0])
	b = float(cell_data[1])
	c = float(cell_data[2])
	# radians are important here, because sin, cos calculate in radians!!!
	alpha = np.radians(float(cell_data[3]))
	beta = np.radians(float(cell_data[4]))
	gamma = np.radians(float(cell_data[5]))
	
	cellparam = '  '.join(cell_data[:])
	#print "Cell constants:", '  '.join(cell_data[:])
	logging.info('\nCell constants: %s\n', cellparam)
	
	hkl = open( options.hkl_file, "r" )  # read in the orginal hkl
	sinthl_zeilen(hkl)
	
	hkl = open( options.hkl_file, "r" )
	outfile = open( options.out_file, "w" )
	scale_zeilen(hkl, len(scalesteps))
	print "finished!"
	