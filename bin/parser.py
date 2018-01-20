#///////////////////////////////////////////////////////////////////
#//            Parses input for myQC 
#//
#//              James H Thorpe, in Group of John Stanton
#//              The University of Florida
#//
#///////////////////////////////////////////////////////////////////

'''
options array
 0:  0		molecule read type	:0- cartesian (default)
 1:  0		calculation type	:0-SCF(default)
 2:  0		basis set 		:0-STO-3G(default)
 3:  0-2	SCF reference		:0-RHF(default),1-UHF,2-ROHF
 4:  0		parallel algorithn	:0-none(default)
 5:  1-?	number of nodes		:integer,1(default)
 6:  any	memory in MB		:1000(default)
 7:  0-3  	verbosity		:0-none,1-some(default),2-detail,3-wtf
 8:  0-12	scf convergence (10^x)	:0-12, 7 (default)
 9:  any	charge			:any, 0(default)
10:  any	multiplicity		:any, 1(default)
11:  0-1	input units		:0-angstroms(default),1-bohr
'''

#=====================================================================
#                       MAIN 

def main():
  import os
  import os.path

  dir_path = os.getcwd()
  print "WE ARE IN :", dir_path

  #check input exists
  if (not os.path.isfile('input')):
    print "parser could not open input file."
    os.system('touch error')
    return 

  #get number of lines
  flines = file_len('input')

  #defaults and setup
  options = [0,0,0,0,0,1,1000,1,7,0,1,0]
  nline = 0
  atoms = []
  xyz = []
  nnuc = 0

  print "=================================="
  print "Input parameters"

  #read input file 
  with open('input') as f:
    #cartesian coordinates
    line = (f.readline()).split()
    if (str(line[0]) == 'CARTESIAN'):
      nline,nnuc =  cartesian(f,nline,nnuc,atoms,xyz)
    else:
      print "Sorry, that input style is not implimented yet"
      os.system('touch error')
      return

    #proccess the rest of the output
    read_options(f,flines-nline,options)

  #check geometry
  if checkgeom(xyz,options[11],nnuc):
    build(atoms,xyz,options,nnuc) 
  else:
    os.system('touch error')

  print "==================================" 
  #write to various output files
  print ""
  if options[7] >= 1:
    print "options list"
    for i in range(0,len(options)):
      print(options[i])

#=====================================================================
#                     FUNCTIONS 

#---------------------------------------------------------------------
#		Get number of lines in file 
#---------------------------------------------------------------------
# THANKS, stackoverflow!
def file_len(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i+1

#---------------------------------------------------------------------
#		Parse geometry of a cartesian molecule
#---------------------------------------------------------------------
def cartesian(f,nline,nnuc,atoms,xyz):

  #assign initial values
  flag1 = 0

  #read in the coordinates
  while not flag1:
    nline += 1 
    st = f.readline()
    line = st.split()
    if line[0] == 'END':
      flag1 = 1
    else:
      print line[0]+" "+str(line[1])+" "+str(line[2])+" "+str(line[3])
      nnuc += 1
      atoms.append(getelem(line[0]))
      xyz.append([float(line[1]),float(line[2]),float(line[3])])

  print ""
  f.readline() 		#blank line at the end
  nline += 1
  return nline,nnuc

#---------------------------------------------------------------------
#		Retrun atomic number	
#---------------------------------------------------------------------
def getelem(st):
  import os
  elem = ['H','He','Li','Be','B','C','N','O','F','Ne'] 
  val = -1
  try:
    val = elem.index(st)+1
  except ValueError:
    os.system('touch error')
    print "Sorry, element "+st+" is not implemented yet."
    raise SystemExit(1) 
  return val 

#---------------------------------------------------------------------
#		Get input options	
#---------------------------------------------------------------------
def read_options(f,rline,options):
  import os
  for i in range(0,rline):
    line = f.readline().split()
    #go to different options as needed  
    if not line:
      continue
    else:
      st = str(line[0])
      if st == 'CALC=':
        getcalc(line[1],options)
      elif st == 'BASIS=':
        getbasis(line[1],options)
      elif st == 'CHARGE=':
        getcharge(line[1],options)
      elif st == 'MULTI=':
        getmulti(line[1],options)
      elif st == 'REF=':
        getref(line[1],options)
      elif st == 'PAR=':
        getpar(line[1],options)
      elif st == 'NODES=':
        getnodes(line[1],options)
      elif st == 'MEMORY=':
        getmem(line[1],options)
      elif st == 'VERB=':
        getverb(line[1],options)
      elif st == 'SCF_Conv=':
        getSCF_Conv(line[1],options)
      elif st == 'UNITS=':
        getUnits(line[1],options)
      else:
        print "parser could not understand this input"
    
  return 
#---------------------------------------------------------------------
#		Get calculation type 
#---------------------------------------------------------------------
def getcalc(st,options):
  import os
  calcs = ['SCF']
  val = -1
  try:
    val = calcs.index(st)
    options[1] = val
    print "Calculation : "+calcs[val]
  except ValueError:
    os.system('touch error')
    print "Sorry, calculation "+st+" is not implemented yet."
    raise SystemExit(1) 
  return 
#---------------------------------------------------------------------
#		Get basis set 
#---------------------------------------------------------------------
def getbasis(st,options):
  import os
  basis = ['STO-3G']
  val = -1
  try:
    val = basis.index(st)
    options[1] = val
    print "Basis : "+basis[val]
  except ValueError:
    os.system('touch error')
    print "Sorry, basis "+st+" is not implemented yet."
    raise SystemExit(1) 
  return 
#---------------------------------------------------------------------
#		Get charge 
#---------------------------------------------------------------------
def getcharge(st,options):
  import os
  charge = 0
  try:
    charge = int(st) 
    options[9] = charge
    print "Charge : "+str(charge)
  except ValueError:
    os.system('touch error')
    print "Your input "+st+" is not a valid charge"
    raise SystemExit(1) 
  return 
#---------------------------------------------------------------------
#		Get multiplicity 
#---------------------------------------------------------------------
def getmulti(st,options):
  import os
  multi = 0
  try:
    multi = int(st) 
  except ValueError:
    os.system('touch error')
    print "Your input "+st+" is not a valid multiplicity"
    raise SystemExit(1) 
  if multi > 0:
    options[10] = multi
    print "Multiplicity : "+str(multi)
  else:
    print "Your input "+st+" is not a valid multiplicity"
    raise SystemExit(1) 
  return 
#---------------------------------------------------------------------
#		Get ref 
#---------------------------------------------------------------------
def getref(st,options):
  import os
  ref = ['RHF','UHF','ROHF']
  val = 0
  try:
    val = ref.index(st) 
    options[3] = val
    print "Ref : "+str(ref[val])
  except ValueError:
    os.system('touch error')
    print "Sorry, your reference "+st+" has not been implimented yet"
    raise SystemExit(1) 
  return 
#---------------------------------------------------------------------
#		Get parallel algorithm 
#---------------------------------------------------------------------
def getpar(st,options):
  import os
  par = ['none','OMP','MPI']
  print "No parallel algorithms implimented yet, sorry!"
  return 
#---------------------------------------------------------------------
#		Get number of nodes 
#---------------------------------------------------------------------
def getnodes(st,options):
  import os
  nodes = 1
  print "No parallel algorithms implimented yet, sorry!"
  return 
#---------------------------------------------------------------------
#		Get memory 
#---------------------------------------------------------------------
def getmem(st,options):
  import os
  mem = 0
  try:
    mem = int(st) 
  except ValueError:
    os.system('touch error')
    print "Your input "+st+" is not a valid memory"
    raise SystemExit(1) 
  if mem > 0:
    options[6] = mem
    print "Memory per CPU (MB) : "+str(mem)
  else:
    print "Your input "+st+" is not a valid memory"
  return 
#---------------------------------------------------------------------
#		Get verbosity 
#---------------------------------------------------------------------
def getverb(st,options):
  verb = ['none','some','detailed','wtf']
  val = 0
  try:
    val = int(st)
  except ValueError:
    print "Verbosity : 0 ("+verb[0]+")"
  if (val > 0) :
    print "Verbosity : "+str(val)+" ("+verb[val]+")"
    options[7] = val
  else:
    print "Verbosity : 0 ("+verb[0]+")"
  return 
#---------------------------------------------------------------------
#		Get SCF convergence 
#---------------------------------------------------------------------
def getSCF_Conv(st,options):
  c = 7
  try:
    c = int(st)
  except ValueError:
    print "SCF_Convergence : 7"
  if (c > 0 and c < 13 and not c == 11):
    print "SCF_Convergence : "+str(c) 
    options[8] = c
  elif (c == 11):
    print "SCF_Convergence (these go to) : "+str(c)
    options[8] = c
  else:
    print "SCF_Convergence : 7" 
  return

#---------------------------------------------------------------------
#		Get Units 
#---------------------------------------------------------------------
def getUnits(st,options):
  units = ['Ang','Bohr']
  val = 0
  try:
    val = units.index(st)
    print "Units : "+units[val]
    options[11] = val
  except ValueError:
    print "Units : "+units[0]
  return
#---------------------------------------------------------------------
#		Check radii of molecule
#---------------------------------------------------------------------
def checkgeom(xyz,units,nnuc):
  import math
  flag = 1
  for i in range(0,nnuc):
    for j in range(i+1,nnuc):
      r = (xyz[i][0]-xyz[j][0])**2.0E0
      r = r + (xyz[i][1]-xyz[j][1])**2.0E0
      r = r + (xyz[i][2]-xyz[j][2])**2.0E0
      if (units == 0 and math.sqrt(r) <= 0.2E0):
        print "These atoms are too close (r = "+str(r)+" A) :", str(i)+", "+str(j) 
        flag = 0
      elif (units == 1 and math.sqrt(r) <= 0.377945198):
        print "These atoms are too close (r = "+str(r)+" bohrs ) :", str(i)+", "+str(j)
        flag = 0
  return flag 
#---------------------------------------------------------------------
#		Check build the molecule 
#---------------------------------------------------------------------
def build(atoms,xyz,options,nnuc):
  import os

  COM = [0.00E0, 0.00E0, 0.00E0]
  mass = [1.000, 4.000, 7.000, 9.000, 11.000, 12.000, 14.000, 16.000, 19.000, 20.000]
  A2B = 1.8897161646320724E0

  units = options[11]
  charge = options[9]
  multi = options[10]
  
  #give new values if input was angstroms
  if units == 0:
    for i in range(0,nnuc):
      xyz[i][0] *= A2B 
      xyz[i][1] *= A2B 
      xyz[i][2] *= A2B 

  #get center of mass
  temp = 0.0E0
  for i in range(0,nnuc):
    COM[0] += mass[atoms[i]-1]*xyz[i][0]     
    COM[1] += mass[atoms[i]-1]*xyz[i][1]     
    COM[2] += mass[atoms[i]-1]*xyz[i][2]     
    temp += mass[atoms[i]-1]

  COM[0] = COM[0] / temp
  COM[1] = COM[1] / temp
  COM[2] = COM[2] / temp
  
  #recenter molecules
  for i in range(0,nnuc):
    xyz[i][0] -= COM[0]
    xyz[i][1] -= COM[1]
    xyz[i][2] -= COM[2]

  #write to nucpos file
  with open('nucpos','w') as n:
    for i in range(0,nnuc):
      n.write(str(atoms[i])+"    "),
      n.write(("%.16E" % float(xyz[i][0]))+"    "),
      n.write(("%.16E" % float(xyz[i][1]))+"    "),
      n.write(("%.16E" % float(xyz[i][2]))+"\n")

  #check multiplicity and charge
  unpr = multi-1
  nelc = sum(atoms)-charge
  nelcA = (nelc-unpr)/2
  nelcB = (nelc-unpr)/2
  nelcA += unpr

  if nelcA+nelcB != nelc:
    print "That charge and multiplicity is not allowed."
    os.system('touch error')
  if (nelc < 0 or nelcA < 0 or nelcB < 0):
    print "You have less than 0 electrons :)"
    os.system('touch error')
  if nelcA != nelcB and options[3] == 0:
    print "You cannot use RHF for open shell molecules!"
    os.system('touch error')

  #write to envdat
  with open('envdat','w') as env:
    env.write(str(nnuc)+"  \n")
    env.write(str(nelcA)+"    "+str(nelcB)+"  \n")
    env.write(str(len(options))+"  \n")
    for i in range(0,len(options)):
      env.write(str(options[i])+"  ")
    env.write("  \n")
    env.write("  \n")
    env.write("#number of nuclei  \n")
    env.write("#number of electrons  \n")
    env.write("#length of options array  \n")
    env.write("options array  \n")

  with open('fmem','w') as fmem:
    fmem.write(("%.16E" % float(options[6]))+"  \n")
  
#---------------------------------------------------------------------
#CAll MAIN
main()

