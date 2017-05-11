#!/usr/bin/env python

#################################################################################
#                                                                               #
#                                                                               #
#                                                                               #
#                                                                               #
#                                                                               #
#                                                                               #
#################################################################################

import sys, os
import argparse as arg
import numpy as np
import math as mt
import string

from numpy import linalg as LA

#import time



def option():

    parser = arg.ArgumentParser(description='Perform selection for both active MM atoms and polarizable environment.')
    parser.add_argument('input1', metavar='G', help='Gaussian .com input',nargs='?', default='gaussian.com')
    parser.add_argument('input2', metavar='A', help='.gau input containing the polarizable environment',nargs='?', default='amoeba.gau')
    parser.add_argument('--polsel',metavar='NPol',  help='performe the selection of the polarizable environment around the QM subsystem',nargs='?', default='yes')
    parser.add_argument('--polrad',metavar='RPol',  help='define the radius within the environment around QM is considered polarizable',nargs='?', default='15')
    parser.add_argument('--polrad2',metavar='MMPol',  help='define the radius within the environment around QM is considered iat least MM',nargs='?', default='20')
    parser.add_argument('--actsel', help='performe the selection of the MM atom to define as active for the MM-gradient calculation')
    parser.add_argument('--polfix', help='number of atoms that always have to be polarizable. Be careful when preparing your .gau file!!') 
    args = parser.parse_args()

    return args


if __name__=='__main__':
  
#    start_time = time.clock()
    args = option()

    file_name_in1    =  args.input1
    file_name_in2    =  args.input2
    Rpol             =  float(args.polrad)
    Rpol2            =  float(args.polrad2)    

    def PolEnvSelection(input1,input2,PolRad,PolRad2):

       QM_Atom = np.array([])
       AMOEBA_Atom = np.array([])

######################################################################## 
#                                                                      #
# - MMP_Atom : array of indexes for atoms connected to a polarizable   #
#   oxygen atom, because inside of the polarizable radius              #
#                                                                      #
# - MM_Atom  : array of indexes for hydrogens bonded to a non          #   
#   polarizable oxygen atom, because found to be outside the polariza- #
#   tion radius                                                        #                
#                                                                      #
########################################################################


       MM_Atom = []
       MMP_Atom = []


#######################################################
#                                                     #
# Reading the Gaussian .com to get the QM coordinates #
#                                                     #
#######################################################


       with open(input1, 'r') as data:
          for line in data:
            splt = np.array(line.split()) 
            if len(splt) == 2 and float(splt[0]).is_integer() and float(splt[0]).is_integer():
                 for line in data:
                    ksplt = np.array(line.split())
                    if len(ksplt) == 4:    
                       QM_Atom = np.append(QM_Atom,map(float,ksplt[1:4]))      
       QM_Mat = np.reshape(QM_Atom, (len(QM_Atom)/3,3))
     
 
       f = open('amoeba.gau.mod', 'w+')
#       header = open(input2,'r')
#       line = header.readline()
#       f.write(line) 
#       header.close()

       str_form1 = '%72s%7d%7d%164s%7d%7s%7d%41s'
       add_i = 0
       with open(input2, 'r+') as data: 
        next(data)
        f.write('%-76s\n' % (' '))
        if args.polfix is not None:
          fixed_set = int(args.polfix)
          for k,line in  enumerate(data):
               if k+1 < fixed_set:
                  f.write(line)
               else:
                  f.write(line)
                  for l,line in  enumerate(data):
                     if (l % 3) == 0:
                        splt = np.array(map(float,line.split()))
                        AMOEBA_Atom = np.array(map(int,splt[6:14])) - fixed_set
                        Dist_O1toMMP = np.subtract(splt[2:5], QM_Mat)
                        A = sorted(np.apply_along_axis(np.linalg.norm, 1, Dist_O1toMMP))
                        b = line[232:238]
                        if A[0] <= PolRad:
                            f.write(str_form1 %(line[0:72],int(line[72:79])- add_i,int(line[79:86]) - add_i,line[86:250],int(line[250:257])- add_i,line[257:264],int(line[264:271])- add_i,line[271:312])) 
#                            print add_i
                            WasOPol = True
                        elif (A[0] <= PolRad2):
                            line = string.replace(line, b, '0.0000')
                            f.write(str_form1 %(line[0:72],int(line[72:79])- add_i,int(line[79:86]) - add_i,line[86:250],int(line[250:257])- add_i,line[257:264],int(line[264:271])- add_i,line[271:312])) 
                            WasOPol = False
                            WasOMM  = True 
                        else:
                               WasOPol = False
                               WasOMM  = False 
#                               print 'oxygen',adding 
                               add_i += 1

                     else:
                         if WasOPol is True:
#                              print line 
#                              print add_i
                              f.write(str_form1 %(line[0:72],int(line[72:79])- add_i,int(line[79:86]),line[86:250],int(line[250:257])- add_i,line[257:264],int(line[264:271])- add_i,line[271:312])) 

                         elif (WasOPol is False) and (WasOMM is True):     
                            b = line[232:238]
                            line = string.replace(line, b, '0.0000')
                            f.write(str_form1 %(line[0:72],int(line[72:79])- add_i,int(line[79:86]),line[86:250],int(line[250:257])- add_i,line[257:264],int(line[264:271])- add_i,line[271:312])) 

                         else:
                            add_i += 1
#                            print '@@4', add_i
                            continue 

        else:
           for l,line in  enumerate(data):
             if l+1 in MMP_Atom:
               f.write(line)
             else:
               splt = np.array(map(float,line.split()))
               AMOEBA_Atom = np.array(map(int,splt[6:14])) 
               Dist_O1toMMP = np.subtract(splt[2:5], QM_Mat)
               A = sorted(np.apply_along_axis(np.linalg.norm, 1, Dist_O1toMMP))
               b = line[232:238]
               if l+1 in MM_Atom:
#                  print 'Already in MM_Atom!!'
                  line = string.replace(line, b, '0.0000')
                  f.write(line)
               elif A[0] <= PolRad:
                  f.write(line)
                  if (l+2 in AMOEBA_Atom) and (l+3 in AMOEBA_Atom): 
                      MMP_Atom.append(l+2)
                      MMP_Atom.append(l+3)
               else:
#
#                  print 'Here A[0] too big'
#
                  line = string.replace(line, b, '0.0000')
                  f.write(line)

                  if (l+2 in AMOEBA_Atom) and (l+3 in AMOEBA_Atom):
                    MM_Atom.append(l+2)
                    MM_Atom.append(l+3)
#                  print 'polarization # %d put to zero' %(l+1)


       str_form  = '%8d%68s'
       f.write('  0  11   9999.999999979400   9999.999999979500   9999.999999979600  405      0      0      0      0      0      0      0      0   1.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000    0.0000    0.3900 5      0      0      0    0    0    0    0    0    0    0    0')
       f.seek(0,0)
       header = open(input2,'r')
       line = header.readline()
#       print int(line[0:9])
#       print line
       f.write(str_form %(int(line[0:9])-add_i +1,line[9:76])) 
       header.close()
       f.close()



#       print time.clock() - start_time, "seconds"
    if args.polsel is not None:
       PolEnvSelection(file_name_in1,file_name_in2,Rpol,Rpol2)
    

