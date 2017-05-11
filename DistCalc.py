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
        
    def PolEnvSelection(input1,input2,PolRad):

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

#
# Reads a gaussian.com
#
#       with open(input1, 'r') as data:
#          for line in data:
#            splt = np.array(line.split()) 
#            if len(splt) == 2 and float(splt[0]).is_integer() and float(splt[0]).is_integer():
#                 for line in data:
#                    ksplt = np.array(line.split())
#                    if len(ksplt) == 4:    
#                       QM_Atom = np.append(QM_Atom,map(float,ksplt[1:4]))      
#
# Reads an amoeba.qm
#
       with open(input1, 'r') as data:
                 for line in data:
                    ksplt = np.array(line.split())
                    QM_Atom = np.append(QM_Atom,map(float,ksplt[1:4]))      
 
       QM_Mat = np.reshape(QM_Atom, (len(QM_Atom)/3,3))
     
 
       f = open('myfile', 'w+')
       header = open(input2,'r')
       line = header.readline()
       f.write(line) 
       header.close()


       with open(input2, 'r+') as data: 
        next(data)
        if args.polfix is not None:
          fixed_set = int(args.polfix)
          for k,line in  enumerate(data):
               if k+1 < fixed_set:
                  f.write(line)
               else:
                  f.write(line)
                  for l,line in  enumerate(data):
#                    if l+1 in MMP_Atom:
#                      f.write(line)
#                    else: 
                     if (l % 3) == 0:
                        splt = np.array(map(float,line.split()))
                        AMOEBA_Atom = np.array(map(int,splt[6:14])) - fixed_set
                        Dist_O1toMMP = np.subtract(splt[2:5], QM_Mat)
                        A = sorted(np.apply_along_axis(np.linalg.norm, 1, Dist_O1toMMP))
                        b = line[232:238]
                        if A[0] <= PolRad:
                            f.write(line)
                            WasOPol = True
                        else:
                            line = string.replace(line, b, '0.0000')
                            f.write(line)
                            WasOPol = False
                     else:
                         if WasOPol is True:
                            f.write(line)
                         else:     
                            b = line[232:238]
                            line = string.replace(line, b, '0.0000')
                            f.write(line)

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


#       print time.clock() - start_time, "seconds"
    if args.polsel is not None:
       PolEnvSelection(file_name_in1,file_name_in2,Rpol)
    


#    if args.actsel is not None:




############################################################
#                                                          #
# Reading the AMOEBA .gau input to get the MM connectivity #
#                                                          #
############################################################

#       with open(input2, 'r') as data:
#           next(data)
#           if args.polfix is not None:
#              fixed_set = int(args.polfix)
#              for k,line in  enumerate(data):
#                 if k+1 < fixed_set:
#                    continue
#                 elif k+1 == fixed_set:
#                    for line in data:
#                      splt = np.array(line.split())
#                      AMOEBA_Atom = np.append(AMOEBA_Atom,map(int,splt[6:14]))
#           else:
#              for line in data:
#                 splt = np.array(line.split())
#                 AMOEBA_Atom = np.append(AMOEBA_Atom,map(int,splt[6:14])) 
#
#       AMOEBA_Mat = np.reshape(AMOEBA_Atom, (len(AMOEBA_Atom)/8, 8))   
