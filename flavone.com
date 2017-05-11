%nproc=16
%mem=7000MW
#p b3lyp 6-21G nosymm mmpol=(conver=8,nofmm,inversion) TD=(Nstates=3)  

ONIOM QM/MMPol HL=Cytosine LL=Guanine FREE (Init Geom MP2)
CHARGES = ff02r1

0 1
C                     33.36     29.261    32.214
C                     34.056    29.345    31.071
C                     33.408    29.326    29.848
C                     32.001    29.348    29.812
C                     31.282    29.382    30.966
C                     31.961    29.312    32.175
H                     33.892    29.274    33.161
H                     35.141    29.366    31.084
C                     34.161    29.356    28.543
H                     30.201    29.482    30.927
H                     31.403    29.442    33.098
C                     31.968    29.288    27.422
C                     33.348    29.322    27.316
C                     31.03     29.28     26.257
C                     31.428    29.309    24.901
C                     29.638    29.286    26.531
C                     30.52     29.34     23.884
H                     32.488    29.298    24.668
C                     28.697    29.233    25.545
H                     29.266    29.199    27.547
C                     29.112    29.336    24.212
H                     30.825    29.406    22.844
H                     27.659    29.128    25.847
H                     28.427    29.193    23.382
O                     31.315    29.338    28.598
O                     35.373    29.397    28.461
O                     34.128    29.41     26.227
H                     35.016    29.198    26.564

