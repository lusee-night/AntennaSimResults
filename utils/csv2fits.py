#!/usr/bin/env python
import fitsio
import numpy as np
import argparse

import argparse

parser = argparse.ArgumentParser(description="""Process antenna files into a fits file. 
            Specify csv file format on input, e.g ElectricFieldPlot_Freq__freq__MHz.csv.""")
parser.add_argument('input', help='input csv file format, __freq__ replaces frequency')
parser.add_argument('-o', dest='output', type=str, help="output file")
parser.add_argument('--freq_start', dest='fstart', default = 1, type=int, help="output file")
parser.add_argument('--freq_end', dest='fend', default = 50, type=int, help="output file")

args = parser.parse_args()
informat = args.input
outfile = args.output


allE = []

for freq in range(args.fstart, args.fend+1):
    fname = informat.replace('__freq__',str(freq))
    print (f"Processing {fname}...", end="")
    dataR = np.recfromcsv(fname, delimiter=',')
    dataI = np.recfromcsv(fname.replace("in_phase","out_of_phase"), delimiter=',')
    thetadeg = dataR['thetadeg']
    phideg = dataR['phideg']
    Ntheta = len(set(thetadeg))
    Nphi = len(set(phideg))
    dtheta = 180/(Ntheta-1)
    dphi = 360/(len(set(phideg))-1)
    thetandx = ((thetadeg-0)/dtheta).astype(int)
    phindx = ((phideg-(-180))/dphi).astype(int)
    E = np.zeros((Ntheta,Nphi,3),np.complex128)+np.nan ## Make them nan to see if they are filled
    Ereal = np.array([dataR['rerexmv'],dataR['rereymv'],dataR['rerezmv']]).T
    Eimag = np.array([dataR['imrexmv'],dataR['imreymv'],dataR['imrezmv']]).T 
    E[thetandx,phindx]=Ereal + 1j*Eimag
    assert(not np.any(np.isnan(E)))
    print (f"    read {Ntheta}x{Nphi} data points.")
    allE.append(E)
    
header = {'from_csv':informat,
          'freq_start':args.fstart,
          'freq_end':args.fend,
          'freq_step':1,
          'phi_start':0,
          'phi_end':360,
          'phi_step':dphi,
          'theta_start':0,
          'theta_end':180,
          'theta_step':2}

out = np.array(allE)
print ('Output array shape',out.shape)
print ('Saving to ',outfile)
fits = fitsio.FITS(outfile,'rw')
fits.write(np.real(out), header=header)
fits.write(np.imag(out), header=header)
print ('Done.')
