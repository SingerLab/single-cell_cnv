#!/usr/bin/env python

import sys
import time
import string
import os
import subprocess
import gzip
from itertools import imap
import operator


def main():

	outdir = sys.argv[1]
	barcodefilename = sys.argv[2]
	infilename1 = sys.argv[3]
	lane = sys.argv[4]

	bcLength = 11	
	
	BARCODE = open(barcodefilename, "r")
	INFILE1 = open(infilename1, "r")
	
	barcode = dict()
	barcode_list = []
	for x in BARCODE:
		arow = x.rstrip().split()
		barcode_list.append(arow)
		thisBarcodeId = arow[0]
		thisBarcode = arow[1]
		if barcode.has_key(thisBarcode):
			print "ERROR: Duplicate barcode."
			return -1
		else:
			barcode[thisBarcode] = thisBarcodeId
	
	#print barcode		


	OUTFILES1 = dict()
	for k, v in barcode.items():
		thisOutfilename = outdir + "/s_" + lane + "_1_sequence." + v + ".txt.gz"
		OUTFILES1[v] = gzip.open(thisOutfilename, 'wb')
	thisOutfilename = outdir + "/s_" + lane + "_1_sequence." + "0" + ".txt.gz"
	OUTFILES1["0"] = gzip.open(thisOutfilename, 'wb')
		
	thisRead1 = None
	thisRead1 = []
	thisRead2 = None
	thisRead2 = []
	
	thisBcId1 = ""
	thisBc1 = ""
	outputBc = ""
	counter = 0
	
	for aline1 in INFILE1:
		aline1 = aline1.rstrip()

		if counter % 4 == 0:

			##  WRITE PREVIOUS READ  ##
			if outputBc != "":
				for j in range(len(thisRead1)):
					OUTFILES1[outputBc].write(thisRead1[j])
					OUTFILES1[outputBc].write("\n")

			##  START NEW READ  ##
			thisBcId = ""
			thisRead1 = None
			thisRead1 = []
			thisRead1.append(aline1)
			
		elif counter % 4 == 1:

			outputBc = "0"
			thisStart = aline1.find("TGTGTTGGG")
			if thisStart > 12:
				thisBc1 = aline1[(thisStart - 11):thisStart]
				if barcode.has_key(thisBc1):
					outputBc = barcode[thisBc1]
			else:
				for i in range(len(barcode_list)):
					if hamming(thisBc1, barcode_list[i][1]) < 3:
						outputBc = barcode_list[i][0]
						break
				
			thisRead1.append(aline1)

		elif counter % 4 == 2:
			thisRead1.append(aline1)

		else:
			thisRead1.append(aline1)

		counter += 1
		
	INFILE1.close()
	BARCODE.close()
	for k in OUTFILES1.keys():
		if OUTFILES1[k]:
			OUTFILES1[k].close()
			

def hamming(str1, str2):
    ne = operator.ne
    return sum(imap(ne, str1, str2))


if __name__=="__main__":
	main()
