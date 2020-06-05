#!/usr/bin/env python

import sys
import time
from operator import itemgetter
import os
import subprocess

def main():

	python_dir = "/YOUR/FILE/PATH/HERE/Python-2.7.1"
	flowcell_dir = "/YOUR/FILE/PATH/DATA/"
	bowtie_dir = "/YOUR/FILE/PATH/HERE/bowtie-0.12.7"
	samtools_dir = "/YOUR/FILE/PATH/HERE/samtools-1.5"
	project_dir = "/YOUR/FILE/PATH/HERE"
	barcoded_dir = “`BarcodeSplit"
	illumina_dir = "/YOUR/FILE/PATH/DATA"
	infilename = project_dir + "/meta/Global.tb.guide.txt"
	mapped_dir = project_dir + "/mapped"
	processed_dir = project_dir + "/processed"
	badbinsfilename = "/YOUR/FILE/PATH/HERE/sequences/hg19.50k.k50.bad.bins.txt"
	gcfilename = "/YOUR/FILE/PATH/HERE/sequences/varbin.gc.content.50k.bowtie.k50.hg19.txt"
	gcfilename5k = "/YOUR/FILE/PATH/HERE/sequences/varbin.gc.content.5k.bowtie.k50.hg19.txt"
	gcfilename20k = "/YOUR/FILE/PATH/HERE/sequences/varbin.gc.content.20k.bowtie.k50.hg19.txt"
	rscript_template = "/YOUR/FILE/PATH/HERE/programs/cbs.varbin.template.asya07.50k.ploidy.txt"
	rscript_template5k = "/YOUR/FILE/PATH/HERE/programs/cbs.varbin.template.asya07.5k.ploidy.txt"
	rscript_template20k = "/YOUR/FILE/PATH/HERE/programs/cbs.varbin.template.asya07.20k.ploidy.txt"
	paa_template = "/YOUR/FILE/PATH/HERE/programs/paa.bash.template03.txt"
	varbin_suffix = ".varbin.50k.txt"
	varbin_suffix5k = ".varbin.5k.txt"
	varbin_suffix20k = ".varbin.20k.txt"
	alpha = 0.02
	nperm = 1000
	undo_SD = 1.0
	min_width = 5

	trim5 = 0
	trim3 = 0
	
	guide = fileToGuide(infilename)
	ncells = len(guide["seq.unit.id"])
	for i in range(ncells):
		#if i > 158:
		#	break
		
		if guide["process"][i] == "0":
			continue
		if guide["trueseq"][i] == "":
			if guide["barcode"][i] == "":
				thisSeqfile = flowcell_dir + "/" + guide["flowcell"][i] + "/s_" + guide["lane"][i] + "_1_sequence.txt.gz"
			else:
				thisSeqfile = barcoded_dir + "/" + guide["flowcell"][i] + "/s_" + guide["lane"][i] + "_sequence." + guide["barcode"][i] + ".txt.gz"
			if os.path.isfile(thisSeqfile):
				pass
			else:
				thisSeqfile = barcoded_dir + "/" + guide["flowcell"][i] + "/s_" + guide["lane"][i] + "_sequence.WB" + guide["barcode"][i] + ".txt.gz"
		else:
			thisSeqfile = illumina_dir + "/" + guide["flowcell"][i] + "/" + guide["trueseq"][i]


		thisSamfile = project_dir + "/mapped/" + guide["seq.unit.id"][i] + ".sam"
		thisRmdupSamfile = project_dir + "/mapped/" + guide["seq.unit.id"][i] + ".rmdup.sam"
		thisBamfile = project_dir + "/mapped/" + guide["seq.unit.id"][i] + ".bam"
		thisSortedBamfile = project_dir + "/mapped/" + guide["seq.unit.id"][i] + ".sorted.bam"
		thisRmdupBamfile = project_dir + "/mapped/" + guide["seq.unit.id"][i] + ".rmdup.bam"
		
		if guide["trueseq"][i] == "":
			if os.path.isfile(thisSeqfile):
				memoryNeeded = str( (os.path.getsize(thisSeqfile) / 400000000 ) + 5) + "G"
			else:
				memoryNeeded = "4G"

			thisCommand = "gunzip -c " + thisSeqfile + " | python solexa.quals01.py"
			#print thisCommand
			p = subprocess.Popen(thisCommand, shell=True, stdout=subprocess.PIPE)		
			thisSolexaQuals = p.communicate()[0].rstrip()
		else:
			thisSolexaQuals = ""
			memoryNeeded = "4G"
			
		
		#print "thisSolexaQuals ", thisSolexaQuals
		
		###  read cbs template

		CBS = open(rscript_template, "r")
		cbs = CBS.read()
		CBS.close()


		varbinFilename = guide["seq.unit.id"][i] + varbin_suffix
		varbinFilename5k = guide["seq.unit.id"][i] + varbin_suffix5k
		varbinFilename20k = guide["seq.unit.id"][i] + varbin_suffix20k

		###  write r script file

		rscript_filename = processed_dir + "/" + guide["seq.unit.id"][i] + ".cbs.r"
		RSCRIPT = open(rscript_filename, "w")
		a = cbs.replace("<INDIR>", "\"" + processed_dir + "\"") 
		a = a.replace("<OUTDIR>", "\"" + processed_dir + "\"") 
		a = a.replace("<BAD.BINS>", "\"" + badbinsfilename + "\"")
		a = a.replace("<VARBIN.GC>", "\"" + gcfilename + "\"") 
		a = a.replace("<VARBIN.DATA>", "\"" + varbinFilename + "\"")
		a = a.replace("<SAMPLE.NAME>", "\"" + guide["seq.unit.id"][i] + "\"")
		a = a.replace("<ALT.SAMPLE.NAME>", "\"" + guide["flowcell"][i] + " lane " + guide["lane"][i] + " " + guide["seq.unit.id"][i] + " " + guide["sample"][i] + " bc" + guide["barcode"][i] + " " + guide["dna.source"][i] + " well " + guide["well"][i] + "\"")
		if guide["mult.min"][i] == "":
			a = a.replace("<MULT_MIN>", "1.5")
		else:
			a = a.replace("<MULT_MIN>", guide["mult.min"][i])
		if guide["mult.max"][i] == "":
			a = a.replace("<MULT_MAX>", "5.5")
		else:
			a = a.replace("<MULT_MAX>", guide["mult.max"][i])
		a = a.replace("<ALPHA>", str(alpha))
		a = a.replace("<NPERM>", str(nperm))
		a = a.replace("<UNDO.SD>", str(undo_SD))
		a = a.replace("<MIN.WIDTH>", str(min_width))

		RSCRIPT.write(a)
		RSCRIPT.close()


		CBS = open(rscript_template5k, "r")
		cbs = CBS.read()
		CBS.close()

		rscript_filename5k = processed_dir + "/" + guide["seq.unit.id"][i] + ".cbs.5k.r"
		RSCRIPT = open(rscript_filename5k, "w")
		a = cbs.replace("<INDIR>", "\"" + processed_dir + "\"") 
		a = a.replace("<OUTDIR>", "\"" + processed_dir + "\"") 
		a = a.replace("<VARBIN.GC>", "\"" + gcfilename5k + "\"") 
		a = a.replace("<VARBIN.DATA>", "\"" + varbinFilename5k + "\"")
		a = a.replace("<SAMPLE.NAME>", "\"" + guide["seq.unit.id"][i] + "\"")
		a = a.replace("<ALT.SAMPLE.NAME>", "\"" + guide["flowcell"][i] + " lane " + guide["lane"][i] + " " + guide["seq.unit.id"][i] + " " + guide["sample"][i] + " bc" + guide["barcode"][i] + " " + guide["dna.source"][i] + " well " + guide["well"][i] + "\"")
		if guide["mult.min"][i] == "":
			a = a.replace("<MULT_MIN>", "1.5")
		else:
			a = a.replace("<MULT_MIN>", guide["mult.min"][i])
		if guide["mult.max"][i] == "":
			a = a.replace("<MULT_MAX>", "5.5")
		else:
			a = a.replace("<MULT_MAX>", guide["mult.max"][i])
		a = a.replace("<ALPHA>", str(alpha))
		a = a.replace("<NPERM>", str(nperm))
		a = a.replace("<UNDO.SD>", str(undo_SD))
		a = a.replace("<MIN.WIDTH>", str(min_width))

		RSCRIPT.write(a)
		RSCRIPT.close()

		CBS = open(rscript_template20k, "r")
		cbs = CBS.read()
		CBS.close()

		rscript_filename20k = processed_dir + "/" + guide["seq.unit.id"][i] + ".cbs.20k.r"
		RSCRIPT = open(rscript_filename20k, "w")
		a = cbs.replace("<INDIR>", "\"" + processed_dir + "\"") 
		a = a.replace("<OUTDIR>", "\"" + processed_dir + "\"") 
		a = a.replace("<VARBIN.GC>", "\"" + gcfilename20k + "\"") 
		a = a.replace("<VARBIN.DATA>", "\"" + varbinFilename20k + "\"")
		a = a.replace("<SAMPLE.NAME>", "\"" + guide["seq.unit.id"][i] + "\"")
		a = a.replace("<ALT.SAMPLE.NAME>", "\"" + guide["flowcell"][i] + " lane " + guide["lane"][i] + " " + guide["seq.unit.id"][i] + " " + guide["sample"][i] + " bc" + guide["barcode"][i] + " " + guide["dna.source"][i] + " well " + guide["well"][i] + "\"")
		if guide["mult.min"][i] == "":
			a = a.replace("<MULT_MIN>", "1.5")
		else:
			a = a.replace("<MULT_MIN>", guide["mult.min"][i])
		if guide["mult.max"][i] == "":
			a = a.replace("<MULT_MAX>", "5.5")
		else:
			a = a.replace("<MULT_MAX>", guide["mult.max"][i])
		a = a.replace("<ALPHA>", str(alpha))
		a = a.replace("<NPERM>", str(nperm))
		a = a.replace("<UNDO.SD>", str(undo_SD))
		a = a.replace("<MIN.WIDTH>", str(min_width))

		RSCRIPT.write(a)
		RSCRIPT.close()

		PAA = open(paa_template, "r")
		paa = PAA.read()
		PAA.close()
		
		mapbash_filename = mapped_dir + "/" + guide["seq.unit.id"][i] + ".paa.bash"
		MAPBASH = open(mapbash_filename, "w")
		a = paa.replace("<BOWTIE_DIR>", bowtie_dir)
		a = a.replace("<PYTHON_DIR>", python_dir)
		a = a.replace("<SAMTOOLS_DIR>", samtools_dir)
		a = a.replace("<SEQUENCE_FILE>", thisSeqfile)
		a = a.replace("<OUTPUT_PREFIX>", guide["seq.unit.id"][i])
		a = a.replace("<TRIM3>", str(trim3))
		a = a.replace("<TRIM5>", str(trim5))
		a = a.replace("<SOLEXA_QUALS>", thisSolexaQuals)
		a = a.replace("<MAPPED_DIR>", mapped_dir)
		a = a.replace("<PROCESSED_DIR>", processed_dir)
		a = a.replace("<PROGRAMS_DIR>", project_dir + "/programs")

		MAPBASH.write(a)
		MAPBASH.close()


		bowtieReportFilename = mapped_dir + "/" + guide["seq.unit.id"][i] + ".bowtie.report.txt"
                qsubFile = mapped_dir + "/" + guide["seq.unit.id"][i] + ".bowtie.qsub"
                qsubFileStdout = mapped_dir + "/" + guide["seq.unit.id"][i] + ".bowtie.qsub.stdout"
                qsubFileStderr = mapped_dir + "/" + guide["seq.unit.id"][i] + ".bowtie.qsub.stderr"
                QSUB = open(qsubFile, "w")
                outline = '#$ -S /bin/bash\n'
                QSUB.write(outline)
                outline = "#$ -l virtual_free=" + memoryNeeded + "\n"
                QSUB.write(outline)
		outline = mapbash_filename + "\n"
                QSUB.write(outline)
		#outline = "/data/software/R/2.15.2/bin/R CMD BATCH --no-restore --no-save " + rscript_filename + " " + processed_dir + "/" + guide["seq.unit.id"][i] + ".cbs.r.out\n"
		outline = "R CMD BATCH --no-restore --no-save " + rscript_filename + " " + processed_dir + "/" + guide["seq.unit.id"][i] + ".cbs.r.out\n"
                QSUB.write(outline)
		#outline = "/data/software/R/2.15.2/bin/R CMD BATCH --no-restore --no-save " + rscript_filename5k + " " + processed_dir + "/" + guide["seq.unit.id"][i] + ".cbs.5k.r.out\n"
		outline = "R CMD BATCH --no-restore --no-save " + rscript_filename5k + " " + processed_dir + "/" + guide["seq.unit.id"][i] + ".cbs.5k.r.out\n"
                QSUB.write(outline)
		#outline = "/data/software/R/2.15.2/bin/R CMD BATCH --no-restore --no-save " + rscript_filename20k + " " + processed_dir + "/" + guide["seq.unit.id"][i] + ".cbs.20k.r.out\n"
		outline = "R CMD BATCH --no-restore --no-save " + rscript_filename20k + " " + processed_dir + "/" + guide["seq.unit.id"][i] + ".cbs.20k.r.out\n"
                QSUB.write(outline)
                QSUB.close()

                time.sleep(1)
                thisCommand = "chmod 755 " + mapbash_filename
                os.system(thisCommand)
                thisCommand = "chmod 755 " + qsubFile
                os.system(thisCommand)
                #thisCommand = "qsub -q all.q@wigclust5 " + qsubFile
                thisCommand = "qsub " + qsubFile
                #thisCommand = qsubFile + " > " + qsubFileStdout + " 2> " + qsubFileStderr + " &"
                os.system(thisCommand)



def fileToGuide(infilename):
	INFILE = open(infilename, "r")
	
	guide = dict()
	x = INFILE.readline()
	colnames = x.rstrip().split("\t")
	for c in colnames:
		c = c.strip()
	ncols = len(colnames)

	a = []
	for x in INFILE:
		arow = x.rstrip().split("\t")
		for b in arow:
			b = b.strip()
		if len(arow) < ncols:
			n = len(arow)
			for i in range(n, ncols):
				arow.append("")
		a.append(arow)

	z = zip(*a)
	if len(z) == ncols:
		for i in range(ncols):
			if guide.has_key(colnames[i]):
				print "ERROR: Duplicate column name."
			else:
				guide[colnames[i]] = z[i]
	else:
		print "ERROR: Data and colnames length not equal."
		
	INFILE.close()
	return guide


if __name__ == "__main__":
	main()
