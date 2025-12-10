#!/usr/bin/env python
# -*- coding: utf-8 -*-

import getopt, os
#import platform, math

import sys

def handleOutputs(output):
    
    if not output:
        outputDir = 'results'
        outputPath = os.path.join(os.getcwd(), outputDir)
    else:
        try:
            output.split('/')[-2]
        except IndexError:
            outputPath = os.path.join(os.getcwd(), output)
        else:
            outputPath = output

    return outputPath    


def parseCommandLineOptions(argv):
	# parse command line options
    try:
        opts, args = getopt.getopt(argv[1:], 'hu', ['help', 'usage', 's=', 'l=', 'o=', 'm=', 'db='])
    except(getopt.GetoptError, err):
        print('\n'+str(err)+'\nFor help use -h or --help\n')
        sys.exit()

    if opts != []:

        output = None
        inputType = None
        pathToInput = None
        metal = None
        db_path = False
        
		# process options
        for o, a in opts:
            if o == '--s':
                inputType = 'sequence'
                pathToInput = a
            elif o == '--l':
                inputType = 'list'
                pathToInput = a
            elif o == '--o':
                output = a		
            elif o == '--m':
                metal = a
            elif o == '--db':
                db_path = a
                if not os.path.exists(db_path):
                   print('\n%s file not exists.\nFor help use -h or --help\n'%str(db_path))
                   sys.exit()
                
            elif o in ('-u', '--usage'):
		# print usage information and exit
                usage()
                sys.exit()
            elif o in ("-h", "--help"):
		# print help information and exit
                printHelpInfo()
                sys.exit()
            else:
                assert False, "unhandled option"
                
        if not pathToInput:
            print('\nInput parameter for input file is mandatory.\nFor help use -h or --help\n')
            sys.exit()
        if not metal:
            print('\nInput parameter for metal is mandatory.\nFor help use -h or --help\n')
            sys.exit()

        optList = []
        for opt in opts:
            optList.append(opt[0])
            
        if (('--s' in optList) and ('--l' in optList)):
            print('\nThe program accept as input data a fasta file OR a directory which contain fasta files.\nFor help use -h or --help\n')
            sys.exit()
            
        if ('--s' in optList) and (not opened(pathToInput)):
            if os.path.isdir(pathToInput):
                print('\n%s is a directory, use --l to give a directory as input data\n\nFor help use -h or --help\n'%pathToInput)
            else:
                print('\nFile %s do not exist or not accessible\n'%pathToInput)
                sys.exit()
            
        if ('--l' in optList) and (not diropen(pathToInput)):
            if os.path.isfile(pathToInput):
                print('\n%s is a file, use --s to give a file as input data\n\nFor help use -h or --help\n'%pathToInput)
            sys.exit()      
            
        if not output:
            name_out1 = pathToInput.split("/")[-1]
            name_out = name_out1.split(".")[0]
            outputDir = 'results_%s_%s'%(metal, name_out)
            output = os.path.join(os.getcwd(), outputDir)
    
    else:
        print("\nInput parameters for input file is mandatory.\nFor help use -h or --help\n")
        sys.exit()

    return inputType, metal, pathToInput, output, db_path


def printHelpInfo():
	summary()
	usage()
	example()
	helpinfo()	
	moreinfo()

	
def summary():
	print('\nSummary:\n--------')
	print('MetalPredator is a program aimed at predict iron-sulfur cluster binding proteins from protein sequence(s).\nMetalPredator is able to process complete proteomes rapidly with high recall and precision.\n')

	
def usage():
	print('\nUsage:\n------')
	print('$./metalPredator.py [input parameter] <file OR directory> [output options] <output directory>\n')

def example():
	print('\nExamples:\n------')
	print('$./metalPredator.py --s ./data/Escherichia_coli_str_K_12.fasta\n')
	print('or\n')
	print('$./metalPredator.py --l ./data/fasta_dir\n')

def helpinfo():

	print('{0}{1}'.format('Options:\n','--------'))
	print('{0:<22}{1:<20}{2}{3:<22}{4:<20}{5}{6:<22}{7:<20}{8}{9:<22}{10:<20}{11}{12:<22}{13:<20}{14}{15:<22}{16:<20}{17}{18:<22}{19:<20}{20}'.format(
			'   --s  <file>',  'input parameter', 'specify the path to a fasta file\n',
			'   --l  <directory>',  'input parameter', 'specify the path to a directory wich contain fasta files\n',
			'   --m  <name>',  'input parameter', 'specify the name of metal (eg:iron_ion)\n',
			'   --o  <number>', 'output option', 'specify the directory to save output files\n',
                        '   --db <directory>', 'input parameter', 'specify the path to a directory wich contain hammer databases\n',                            
			'   -h   --help',     'flag', 'print help information\n', 
			'   -u   --usage',    'flag', 'print usage summary'))


def moreinfo():
	print('\nMore info:\n----------')
	print('Input parameter is mandatory.\nInput options, options and output directory are non-mandatory.\nIf output directory is not supplied, results will be stored in a current working directory.\n')
 

def opened(fileName):
	try:
		f = open(fileName)
		f.close()
		return True    
	except IOError:
		return False
    
    
def diropen(dirName):
    
    flag = False
    if os.path.isdir(dirName):
        for file_name in os.listdir(dirName):
            flag = True
            file_path = os.path.join(dirName, file_name)
            if opened(file_path):
                return True
            else:
                print('\nFiles in directory %s do not accessible\n'%dirName)
                return False
            
        if not flag:
            print('\nDirectory %s is empty\n'%dirName)
            return False
    else:
        print('\nDirectory %s do not exist\n'%dirName)
        return False


def askToOverwriteIt():
	
	answer = raw_input()

	if answer == 'y':
		flag = True
	elif answer == 'n':
		print('\nThe process was aborted.\n')
		flag = False
	else:
		print('\nEnter "y" or "n".\n')
		flag = False
		askToOverwriteIt()
	
	return flag
