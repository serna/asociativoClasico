import settings as sets # reads the content of the file config.py
import shutil as shu # to use the function of copy files 
import sys
import os # for operations in the OS
import os.path
import time
import numpy as np
def mkdir(dir_name):
    """ Creates a directory named 'dir_name' """
    try:
        os.makedirs(dir_name)
        return 1
    except OSError:
        if os.path.exists(dir_name):
            #pass
            return 0
        # let exception propagate if we just can't
def createSettings(N,T,RHO,lb,rc,rd,epsilon,nPass,path):
	""" Creates the configuration file for a particular simulation 
    
		This function creates the configuration file "settings.dat" that
		defines the conditions of simulation, the file is created in 
		the directory defined by the argument path.
	"""
	if path == "":
		fo=open("settings.dat","w")
	else:
		if not os.path.exists(path):
			mkdir(path)
			shu.copy('cesar.out',path+'/') # copy the file 
		if not os.path.exists(path+'/'+"cesar.out"):
			shu.copy('cesar.out',path+'/') # copy the file 
		fo=open(path+'/'+"settings.dat","w")
	fo.write("Particles(N): "+str(N))
	fo.write("\nTEMP: "+str(T))
	fo.write("\nRHO: "+str(RHO))
	fo.write("\nlb: "+str(lb))
	fo.write("\nDISPL: 0.5")
	fo.write("\nNMOVE: "+str(nPass*N))
	fo.write("\nNSUB: "+str(int(N*2.0/0.4)))
	fo.write("\nLoad_prev: 0")
	fo.write("\nSource_file: final_conf.dat")
	fo.write("\nDISPLONESITE: 0.5")
	fo.write("\nrc: " + str(rc))
	fo.write("\nrd: " + str(rd))
	fo.write("\nepsilon: " + str(epsilon))
	fo.close()
def dirNames():
	""" Return a list with all the necessary directories for the simulation
	
		This function return a list with all the paths necessary for the set of
	simulations specified in the file settings.py, it will return all the 
	possible combinations between the parameters of the set of simulations
	"""
	dirList = []
	for N in sets.N:
		for T in sets.T:
			for RHO in sets.RHO:
				for lb in sets.lb:
					for rc in sets.rc:
						for rd in sets.rd:
							for epsilon in sets.epsilon:
								path = 'N'+str(N)+'/TEMP'+str(T)+'/RHO'+str(RHO)+'/LB'+str(lb)+'/RC'+str(rc)+'/RD'+str(rd)+'/EPS'+str(epsilon)
								dirList.append(path)
	return dirList
def takeCharAway(string):
	""" Return the numbers of the chain 'string'
		
		This function return only the numbers inside a string in the 
		same order that them appears, the char '-' is considered a number
		in this function
	"""
	lista = []
	for car in string:
		lista.append(car)
	cnt = 0
	for car in lista:
		if car.isdigit()==False and car!="-":
			cnt += 1
		else:
			break
	newCad = ""
	for i in range(cnt,len(string)):
		newCad += string[i]
	return newCad
def breakDir(path):
	""" From the directory given in path it return the corresponding values to set a simulation
		
		Once the directories have been created it is neccesary to fill it with the corresponding bin and 
		settings.dat files, but it is necessary to know what are the corresponding values for the simulation
		these can be read from the name of the directory, so this function does that thing, read the directory 
		the values for the simulation
	"""
	N = 0
	T = 1.0
	RHO = 0.5
	lb = 0.1
	lista = path.split("/")
	onlyNumbers = []
	#print path.split("/")
	for item in lista:
		onlyNumbers.append(takeCharAway(item))
	N = int(onlyNumbers[0])
	T = float(onlyNumbers[1])
	RHO = float(onlyNumbers[2])
	lb = float(onlyNumbers[3])
	rc = float(onlyNumbers[4])
	rd = float(onlyNumbers[5])
	epsilon = float(onlyNumbers[6])
	return N,T,RHO,lb,rc,rd,epsilon
def buildDirectories(lista):
	""" Creates the set of directories in the list 'lista'
	
		This function check if the directory already exist, in that case, nothing is done
	"""
	newDirectories = []
	cnt = 0
	for directory in lista:
		#print directory
		#print breakDir(directory)
		if mkdir(directory) == 1: # the directory has been succesufully created
			newDirectories.append(directory)
			createSettings(*breakDir(directory),nPass = sets.thermalPass, path = directory)
			print "Directory " + directory + " builted"
			myfile=directory+'/cesar.out'
			shu.copy('cesar.out',directory+'/') # copy the file 
			cnt += 1
		else: # the directory already exists
			print "The directory " + directory + " already exist, any change done on this directory"
	print str(cnt) + " new directories created"
def allSubdirs():
	""" Return a list of all the subdirs"""
	dirs = []
	for N in sets.N:
		directory = "N"+str(N)+"/"
		lista = [x[0] for x in os.walk(directory)]
		for item in lista:
			if len(item)>24:
				dirs.append(item+"/")
	return dirs
def printAllSubdirs():
	dirs = []
	for N in sets.N:
		directory = "N"+str(N)+"/"
		lista = [x[0] for x in os.walk(directory)]
		for item in lista:
			if len(item)>24:
				dirs.append(item+"/")
	cnt  = 1
	for item in dirs:
		print cnt,"\t",item
		cnt += 1

def createScripts():
	""" Create scripts to run the set of simulations

		This function updates the 'scripts.sh' each time it is called, when a set of simulations
		are set and for some reason the whole set has not finished it is posibble to run each script
		again and each script will command the uncompleted runs, also if there are other simulation
		not completed or added in a subsecuent time this function will add the corresponding instructions
	"""
	script = []
	for i in range(sets.nProcessors):
		script.append("script"+str(i+1)+".sh")
	for item in script: # create scrips.sh files
		fo = open(item,"w")
		fo.close()
	# first checking the directories defined by ''settings.py'
	allDirs = dirNames()
	lista = [] # save the list of all the pending simulations
	for currentDir in allDirs:
		try:
			ff = open(currentDir+"/status.dat")
		except IOError:
			# if cannot open the status.dat file then this simulation has to be re-scheduled
        		lista.append(currentDir)
		else:
		        # looking for the current progress of the simulation, if it is not 100% then it will be scheduled
		        ff.readline()
			linea = ff.readline()
			if int(linea.split()[0])!=100:
				lista.append(currentDir)
			ff.close()
	for i in range(len(lista)):
		currentScript = i%(len(script))
		fo = open(script[currentScript],"a")
		fo.write("cd " + lista[i] + "\n")
		#########################
		# here is where an improved can be done, by taking the smulation from the last stage it was 
		# and not doing again from the very start
		N,T,RHO,lb,rc,rd,epsilon = breakDir(lista[i])
		createSettings(N,T,RHO,lb,rc,rd,epsilon,sets.thermalPass,lista[i])
		#########################
		fo.write("./cesar.out >> out.txt\n")
		fo.write("./cesar.out " + str(N*sets.averagePass) + " >> out.txt\n")
		n = lista[i].count("/") # count the number of subdirs
		for i in range(n+1):
			fo.write("cd ..\n")
		fo.close()
		# create a script that runs all the scripts
		fo = open("runAll.sh","w")
		n = 1
		for item in script:
			fo.write("./" + item + " >> out" + str(n) + ".txt&\n")
			n+=1
		fo.close()
		# finally change to ejecution mode all the scripts
		for element in script:
			os.chmod(element, 0755) # change to execution mode all the scripts
		os.chmod("runAll.sh", 0755)

		
	#allDirs = allSubdirs():
	# HETE I WANT TO REMOVE FROM ALL SUBDIRS ALL THE PREVIOUS CONSIDERED, IN ORDER TO DO THAT I HAVE TO USE A REMOVE
	# METHOD FOR LIST, AND THEN SPLIT THE TOTAL OF REMAINING SIMULATION IN ALL THE SH FILES
	
def copyFilesTo(lista):
	"""This function copy al the bin files to the correspond
		
	"""
	for directory in lista:
		myfile=directory+'/cesar.out'
		if os.path.isfile(myfile):
			os.remove(myfile)
		else:    ## Show an error ##
			print("Error: %s file to delete not found" % myfile)
		shu.copy('cesar.out',directory+'/') # copy the file 
def getMaximum_gr(fileName):
	ff = open(fileName,"r")
	ff.readline() # ignore the irst line
	x = []
	y = []
	for line in ff:
		x.append(float(line.split()[0]))
		y.append(float(line.split()[1]))
	ff.close()
	index = 0
	maxVal = 0.0
	for i in range(len(x)):
		if maxVal < y[i]:
			index = i
			maxVal = y[i]
	return x[index],y[index]
#printAllSubdirs()
#
if sys.argv[1]=="-build":
	lista = dirNames() # create a list with all the directories of the simulations indicated in settings.py
	#print dirNames()
	buildDirectories(lista) # build all the directories and copy the corresponding files 
	createScripts() # create the scripts to run the set of simulation automatically
elif sys.argv[1] == "-status":
	# Here the program will print the status of the current 
	print "Under construction"
	lista = dirNames() #lista de todos los directorios
	cnt = 1
	for directory in lista:
		if os.path.isfile(directory+"/status.dat"):
			
			ff = open(directory+"/status.dat",'r')
			ff.readline()
			linea = ff.readline()
			avance = linea.split()[0]
			tiempo = linea.split()[2]
			ff.close()
			print cnt,directory,"\t",avance,"\t",tiempo
		else:
			print cnt,directory,"\t","0"	
		cnt += 1

elif sys.argv[1] == "-getFile":
	if len(sys.argv) < 2:
		print "You must indicate the file you want to get from the set of simulations"
		sys.exit()
	else:
		lista = dirNames()
		mkdir("files") # build the directory where the files are going to be copied
		fileName = sys.argv[2]
		for directory in lista:
			copyName = directory.replace("/", "")
			fileToCopy = directory + "/" + fileName 			
			if os.path.isfile(fileToCopy):
				# The file exist
				print "Copying the file " + fileToCopy
				shu.copy(fileToCopy,"files/" + copyName + fileName)
			else:
				print "The file " + fileToCopy + "do not exist"
elif sys.argv[1] == "-restart": # if something goes wrong, like a fail in the electric service, then it run all simulation from where it was before the fail
	createScripts() # create the scripts to run the set of simulation automatically
elif sys.argv[1]=="-getInfo":
	# first checking the directories defined by ''settings.py'
	allDirs = dirNames()
	lista = [] # save the list of all the pending simulations
	gg = open("files/infoSummary.txt","w")
	title = ("N\tT\trho\tlb\trc\t\trd\teps\tmeanXfree")
	gg.write(title)
	print title
	cntGral =1
	missingSim = 0
	for currentDir in allDirs:
		try:
			ff = open(currentDir+"/status.dat","r")
		except IOError:
			print cntGral , "Nothing to print " , currentDir +"/status.dat"	
			# if cannot open the status.dat file then there is any information about this simulation
			cntGral += 1
			missingSim += 1
		else:   # check the current status of the simulation
			ff.readline()
			line = ff.readline()
			doneHere = line.split()[0]
			computeUncertainty = True
			if int(line.split()[0])!=100:
				lista.append(currentDir)
				missingSim += 1
				computeUncertainty = False
			ff.close()
			ff=open(currentDir+"/out.txt","r")
			fechaMod = os.path.getmtime(os.path.join(currentDir, "out.txt"))
			fechaMod = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(fechaMod))
			lineFinal1 = "" # the line previous to the final line
			lineFinal2 = "" # the line previous to lineFinal1
			finalLine = "" # the last line
			cnt = 0
			for line in ff:
				lineFinal2 = lineFinal1
				lineFinal1 = finalLine
				finalLine = line
				cnt+=1
			ff.close()
			N,T,RHO,lb, rc, rd, eps = breakDir(currentDir)
			if len(lineFinal2.split())>7:
				if computeUncertainty==True:
					ff=open(currentDir+"/out.txt","r")
					cnt1 = 0
					sumx = 0.0
					sumx2 = 0.0
					for line in ff:
						cnt1+=1
						if cnt1>=cnt-101 and cnt1<cnt-1:
							#print line
							val = float(line.split()[7])
							sumx += val
							sumx2 += val*val
					sigma = np.sqrt(sumx2/100.0-sumx*sumx/10000.0)
					ff.close()
					info = "\n{0}\t{1}\t{2}\t{3}\t{4}  \t{5}\t{6}\t{7}\t{8}".format(N,T,RHO,lb, rc, rd, eps,lineFinal2.split()[7],sigma)
					
				else:
					info = "\n{0}\t{1}\t{2}\t{3}\t{4}  \t{5}\t{6}\t{7}".format(N,T,RHO,lb, rc, rd, eps,lineFinal2.split()[7])
				print cntGral,doneHere, info[1:],fechaMod
				gg.write(info)
			cntGral +=1
			

	gg.close()
	print "There are missing ",missingSim , " simulation"
	timePerSimulation = 2 # this time has to be in hours
	print "Assuming 2 hours for each simulation, estimated total time (hours): ", round(timePerSimulation*missingSim/sets.nProcessors,2), ", days:", round(timePerSimulation*missingSim/24.0/sets.nProcessors,2)
elif sys.argv[1] == "-getSummary_gr":
	lista = dirNames() # get the list of all the directories where the study will be done
	ff = open("files/summary_gr.dat","w")
	print "Creating file summary_gr.dat"
	for directory in lista:
		fileName = "files/" + directory.replace("/","")
		fileName += "radial.dat"
		N,T,RHO,lb,rc,rd,epsilon= breakDir(directory) # get the configuration of the curret simulation 
		if os.path.isfile(fileName): # if the file exist then compute the maximum of the radial distribution function
			print "Getting data form file " + fileName
			x,y = getMaximum_gr(fileName)
			ff.write(str(RHO) + "\t" + str(lb) + "\t" + str(x) + "\t" + str(y) + "\n")
	ff.close()
else:
	print "You have to type an option!"		
