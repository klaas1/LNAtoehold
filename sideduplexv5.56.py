
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib as mpl
mpl.rcParams['mathtext.default'] = 'regular'
#mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sequtils as seq
#np.set_printoptions(precision=2)

#Boolean, export images to public_html (False) or local folder (True)
local = True
if local:
	path = "."
else:
	path = "../../public_html"

#Gas constant in kcal/K/mol
R = 1.987204118e-3

#5'-3' Sequence of the probe, LNA bases indicated by +N
#probe ="G+CGGC+CC+ACC+TG+CTGGT+A+CG" #human FSDH-KpnI probe
#probe = "T+A+TA+GG+GAATA+TT+AAGCT" #LNA probe w/o extra
probe = "T+A+TA+GG+GAATA+TT+AAGCTG" #LNA probe
#probe = "G+CA+C+A+T+ATA+CACC+ATGC" #LNA probe LINE1 NcoI
#probe = "T+TGGAG+T+T+GC+TC+T+TC+TCGAC" #LNA probe LINE1 XhoI
#probe = "T+A+TA+GG+GAATA+TT+AAGCT" #LNA probe w/o 3'stack
#probe = "T+A+TA+GG+GAATA+TT+A" #LNA probe w/o toehold
#probe = "TATAGGGAATATTAAGCT" #LNA probe without LNA
#probe = "CGGCAA+GCCGC+A+TGGC" #Sideduplex probe
#probe = "GGGCGGCGACCTC" #Lambda
#probe = "GGGCGGCGACCT" #Lambda w/o extra bases
#probe = "GCCGCGGCCGCCCG"
#probe = "G+GT+CGTT+CG+CTC+CA+AGCT+GGG" #Arthur probe BseYI
#probe = "G+TG+GG+TC+TCGC+GGTG" #Arthur probe BsaI
#probe = "C+AA+TT+TG+TGG+AA+TTC+TCGAC" #LacO Repeat XhoI
#probe = "A+TTT+G+TGGAA+TTC+TCGAC" #LacO Repeat XhoI
#probe = "G+AGTCG+ACGC+ATGC+A+AGCT" #LacO repeat hindIII
#probe = "T+GGAA+TTC+TC+GA+C+TC+TAGC" #lacO XbaI
#read probe from txt file
#probelist = open('probe.txt', 'r')
#probe = probelist.readline() # reads first probe, "for line in probelist" reads all
#3'-5' Complementary sequence of the target with toehold
#for non-complementary sequences, define target manually
target = ""
for i in probe:  
    if(i == "T"):
        target = target + "A"
    if(i == "A"):
        target = target + "T"
    if(i == "G"):
        target = target + "C"
    if(i == "C"):
        target = target + "G"
bp=len(target)        
#target = "GCCGAAACGGCGTAC" #override autotarget for mismatches
#target = "ATATCGCTTATAATTCGACG"
#print 'target =', target
print(bp,'bp')
#Number of 5'-toehold bases in the target
toehold = 4
probeextra = 1 #Are there bases to stack to (if probe is a hairpin, probably yes) Add the last base of the hairpin to the probe sequence
if probeextra:
    target = target[:-probeextra]
    bp = bp-probeextra

parameters = np.genfromtxt('NNparamsv3.csv',delimiter=",",dtype=None,skip_header=1) #loads parameters
NNparameters = {label[0]: [float(label[1]),float(label[2]),float(label[3])] for label in parameters}

#Function to calculate the deltaG at 37 degrees C
def deltaG37(seq,seq2,init):
	prevc = ''
	prevd = ''
	initc = ''
	if len(seq) > 1:
		initc = seq[0]	
	terms = []
	deltag = 0
	for c,d in zip(seq,seq2):
		if prevc != '':
			terms.append(prevc+c+"/"+prevd+d)
		prevc = c
		prevd = d
	for t in terms:
		print( t)
		#print NNparameters[t][2]
		deltag = deltag + NNparameters[bytes(t, encoding='utf-8')][2]
		
	#Initiation terms for the deltaG
	#WARNING: this is not correct for the LNA bases
	if init:
		for prevc in [prevc,initc] :
			if prevc == 'C' or prevc == 'G' or prevc == "+C" or prevc == "+G":
				deltag += 0.98
				#print "CG initiation: ", 0.98
			if prevc == 'A' or prevc =='T' or prevc == "+A" or prevc == "+T":
				deltag += 1.03
					#print "AT initiation: ", 1.03	
	return deltag

#Function to calculate the entropy at 37 degrees C
def entropy(seq,seq2,init):
	prevc = ''
	prevd = ''
	initc = ''
	if len(seq) > 1:
		initc = seq[0]	
	terms = []
	entropy=0
	for c,d in zip(seq,seq2):
		if prevc != '':
			terms.append(prevc+c+"/"+prevd+d)
		prevc = c
		prevd = d
	for t in terms:
	#	print t
	#	print NNparameters[t][1]
		entropy = entropy + NNparameters[bytes(t, encoding='utf-8')][1]
					
#	for t in terms:
#		entropy= entropy + NNparameters[t][1]
#		#print t, 
#		#print NNentropy[t]
#
	#Initiation terms for the entropy
	#WARNING: this is not correct for the LNA bases 
	if init:

		for prevc in [prevc,initc] :
			if prevc == 'C' or prevc == 'G' or prevc == "+C" or prevc == "+G":
				entropy += -2.8
				#	print "CG initiation: ", -2.8
			if prevc == 'A' or prevc =='T' or prevc == "+A" or prevc == "+T":
				entropy += 4.1
					#		print "AT initiation: ", 4.1
	return entropy

#Function to calculate the enthalpy at 37 degrees C
def enthalpy(seq,seq2,init):
	prevc = ''
	prevd = ''
	initc = ''
	if len(seq) > 1:
		initc = seq[0]	
	terms = []
	enthalpy=0
	for c,d in zip(seq,seq2):
		if prevc != '':
			terms.append(prevc+c+"/"+prevd+d)
		prevc = c
		prevd = d
	for t in terms:
		#print t
		#print NNparameters[t][0]
		enthalpy = enthalpy + NNparameters[bytes(t, encoding='utf-8')][0]
				
	#Initiation terms for the entropy
	#WARNING: this is not correct for the LNA bases 
	if init:
		for prevc in [prevc,initc] :
			if prevc == 'C' or prevc == 'G' or prevc == "+C" or prevc == "+G":
				enthalpy += 0.1
			#	print "CG initiation: ", 0.1
			if prevc == 'A' or prevc =='T' or prevc == "+A" or prevc == "+T":
				enthalpy += 2.3
			#	print "AT initiation: ", 2.3	
	return enthalpy


#Function to find possible matches between two sequences and returns the
#lowest deltaG
def findMatches(seq1,seq2):
	dG = 0	
	seqN = []
	cB = 0
	cC = 0
	matches = []
	#Loop over all starting positions in seq1 for a match
	for i,c in enumerate(seq1):	
		#Loop over all lengths of matches
		for j in np.arange(1,len(seq1)):
			#Loop over all starting positions of the match in seq2
		#	for k,e in enumerate(seq2):
                 for k in np.arange(0,len(seq2)):		
            		#print "\n"
				#print seq1[i:i+j]
				#print seq2[k:k+i+j]
				#print "\n"
                     comp = seq.checkComp(seq1[i:i+j],seq2[k:k+j])
                     if comp:# and j+k >4:
#TODO:				Check for bulge possibility
                         dGtemp = deltaG37(seq1[i:i+j],seq2[k:k+j],True)
					#print "Yeah, ",dGtemp, " ", seq1[i:i+j], " ", seq2[k:k+j]
                         matches.append(dGtemp)
                         if dGtemp < dG:
                             dG = dGtemp
                             #print "DG:", dG," ", seq1[i:i+j]
                             seqN = seq1[i:i+j]
                             cB = i
					#cE = j
					#seqT = seq2[k:k+i+j]
                             cC = k
	return dG,seqN,cB,cC,matches

def calcKd(dG):
	return 1.0/np.exp(-dG/(R*(273.15+37)))

def plotStructure(seq1,seq2,seq3,i,j,k):
    plt.clf()
    color = 'black'
    plt.text(0,0.35,"Target: "+str(j)+" " + " " +str(k) + " "+ str(len(seq3)))
    plt.text(0,0.3,"Bulge:"+''.join(seq2[-i:][:-toehold][:-k-len(seq3)]))
    #plt.text(0,0.2,"Hybrid: "+''.join(seq.makeComp(seq3[::-1]))+ "      dG: " + str(deltaG37(seq3)))
    if toehold > 0:
        plt.text(0,0.1,"End: "+''.join(seq2[-i:][:-toehold][-k:]))
    #print seq1
    #print seq2
    #print seq2[-i:]
    #print seq2[-i:][:-toehold]
    #print seq2[-i:][:-toehold][:-k-len(seq3)]
    plt.text(0.8,0.35,"Probe")
    plt.text(.8,0.3,''.join(seq1[:-i][j+len(seq3):]))#Bulge
    plt.text(.8,0.2,''.join(seq3))#Hybrid
    plt.text(.8,0.1,''.join(seq1[:-i][:j]))#End
    #print seq1
    #print seq1[:-i]
    #print seq1[:-i][j+len(seq3):]
    #Plot doublestranded part
    for n,c in enumerate(seq2[:-i]+seq1[-i:]):
        if n > len(seq2[:-i])-1:
            color="blue"
        if c.startswith("+"):
            plt.text(0.05*n,0.5,c[1:],color="red")
        else:
            if not (n > len(seq2[:-toehold])-1 and n < len(seq2)-i):
                plt.text(0.05*n,0.5,c,color=color)
        plt.text(0.05*n,0.46,seq.makeComp([c])[0],color="black")
        #Plot target tail
    for o,c in enumerate(seq2[-i:-toehold]):
        if c.startswith("+"):
            plt.text(0.05*(n-i),0.55+0.04*o,c[1:],color="red")
        else:
            plt.text(0.05*(n-i),0.55+0.04*o,c,color='black')
        #Plot probe tail
    for p,c in enumerate(seq1[:-i][::-1]):
        if c.startswith("+"):
                plt.text(0.05*(n-i+1),0.55+0.04*p,c[1:],color="red")
        else:
                plt.text(0.05*(n-i+1),0.55+0.04*p,c,color='blue')
    plt.axis('off')
    plt.savefig(path+"/structure_"+str(i)+".png")
		
#Convert the probe and target sequences to list format
probe = seq.toList(probe)
#extra = []
#if probeextra != 0:
#        extra = probe[-probeextra:]
#	probe = probe[:-probeextra]

target = seq.toList(target)

Gibbs_Sideduplex = []
GibbsN = []
GibbsO = [] # The other sideduplex energies
target_energies = []
probe_energies = []
target_entropies = []
probe_entropies = []
target_enthalpies = []
probe_enthalpies = []

x = []

print('\n')
print("========================================")
print("Probe: ", probe)
print(" ")
print( "Target: ", seq.makeComp(target))
print( " ")
print( "========================================" )
print( '\n')

toeholdZero = ( toehold == 0 )

probe_invading = seq.toList(probe)
extra_seq = []
if probeextra:
    extra_seq = seq.toList(probe[-probeextra:])
    probe_invading = seq.toList(probe[:-probeextra])

toehold_seq = []
invaded_target = seq.toList(target)
if toehold:
    toehold_seq = seq.toList(target[-toehold:])
    invaded_target = seq.toList(target[:-toehold])

print( "Probe " , probe_invading, " ", extra_seq)
print( "Target " , invaded_target, " ", toehold_seq)

#Loop to do the strand invasion
for i in np.arange(0,len(probe_invading)+1):

	energy_probe=deltaG37(extra_seq, seq.makeComp(extra_seq),True)
	entropy_probe =entropy(extra_seq, seq.makeComp(extra_seq),True)
	enthalpy_probe =enthalpy(extra_seq, seq.makeComp(extra_seq),True)	

	#energy_probe= deltaG37(extra_seq, seq.makeComp(extra_seq))
	#entropy_probe = entropy(extra_seq, seq.makeComp(extra_seq))
	#enthalpy_probe = enthalpy(extra_seq, seq.makeComp(extra_seq))	

	#print "Probe "  , i
	basesinc = 0
	if i > 0:
		probeseq = probe_invading[-i:]+extra_seq
		probeseq_comp = (invaded_target+toehold_seq)[-i:]+seq.makeComp(extra_seq)
		energy_probe = deltaG37(probeseq,probeseq_comp,False)
		entropy_probe= entropy(probeseq,probeseq_comp,False) 
		enthalpy_probe= enthalpy(probeseq,probeseq_comp,False)

		basesinc = len(probe_invading[-i:])	

	probe_energies.append(energy_probe)
	probe_entropies.append(entropy_probe)
	probe_enthalpies.append(enthalpy_probe)

	#print "END PROBE", "\n"

	#print "Target"
	if i > toehold:
		energy_target = deltaG37(seq.makeComp(invaded_target[:-(i-toehold)]),invaded_target[:-(i-toehold)],False)
		entropy_target = entropy(seq.makeComp(invaded_target[:-(i-toehold)]),invaded_target[:-(i-toehold)],False)
		enthalpy_target = enthalpy(seq.makeComp(invaded_target[:-(i-toehold)]),invaded_target[:-(i-toehold)],False)
	else:
		energy_target = deltaG37(seq.makeComp(invaded_target),invaded_target,False)+ deltaG37(seq.makeComp(invaded_target[-1]),invaded_target[-1],True)
		#print "OOHH ", deltaG37(seq.makeComp(invaded_target[-1]),invaded_target[-1],True)
		entropy_target = entropy(seq.makeComp(invaded_target),invaded_target,False)+ entropy(seq.makeComp(invaded_target[-1]),invaded_target[-1],True)
		enthalpy_target = enthalpy(seq.makeComp(invaded_target),invaded_target,False)+ enthalpy(seq.makeComp(invaded_target[-1]),invaded_target[-1],True)

	target_energies.append(energy_target)
	target_entropies.append(entropy_target)
	target_enthalpies.append(enthalpy_target)

	#print "END TARGET","\n"	
	sideduplex = 0
	sd1 = []
	sd2 = []
	j= k = 0
	matches = []	
	if i > toehold:	
		sideduplex,sd1,j,k,matches = findMatches(probe_invading[:-i],seq.reverse(seq.makeComp(target[-i:])))	
		
	total_energy = energy_probe+energy_target+sideduplex
	total_entropy = entropy_probe+entropy_target
	total_enthalpy = enthalpy_probe+enthalpy_target
	
	Gibbs_Sideduplex.append(total_energy)
	GibbsN.append(energy_probe+energy_target)
	GibbsO.append(matches)
	#print "Free energy (",basesinc," bases incorporated): ", 
	#print energy , "kcal/mol"
	#print "Free energy (no sideduplex): ",
	#print pr+tr
	#print '\n'
	#if energy < pr+tr:		
	plotStructure(probe,seq.makeComp(target),sd1,i,j,k)
	x.append(i)
	#raw_input("Press Enter to continue...")

###ENTROPY#####
print(total_enthalpy)
total_enthalpy = total_enthalpy-(probe_enthalpies[0]+target_enthalpies[0])
print( total_enthalpy)
print( total_entropy)
total_entropy = total_entropy-(probe_entropies[0]+target_entropies[0])
print( total_entropy)
#print "Tm: ", total_enthalpy/(total_entropy/1000.0+R*np.log(1.0e-6))-273.15
#Tm = total_enthalpy/(total_entropy/1000.0+R*np.log(0.156e-9))-273.15
#Tm = total_enthalpy/((total_entropy/1000.0)+R*np.log(0.25e-6 - 0.25e-6/10.0))-273.15
#Tm = total_enthalpy/((total_entropy/1000.0)+R*np.log(0.25e-6))-273.15
#Tm = total_enthalpy/((total_entropy/1000.0)+R*np.log(6e-13/4.0))-273.15
Tm = total_enthalpy/((total_entropy/1000.0)+R*np.log(0.370e-9))-273.15
#Tm = total_enthalpy/((total_entropy/1000.0)+R*np.log(0.25e-6))-273.15
print( "Tm: ", Tm)
Tm2 = 1.0/(1.0/Tm + ((4.29*12.0/14.0 - 3.95)*np.log(2) + 0.940 * (np.log(2))**2.0)* 1.0e-5)
print( "Tm salt correction: ", Tm2)

print( "Gibbs Free Energy (H,S): ", total_enthalpy-(273.15+37.0)*total_entropy/1000.0, "kcal/mol")
###ENTROPY#####
	
#Calculating kD
print( '\n', "########################")
print( "Calculating Kd: ")
energy_target = GibbsN[0]#0

for objects in GibbsN:
	GibbsN.insert(0,GibbsN.pop()-energy_target) #removes target energy, changing G to deltaG

#if toeholdZero:
#energy_target = deltaG37(seq.makeComp(invaded_target), invaded_target,True) -1.03
#energy_target += deltaG37(extra_seq,seq.makeComp(extra_seq))
#else:	
#	energy_target = deltaG37(seq.makeComp(invaded_target[:-toehold]))
#	energy_target += deltaG37(extra_seq)
print( "Target: ", energy_target, " kcal/mol")
print( "Target+ Probe: ", total_energy, " kcal/mol")
print("Gibbs Free Energy (G): ",total_energy-energy_target, "kcal/mol")
Kdcalc = calcKd(-np.abs(energy_target-total_energy))
print("Kd (no SideDuplex)", Kdcalc/1.0e-9, " nM")
#print "Breaking force: ", (energy_target-total_energy)/(bp*0.332*0.1439), " pN"
print( "Breaking force: ", (energy_target-total_energy)/(bp*(1500.0/4650.0)*0.1439), " pN")

GibbsX = []
GibbsY = []
for i,t in enumerate(GibbsO):
	for j in t:
		GibbsX.append(i)
		GibbsY.append(j+GibbsN[i])

#Calculate Gibbs free energy for SideDuplex
#print np.min(GibbsY)
#print energy
Kdmin = Kdcalc
if GibbsY:
    Kdmin = calcKd(-np.abs(np.min(GibbsY)-energy_target))/1.0e-9
print( "Kd lowest energy: ", Kdmin)


#Plotting
plt.clf()
#plt.plot(x,Gibbs_Sideduplex,marker="o",label="With sideduplex")
plt.plot(x,GibbsN,marker="o",c='red',label="Hybrid")
plt.scatter(GibbsX,GibbsY,c='grey',label="Alternative Structures")
plt.legend()
plt.title("Kd: "+str(np.around(Kdcalc/1.0e-9, decimals=1))+ " nM"+" Kdmin:" + str(np.around(Kdmin, decimals=1))+ " nM")
plt.xticks(np.arange(-1,(bp+1),1))
plt.xlim(-0.5,bp+1)
plt.xlabel("Bases incorporated")
plt.ylabel("Gibbs free energy (kcal/mol)")
plt.grid(True,axis='x')
ax = plt.gca()
ax2 = ax.twinx()
ax2.set_ylim((ax.get_ylim()[0]*1.623,ax.get_ylim()[1]*1.623))
ax2.set_ylabel("Gibbs free energy (kT) T=37C ")
ax2.set_ylabel("Gibbs free energy (kT)")
plt.savefig(path+"/gibbs.png")

print(GibbsN[12])
#Plotting
plt.clf()
plt.plot(x,Gibbs_Sideduplex,marker="o",label="With sideduplex")
plt.plot(x,np.array(target_entropies)+np.array(probe_entropies),marker="o",c='red',label="Entropy")
plt.plot(x,np.array(target_enthalpies)+np.array(probe_enthalpies),marker="o",c='blue',label="Enthalpy")
plt.scatter(GibbsX,GibbsY,c='yellow',label="Sideduplex")
plt.legend()
plt.title("Kd: "+str(np.around(Kdcalc/1.0e-9, decimals=1))+ " nM"+" Kdmin:" + str(np.around(Kdmin, decimals=1))+ " nM")
plt.xticks(np.arange(-1,(bp+1),1))
plt.xlim(-0.5,bp+1)
plt.xlabel("Bases incorporated")
plt.ylabel("Gibbs free energy (kcal/mol)")
plt.grid(True,axis='x')
ax = plt.gca()
ax2 = ax.twinx()
ax2.set_ylim((ax.get_ylim()[0]*1.623,ax.get_ylim()[1]*1.623))
ax2.set_ylabel("Gibbs free energy (kT) T=37C ")
ax2.set_ylabel("Gibbs free energy (kT)")
plt.savefig(path+"/ent.png")
#
plt.clf()
plt.scatter(x,probe_energies,c='green',label="probe FE")
plt.scatter(x,target_energies,c="yellow",label="target FE")
plt.legend()
plt.xticks(np.arange(-1,(bp+1),1))
plt.xlabel("Bases incorporated")
plt.ylabel("Gibbs free energy (kcal/mol)")
plt.savefig(path+"/gibbs2.png")

#Plotting Boltzmann
plt.clf()

GibbsY = np.array(GibbsY) - GibbsN[0]
GibbsN = np.array(GibbsN) - GibbsN[0]
Gibbs_Sideduplex = np.array(Gibbs_Sideduplex) - Gibbs_Sideduplex[0]

boltzN = np.exp(-GibbsN/(R*(273.15+37)))/np.sum(np.exp(-GibbsN/(R*(273.15+37))))
boltz = np.exp(-Gibbs_Sideduplex/(R*(273.15+37)))/np.sum(np.exp(-Gibbs_Sideduplex/(R*(273.15+37))))

GibbsTotal = np.concatenate((GibbsN,GibbsY))
boltzTotal = np.exp(-GibbsTotal/(R*(273.15+37)))/np.sum(np.exp(-GibbsTotal/(R*(273.15+37))))

for i,b in enumerate(boltzTotal[len(probe_invading)+1:]): 
	boltzTotal[GibbsX[i]] += b

#plt.plot(x,boltz,marker="o",label="With sideduplex", color='y')
#plt.plot(x,boltzN,marker="o",c='red',label="No sideduplex")
#plt.scatter(x+GibbsX,boltzTotal)
#plt.plot(x,boltzTotal[:len(probe_invading)+1],marker="o",c="y",label="Sideduplex")
b2 = plt.bar(x, boltzTotal[:len(probe_invading)+1],   1, color='grey', label="Naive Hybrid",align="center")
b1 = plt.bar(x, boltzN,   0.7, color='black', label="Alternative structures included",align="center")

plt.legend([b1[0],b2[0]],["Hybrid", "Alternatives included"],loc="upper left", fontsize = "medium")
plt.title("Kd: "+str(np.around(Kdcalc/1.0e-9, decimals = 1))+ " nM"+" Kdmin:" + str(np.around(Kdmin, decimals=1)) + " nM")
plt.xticks(np.arange(0,(bp+1),1))
plt.xlim(-0.5,(bp+1))
plt.yticks(np.arange(0,1,0.1))
plt.yscale("log")
plt.xlabel("Bases incorporated")
plt.ylabel("Probability")
plt.savefig(path+"/boltz.png")

results = open(path+"/results.html",'w')
results.write("<html>")
results.write("<head>")
results.write("</head>")
results.write("<body>")
results.write("<H1>Probe:</H1>")
results.write(''.join(probe))
results.write("<H1>Target:</H1>")
results.write("Toehold: " + str(toehold))
results.write("</br>")
results.write(''.join(target))
results.write("</br>")

results.write("<a href='./gibbs.png'><img src='./gibbs.png' alt='Gibbs Free Energy plot' style='width:608px;height:456px'></a>")
results.write("<a href='./boltz.png'><img src='./boltz.png' alt='Boltzmann plot' style='width:608px;height:456px'></a>")

results.write("<table style='width:100%'>")
results.write("<tr>")
results.write("<tr><th>Bases incorporated</th><th>Probe</th><th>Gibbs Free energy change</th></tr>")

for i,en in enumerate(GibbsN):
        results.write("<tr>")
        results.write("<td>"+str(i)+"</td>")
        results.write("<td><a href='./structure_"+str(i)+".png'>"+"Structure "+str(i)+"</a></td>")
        results.write("<td><a href='./structure_"+str(i)+".png'><img src='./structure_"+str(i)+".png' alt='Structure' style='width:304px;height:228px'></a></td>")
        results.write("<td>"+str(en)+"</td>")
        results.write("</tr>")
results.write("</table>")

results.write("</body>")
results.write("</html>")

#print boltzTotal[:19]
