from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
plt.rc('font', family='serif', serif='Times')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=12, direction='in')
plt.rc('ytick', labelsize=12, direction='in')
plt.rc('axes', labelsize=12)
###############################################################

##############################
### SELECT PLOT
##############################
GvsH          = False
PBESOL        = False
fracH         = True
BEEFexpansion = False
###############################################################

# width was measured in inkscape
width = 8
height = width / 1.618

colorWheel = ['r','rebeccapurple','mediumslateblue','b']

###############################################################
# G vs H figures
###############################################################
if GvsH:
	### DATA
	xPd = np.array([0,0.25,0.5,0.75,1,1.125,1.25,1.375,1.5,1.625,1.75,1.875,2])
	xPdAu=np.array([0,0.5,1,1.5])

	                 #0 V   -0.3V   -0.6V   -0.9V
	yPd = np.array([[0,		0,		0,		0]
					,[-0.20,	-0.80,	-1.40,	-2.00]
					,[-0.34,	-1.54,	-2.74,	-3.94]
					,[-0.33,	-2.13,	-3.93,	-5.73]
					,[0.36,	-2.04,	-4.44,	-6.84]
					,[1.03,	-1.67,	-4.37,	-7.07]
					,[1.70,	-1.30,	-4.30,	-7.30]
					,[2.48,	-0.82,	-4.12,	-7.42]
					,[2.89,	-0.71,	-4.31,	-7.91]
					,[5.04,	1.14,	-2.76,	-6.66]
					,[6.04,	1.84,	-2.36,	-6.56]
					,[8.72,	4.22,	-0.28,	-4.78]
					,[11.86,7.06,	2.26,	-2.54]])




	yPdAu=np.array([[0.00,	0.00,	0.00,	0.00],
					[0.32,	0.02,	-0.28,	-0.58],
					[0.71,	0.11,	-0.49,	-1.09],
					[1.86,	0.96,	0.06,	-0.84]])
	### FIGURE

	labels=["  0 V","-0.3V",'-0.6V','-0.9V']

	fig, (ax1,ax2) = plt.subplots(nrows=1,ncols=2)
	#fig.subplots_adjust(left=.2, bottom=.2, right=.2, top=.2)
	fig.subplots_adjust(hspace=0.6,wspace=0.4)

	# Combine into list
	ax   = [ax1,ax2]
	x    = [xPd,xPdAu]
	y    = [yPd,yPdAu]
	title = ['Pd hydride stability','PdAu hydride stability']

	### Plot
	for a in range(2):
		for i in range(4):
			c = colorWheel[i]
			xInterp = np.linspace(0,max(x[a]),num=30,endpoint=True)
			f = interp1d(x[a], y[a][:,i], kind='cubic')
			ax[a].scatter(x[a],y[a][:,i],s=50,color=c,marker='.',linewidth='1')
			ax[a].plot(xInterp,f(xInterp),color=c,linestyle='-',label=labels[i])

		ax[a].set_ylabel('$\Delta G_{form}$ (eV)')
		ax[a].set_xlabel('H : Pd')
		ax[a].set_xlim(0, max(x[a]))
		ax[a].set_title(title[a])
		legend = ax[a].legend(loc='upper left', shadow=False,fontsize=12)
		legend.get_frame().set_facecolor('#FFFFFF')



	fig.set_size_inches(width, height)
	fig.savefig('GvsH.pdf')
	plt.show()


########################################################################
########################################################################
if PBESOL:
	### DATA
	exptPdAu=np.array([[1,	4.0782],
			[0.9129,	4.0618],
			[0.754,	4.0318],
			[0.668,	4.0156],
			[0.6176,	4.0061],
			[0.5576	,3.9948],
			[0.4476	,3.974],
			[0.3507	,3.9557],
			[0.2647	,3.9395],
			[0.188,	3.9251],
			[0.1526,	3.9184],
			[0.119,	3.912],
			[0.087,	3.906],
			[0.0566,	3.9003],
			[0,	3.8896]])
				# From Maeland and Flanagan, Canadian Journal of Physics 42 (1964)

	exptPdAux = exptPdAu[:,0];exptPdAuy = exptPdAu[:,1];
	theoryPdAux=[0,   0.25,0.375,0.5, 0.625,0.75,1]
	theoryPdAuy=[3.89,3.94, 3.96,3.99, 4.01,4.04,4.09]

	### FIGURE
	fig, (ax) = plt.subplots(nrows=1,ncols=1)

	ax.scatter(theoryPdAux,theoryPdAuy,s=50,color='b',marker='.',linewidth='1',label='PBE-SOL')
	ax.plot(exptPdAu[:,0],exptPdAu[:,1],color='r',linestyle='-',label='Experiment')
	legend = ax.legend(loc='upper left', shadow=False,fontsize=12)
	legend.get_frame().set_facecolor('#FFFFFF')
	ax.set_ylabel('Lattice Parameter (A)')
	ax.set_xlabel('\%Au')
	ax.set_xlim(0, 1)
	ax.set_title('PdAu alloy lattice parameter')
	fig.set_size_inches(width, height)
	fig.savefig('PBESOL.pdf')
	plt.show()

########################################################################
########################################################################
if fracH:
	### DATA	
	pdx = [0,0.25,.5,.75,.875,1,1.125]
	pdy = [3.89,3.94,3.99,4.04,4.06,4.08,4.12]

	pdaux=[0,.25,.5,.75]
	pdauy=[3.99,4.01,4.05,4.10]
	
	### FIGURE
	fPd = interp1d(pdx, pdy, kind='cubic')
	xIntPd=np.linspace(0,1.125,30)
	fPdAu = interp1d(pdaux, pdauy, kind='cubic')
	xIntPdAu=np.linspace(0,.75,30)


	fig, (ax) = plt.subplots(nrows=1,ncols=1)

	ax.scatter(pdx,pdy,s=50,color='b',marker='.',linewidth='1')
	ax.scatter(pdaux,pdauy,s=50,color='r',marker='.',linewidth='1')
	ax.plot(xIntPd,fPd(xIntPd),color='b',linestyle='-',label='Pd (PBEsol)')
	ax.plot(xIntPdAu,fPdAu(xIntPdAu),color='r',linestyle='-',label='PdAu (PBEsol)')
	ax.plot((0,1.2),(4.1,4.1),color='b',ls='--',label='Pd, -0.9 V')
	ax.plot((0,1.2),(4.04,4.04),color='r',ls='--',label='PdAu, -0.9 V')
	ax.arrow(.44,4.04,0,-0.155,head_width=0.01, head_length=0.005, fc='k', ec='k')
	ax.arrow(1.08,4.1,0,-0.215,head_width=0.01, head_length=0.005, fc='k', ec='k')
	legend = ax.legend(bbox_to_anchor=(0.5, 0.5), shadow=False,fontsize=12)
	legend.get_frame().set_facecolor('#FFFFFF')
	ax.set_ylabel('Lattice Parameter (A)')
	ax.set_xlabel('H : Pd')
	ax.set_xlim(0, 1.2)
	ax.set_ylim(3.88, 4.15)
	ax.set_title('Determining quantity of intercalated H')
	fig.set_size_inches(width, height)
	fig.savefig('fracH.pdf')
	plt.show()
########################################################################
##############################################
if BEEFexpansion:
	### DATA
	pdx       = np.array([0,0.25,0.5, 0.75,0.875,1,  1.125])
	pdbeefy   = np.array([0,0.05,0.10,0.13,0.15, 0.17,0.22])
	pdpbesoly = np.array([0,0.06,0.10,0.15,0.17, 0.19,0.24])
	fPd       = interp1d(pdx, pdpbesoly, kind='cubic')
	xIntPd    = np.linspace(0,1.125,30)


	pdaubeefx   = np.array([0,   0.5,   1,  1.5])
	pdaubeefy   = np.array([0.00,0.05,0.10,0.16])
	pdaupbesolx = np.array([0, 0.25, 0.5,   1])
	pdaupbesoly = np.array([0,0.02,0.06,0.12])
	fPdAu       = interp1d(pdaupbesolx, pdaupbesoly, kind='cubic')
	xIntPdAu    = np.linspace(0,1,30)

	### FIGURE
	fig, (ax) = plt.subplots(nrows=1,ncols=1)

	ax.scatter(pdx,pdbeefy,s=50,color='b',marker='.',linewidth='1',label='Pd (BEEF-vdW)')
	ax.scatter(pdx,pdpbesoly,s=50,color='b',marker='x',linewidth='1',label='Pd (PBEsol)')
	ax.plot(xIntPd,fPd(xIntPd),color='b',linestyle='-')
	ax.scatter(pdaubeefx,pdaubeefy,s=50,color='r',marker='.',linewidth='1',label='PdAu (BEEF-vdW)')
	ax.scatter(pdaupbesolx,pdaupbesoly,s=50,color='r',marker='x',linewidth='1',label='PdAu (PBEsol)')
	ax.plot(xIntPdAu,fPdAu(xIntPdAu),color='r',linestyle='-')
	legend = ax.legend(loc='upper left', shadow=False,fontsize=12)
	legend.get_frame().set_facecolor('#FFFFFF')
	ax.set_ylabel('Lattice Expansion (A)')
	ax.set_xlabel('H : Pd')
	ax.set_xlim(0, 1.2)
	ax.set_ylim(0, 0.25)
	ax.set_title('Comparing predicted lattice expansions')
	fig.set_size_inches(width, height)
	fig.savefig('BEEFexpansion.pdf')
	plt.show()
