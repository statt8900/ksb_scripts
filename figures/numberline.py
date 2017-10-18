import matplotlib.pyplot as plt
from math import cos,sin
#Emin,Emax = -.15,1.3
x = 10
def circ(ind):
	ang = 3. - ind/1.9
	return (3*cos(ang),3*sin(ang))
#def f(x): 
#	print x,Emin,(x - Emin)*(10)/(Emax-Emin)
#	return (x - Emin)*(10)

# set up the figure
fig = plt.figure()
ax 	= fig.add_subplot(111)
ax.set_xlim(-3,13)
ax.set_ylim(0,10)

# draw lines
xmin,xmax 	= -.2,1.3
y 			= 5
height 		= 1

plt.annotate(s='',xy=(x*xmin,y),xytext=(x*xmax,y),arrowprops=dict(arrowstyle='<->'),zorder=-1)
plt.vlines(0, y - height / 2., y + height / 2.)
plt.vlines(x*1.23, y - height / 2., y + height / 2.)

# draw a point on the line
px = [-.11,.02,.07,.08,.09,.16]
rx = [r'$CO_2/$'+r for r in [r'$CO$',r'$CH_3OH$',r'$C_2H_4$',r'$C_2H_5OH$',r'$C_3H_5OH$',r'$CH_4$']]
for i,p in enumerate(px): 
	plt.plot(x*p,y, 'ko', ms = 10, mfc = 'b')

	# add an arrow
	xpos,ypos = circ(i)
	plt.annotate(rx[i], (x*p,y), xytext = (x*p+xpos, y+ypos), 
	              arrowprops=dict(facecolor='black', shrink=0.2), 
	              horizontalalignment='right').draggable()

# add numbers
plt.plot(0,y, 'ko', ms = 10, mfc = 'r')
plt.plot(x*1.23,y, 'ko', ms = 10, mfc = 'g')

plt.text(x*0 , y-1, r'$H_2/H^+$', horizontalalignment='center')
plt.text(x*0 , y-1.8, r'0 V $vs$ RHE', horizontalalignment='center')
plt.text(x*1.23, y-1, r'$O_2/H_2O$', horizontalalignment='center')
plt.text(x*1.23, y-1.8, r'1.23 V $vs$ RHE', horizontalalignment='center')

plt.axis('off')
plt.show()