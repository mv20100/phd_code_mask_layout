from lxml import etree
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

pp = PdfPages('statistics.pdf')
tree = etree.parse("output.xml")

#Plot surface vs. volume curve
extendedCavs = tree.findall(".//SpecialExtendedCavs2")
fig = plt.figure(num=None, figsize=(20, 20))
ax = fig.add_subplot(111)
ax.set_title('Volume vs. Area')
ax.set_xlabel('Inner cavity area')
ax.set_ylabel('Inner cavity perimeter')
for extendedCav in extendedCavs:
	totalArea = int(float(extendedCav.get("totalArea")))
	totalPerimeter = int(float(extendedCav.get("totalPerimeter")))
	label = extendedCav.getparent().tag
	ax.plot(totalArea, totalPerimeter, 'ro')
	ax.annotate(
		label, 
		xy = (totalArea, totalPerimeter), xytext = (40, 20),
		textcoords = 'offset points', ha = 'right', va = 'bottom',
		bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
		arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
plt.draw()
pp.savefig()

#Plot heater and sensor resistance
heaters = tree.findall(".//Heater")
labels = []
heatResistances = []
sensResistances = []
for heater in heaters:
	heatResistances.append(int(float(heater.getchildren()[0].get("resistance"))))
	sensResistances.append(int(float(heater.getparent().find("Sensor").getchildren()[0].get("resistance"))))
	labels.append(heater.getparent().tag)

fig = plt.figure(num=None, figsize=(20, 20))
ax1 = fig.add_subplot(211)
ax1.bar(range(len(heatResistances)), heatResistances, align='center')
ax1.set_ylabel('Heater resistance')
ax1.get_xaxis().set_visible(False)
ax2 = fig.add_subplot(212)
ax2.bar(range(len(sensResistances)), sensResistances, align='center')
ax2.set_ylabel('Sensor resistance')
xrange = np.array(range(len(labels)))
plt.xticks(xrange.tolist(), labels, size='small',rotation=45,ha="right")
# plt.gcf().subplots_adjust(bottom=0.15)
plt.tight_layout()
plt.draw()
pp.savefig()
pp.close()
# plt.show()