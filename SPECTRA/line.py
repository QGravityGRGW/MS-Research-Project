import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.io import ascii
import air_vac
def line(path=None,z=None):
	filename=sorted(glob.glob(path))
	list_1=[770.409,780.324,937.814,949.742,977.030,989.790,991.514,991.579,1025.722,1031.912,1037.613,1066.660,1215.670,1238.821,1242.804,1260.422,1264.730,1302.168,1334.532,1335.708,1393.755]
	list_2=[2320.951,2323.500,2324.690,2648.710,2733.289,2782.700,2795.528,2802.705,2829.360,2835.740,2853.670,2868.210,2928.000,2945.106,3132.794,3187.745,3203.100,3312.329,3345.821,3425.881,3444.052,3466.497,3466.543,3487.727,3586.320,3662.500,3686.831,3691.551,3697.157,3703.859,3711.977,3721.945,3726.032,3728.815,3734.369,3750.158,3758.920,3770.637,3797.904,3835.391,3839.270,3868.760,3888.647,3889.064,3891.280,3911.330,3967.470,3970.079,4026.190,4068.600,4071.240,4076.349,4101.742,4143.761,4178.862]
	lmg2=2799.117
	lmg21=2796.352
	lmg22=2803.531
	lo21=3726.032
	lo22=3729.875
	lha=6564.61
	lhb=4862.68
	lo31=4932.603 	
	lo32=4960.295
	lo33=5008.240
	lk=3934.777
	lh=3969.588
	lg=4305.61
	lmg=5176.7
	lna=5895.6
	lca21=8500.36
	lca22=8544.44
	lca23=8664.52
	wave=[]
	flux=[]
	for file in filename:
		a1=ascii.read(file)
		w,f=air_vac.convert(file)
		wave.append(w)
		flux.append(f)
	for i in range(len(list_2)):
		plt.axvline(list_2[i]*(1+z),color='C'+str(i))
		plt.text(list_2[i]*(1+z),-70,str(list_2[i]))
	for i in range(len(wave)):
		plt.step(wave[i],flux[i],lw=0.6,color='black')
		plt.axvline(lmg2*(1+z),color='red')
		plt.text(lmg2*(1+z),-40,'MgII_e_2799')
		#plt.axvline(lmg21*(1+z),color='aqua')
		#plt.text(lmg21*(1+z),-40,'MgII_a1')
		#plt.axvline(lmg22*(1+z),color='aqua')
		#plt.text(lmg22*(1+z),-40,'MgII_a2')
		plt.axvline(lmg*(1+z),color='yellow')
		plt.text(lmg*(1+z),-40,'Mg_a')
		plt.axvline(lo21*(1+z),color='cyan')
		plt.text(lo21*(1+z),-40,'OII_e_3727')
		plt.axvline(lo22*(1+z),color='blue')
		plt.text(lo22*(1+z),-30,'OII_e_3729')
		plt.axvline(lha*(1+z),color='pink')
		plt.text(lha*(1+z),-40,'HIa_e')
		plt.axvline(lhb*(1+z),color='green')
		plt.text(lhb*(1+z),-40,'HIb_e')
		plt.axvline(lo31*(1+z),color='orange')
		plt.text(lo31*(1+z),-40,'OIII_e1')
		plt.axvline(lo32*(1+z),color='magenta')
		plt.text(lo32*(1+z),-40,'OIII_e2')
		plt.axvline(lo33*(1+z),color='maroon')
		plt.text(lo33*(1+z),-40,'OIII_e3')
		plt.axvline(lh*(1+z),color='aqua')
		plt.text(lh*(1+z),-40,'H_a')
		plt.axvline(lk*(1+z),color='indigo')
		plt.text(lk*(1+z),-40,'K_a')
		plt.axvline(lg*(1+z),color='violet')
		plt.text(lg*(1+z),-40,'G_a')
		plt.show()
	return None
