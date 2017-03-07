
###  Here I just do some basic initializations 

from scipy.stats import gaussian_kde
import pylab
import matplotlib.pyplot as plt
import h5py
import numpy as np
import requests


import requests
baseUrl = 'http://www.illustris-project.org/api/'
headers = {"api-key":"452a77edd4324ec835da440fc3fdc50b"}


def get(path, params=None):
     # make HTTP GET request to path
     r = requests.get(path, params=params, headers=headers)

     # raise exception if response code is not HTTP SUCCESS (200)
     r.raise_for_status()

     if r.headers['content-type'] == 'application/json':
         return r.json() # parse json responses automatically

     if 'content-disposition' in r.headers:
         filename = r.headers['content-disposition'].split("filename=")[1]
         with open(filename, 'wb') as f:
             f.write(r.content)
         return filename # return the filename string

     return r

#### Get full meta data

redshift = 2.2204460492503099e-16

scale_factor = 1.0 / (1+redshift) 
little_h = 0.704    

#### critical strings created for future use in lists
halo_dm_string = ([])
radiistring = ([])
densitystring = ([])
redshift_string = ([])
filename_list =([])
filename1_list = ([])
halo_id = ([])
parents = ([])
dmhalomasses = ([])
progenitors = ([])
dark_mass_percent = ([])



# first convert log solar masses into group catalog unit

mass_min = 10**12/ 1e10 * 0.704

mass_max = 10**12.3/ 1e10 * 0.704


#mass_max = 10**16/ 1e10 * 0.704
# form the search_query string by hand for once
search_query = "?mass__gt=" + str(mass_min) + "&mass__lt=" + str(mass_max)
print search_query



url = "http://www.illustris-project.org/api/Illustris-2-Dark/snapshots/135/subhalos/" + search_query

subhalos = get(url)
print subhalos['count']
#### OVERRIDE DEFAULT LENGTH OF 100
subhalos1 = get((url),{'limit': subhalos['count'], 'order_by':'-mass'}) 
#subhalos1 = get((url),{'limit': , 'order_by':'-mass_dm'})


c = [ subhalos1['results'][i]['id'] for i in range(subhalos['count']) ]
#c = [ subhalos1['results'][i]['id'] for i in range(5) ]



def running_sum(a):
  tot = 0
  for item in a:
    tot += item
    yield tot


params = {'dm':'Coordinates'}
url_parents = ([])

#c = [1]
print len(c)



n = 0

while n < len(c):
	import h5py
	import numpy as np

	id = c[n]

	
	print id
	url = "http://www.illustris-project.org/api/Illustris-2-Dark/snapshots/135/subhalos/" + str(id)
	print url
	
	sub = get(url) # get json response of subhalo properties
	saved_filename = get(url + "/cutout.hdf5",params) # get and save HDF5 cutout file
	dmhalomass = sub['mass']

	progenitor = sub['prog_sfid']	
	parent = sub['parent']	

	parents.insert(n,parent)
	progenitors.insert(n,progenitor)
	dmhalomasses.insert(n,dmhalomass)
	redshift_string.insert(n,redshift)	
	
	halo_id.insert(n,id)
		
	
	Omega_lambda = 0.7
	Omega_M0 = 0.3
	Ho = 2.2683 * 10 ** -18
	Ho_squared = (Ho)**2 
	redde = (1 + redshift)**3
	
	H_square = Omega_lambda + redde*(Omega_M0)
	H_squared = H_square * Ho
	G = 6.67408 *10 ** -11
	pi = 3.14159265359
	P_critical = 3 * H_squared / (8 * pi * G)
	P_200 = P_critical * 200
	print P_200
	# prepare dict to hold result arrays
	fields = ['snap','id','mass']
	radialstarmass_density =([])
	radialstarmass = ([])
	radial_distance = ([])
	dm_density = ([])	
	
	



	
	with h5py.File(saved_filename) as f:
	
		dx = f['PartType1']['Coordinates'][:,0] - sub['pos_x']
		dy = f['PartType1']['Coordinates'][:,1] - sub['pos_y']
		dz = f['PartType1']['Coordinates'][:,2] - sub['pos_z']
		
		rr2 = np.sqrt(dx**2 + dy**2 + dz**2)
		rr2 *= scale_factor/little_h # ckpc/h -> physical kpc
		print len(rr2)
		rrr2 = sorted(rr2)
		darkmattermass_density = ([])
		radial_distance2 = ([])
		coun = 0
		iter = len(rrr2) / (3)
		remainder = len(rrr2) % (3)
		dark_matters = ([])
		outer_radius2 = ([])
		while coun < iter:
			top5 = rrr2[:3]
			outer_radius2.insert(coun,top5[2])
			totaldm_mass = len(top5) * 0.0035271
			dark_matters.insert(coun,totaldm_mass)
			del rrr2[:3]
			coun = coun + 1	
		outer_radius2.insert(0,0)
		
		
		m = 0
		dark_mass_density = ([])
		while m < iter :
			if remainder > 0 :
				volume2 = (4.0/3.0)*3.14*((outer_radius2[m + 1])**3) - (4.0/3.0)*3.14*((outer_radius2[m])**3)
			elif remainder == 0:
				volume2 = (4.0/3.0)*3.14*((outer_radius2[m + 1])**3) - (4.0/3.0)*3.14*((outer_radius2[m])**3)	
			density = (dark_matters[m]) / (volume2)
			dark_mass_density.insert(m,density)
			m = m + 1
		
		outer_radius2.pop(0)
		
		
		print density
	for x in dark_mass_density:
		if x > P_200:
			dm_density.append(x)
			
	print "printing Hsquared"
	
	total_dm_mass = list(running_sum(dark_matters))
	del total_dm_mass[0]
	y1 =  dm_density 
	g = len(y1)
	r1 = outer_radius2[:g]
	unofficial = len(y1) * 13 *  0.0035271
	snapnum = 130
	
	filename = "darkmatterdensityhalo" + str(snapnum) + "redshift" + str(id) + ".txt"
	filename1 = "darkmatterradii"  + str(snapnum) + "redshift" + str(id) + ".txt"
	filename6 = "parenthalos.txt"
	filename7 =  "darkmattermass.txt"
	filename8 = "progenitor.txt"
	filename9 = "descendant.txt"
	filename10 = "halo_id.txt"
	filename11 = "redshift.txt"
	filename12 = "massatradiusofhalo"  + str(snapnum) + "redshift" + str(id) + ".txt"
	filename13 = "mass_at_radii.txt"
	filename16 = "dark_mass_percent.txt"
	filename15 = "totaldarkmatter.txt"
	def format(value):
		return "%.18f" % value
	#del y1[0]
	y2 = [format(x) for x in y1]
	#del r1[0]
	
	def formats(values):
		return "%.8f" % values	
	r2 = [formats(x) for x in r1]
	
	print len(y2),  len(r2)
	import csv
	resultFile = open(filename,'wb')
	resultFile1 = open(filename1,'wb')
	resultFile2 = open(filename12,'wb')
	wr = csv.writer(resultFile, dialect='excel')
	print "now printing ************************************* writerows"
	wr.writerows(y2)
	resultFile.close()
	wr1 = csv.writer(resultFile1, dialect='excel')
	wr1.writerows(r2)
	resultFile1.close()
	radiistring.insert(n,filename1)
	densitystring.insert(n,filename)
	halo_dm_string.insert(n,filename12)
	total_dm_masses= [format(x) for x in total_dm_mass]
	wr2 = csv.writer(resultFile2, dialect='excel')
	wr2.writerows(total_dm_masses)
	resultFile2.close()
	print len(outer_radius2)
	n = n + 1


print "printing redshift"
print redshift_string
# Finally, after this iterative process has produced all the necessary textfiles,
# produce a final text file containing the list of all the text files produced so that
# the non linear regression program can receive this list and know which text files to
# call when fitting the data.
print "printing len of dmdensity"
print len(dm_density)
#print filename_list
#print filename1_list
parent = [int(x) for x in str(parent)]

dmhalomass = [float(dmhalomass)]
print parent,dmhalomass
filename2 = "halos.txt"
filename3 = "radii.txt"
filename4 = "listofradii.txt"
filename5 = "listofhalos.txt"

#######   Here I write filename_list to it's own text file which I will then import into regression so that I know which files to do regression on

data = open(filename2, "w")

for c in filename_list:
	data.write(c)

data.close()
	

data = open(filename3, "w")

for c in filename1_list:
	data.write(c)

data.close()


data = open(filename4,"w")

for c in radiistring:
	data.write("%s\n" % c)
	
data.close()


data = open(filename5,"w")

for c in densitystring:
	data.write("%s\n" % c)
	
data.close()



data = open(filename6,"w")

for c in parents:
	data.write("%s\n" % c)
	
data.close()



data = open(filename7,"w")

for c in dmhalomasses:
	data.write("%s\n" % c)
	
data.close()




data = open(filename8,"w")

for c in progenitors:
	data.write("%s\n" % c)
	
data.close()

print "printing dmhalomass, and progenitor"
print progenitors
print dmhalomasses
print parents

data = open(filename9,"w")

#for c in descendant:
#	data.write("%s\n" % c)
	
#data.close()




data = open(filename10,"w")

for c in halo_id :
	data.write("%s\n" % c)
	
data.close()


	

data = open(filename11,"w")

for c in redshift_string :
	data.write("%s\n" % c)
	
data.close()



data = open(filename13,"w")

for c in halo_dm_string :
	data.write("%s\n" % c)
	
data.close()




data = open(filename16,"w")

for c in dark_mass_percent :
	data.write("%s\n" % c)
	
data.close()






	
	
 
		    
