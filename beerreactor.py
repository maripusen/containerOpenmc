import openmc 
import matplotlib.pyplot as plt
import numpy as np
import random



#Defining uranium oxide material, set to real life density

urox = openmc.Material(1, "uo2")
urox.add_nuclide("U235", 0.02)
urox.add_nuclide("U238", 0.98)
urox.add_nuclide("O16", 2.0)
urox.set_density("g/cm3", 10.97)


#defining water-moderator material
water = openmc.Material(2, "h2o")
water.add_nuclide("H2", 2.0)
water.add_nuclide("O16",1)
water.add_s_alpha_beta('c_D_in_D2O')
water.set_density("g/cm3", 1.1)

#Defining a zircaloy 
zirc = openmc.Material(3, "Zr")
zirc.add_element("Zr", 1.0)
zirc.set_density("g/cm3", 6.49)

#defining helium gas
helium = openmc.Material(4,"He")
helium.add_element("He",1.0)
helium.set_density("g/cm3",0.000178)

#defining a graphite moderator
graphite = openmc.Material(5, "Graph")
graphite.add_element("C", 1.0)
graphite.set_density("g/cm3", 2.09)

weird = openmc.Material(6, "Weird")
weird.add_element("C", 1.0)
weird.set_density("g/cm3", 1.0)



#Defining the major components of beer
alcohol = openmc.Material(7,"ethanol")
alcohol.add_element("C", 2)
alcohol.add_nuclide("O16", 1)
alcohol.add_element("H", 6)
alcohol.set_density("g/cm3", 0.789)

water_b = openmc.Material(8,"water")
water_b.add_element("H",2)
water_b.add_nuclide("O16",1)
water_b.set_density("g/cm3", 1)


#Carbonation materials
carbonation_solution = openmc.Material(9, "CO2")
carbonation_solution.add_element("C",1)
carbonation_solution.add_nuclide("O16",2)
carbonation_solution.set_density("g/cm3",0.0057)

carbonation_gas = openmc.Material(10, "CO2")
carbonation_gas.add_element("C",1)
carbonation_gas.add_nuclide("O16",2)
carbonation_gas.set_density("kg/m3", 0.001965)


#Humulone
humulone = openmc.Material(11, "humulone")
humulone.add_element("C",21)
humulone.add_element("H",30)
humulone.add_nuclide("O16",5)
humulone.set_density("g/cm3",1.157)

#Sugars
glucose = openmc.Material(12,"sugar")
glucose.add_element("C", 6)
glucose.add_nuclide("O16", 6)
glucose.add_element("H",12)
glucose.set_density("g/cm3",1.56)

#Minerals
#https://link.springer.com/article/10.1007/s00217-018-3098-0
calcium = openmc.Material(13, "calcium")
calcium.add_element("Ca",1)
calcium.set_density("g/cm3", 1.55)

iron = openmc.Material(14, "iron")
iron.add_element("Fe",1)
iron.set_density("g/cm3",7.87)

copper = openmc.Material(15, "copper")
copper.add_element("Cu",1)
copper.set_density("g/cm3",8.98)

zinc = openmc.Material(16, "zinc")
zinc.add_element("Zn",1)
zinc.set_density("g/cm3",7.13)

chlorine = openmc.Material(17, "chlorine")
chlorine.add_element("Cl",1)
chlorine.set_density("g/cm3",0.003214)


potassium = openmc.Material(18, "potassium")
potassium.add_element("K",1)
potassium.set_density("g/cm3", 0.856)

magnesium = openmc.Material(19, "magnesium")
magnesium.add_element("Mg",1)
magnesium.set_density("g/cm3", 1.738)

manganese = openmc.Material(20,"manganese")
manganese.add_element("Mn",1)
manganese.set_density("g/cm3",7.21)

sodium = openmc.Material(21, "sodium")
sodium.add_element("Na",1)
sodium.set_density("g/cm3",0.97)

silicon = openmc.Material(22,"silicon")
silicon.add_element("Si",1)
silicon.set_density("g/cm3", 2.39)

phosphorus = openmc.Material(23,"Phosphorus")
phosphorus.add_element("P",1)
phosphorus.set_density("g/cm3",1.82)

mins = [calcium, iron, copper, zinc, chlorine, potassium, magnesium, manganese, sodium, silicon, phosphorus]
mats = [alcohol, water_b, carbonation_solution,humulone, glucose, calcium, chlorine, potassium, phosphorus]
densities = [0.15, 0.80-0.00068, 0.0025, 0.0001, 0.0474,0.0003, 0.00003, 0.0001, 0.00025 ]


beer = openmc.Material.mix_materials(mats, densities, "vo")

def randxy():
    x = random.uniform(-0.7,0.7)
    y = random.uniform(-0.7,0.7)
    return x,y

def bubbles(num, radius):
    regs = []
    bls = []
    for i in range (0,num+1):
        rad = random.uniform(2.5e-4, 1.0e-4)
        k = 0
        while k < 1:
            x,y = randxy()
            if x**2 + y**2 > radius +0.00003:
                k +=1
        
        z = random.uniform(1,-1)
        reg = openmc.Sphere(x0 = x, y0 = y, z0 = z,r = rad)
        reg = -reg 
        nucell = openmc.Cell(cell_id = i + 6,name = "cell %.i" % i)
        nucell.fill = carbonation_gas
        nucell.region = reg
        regs.append(reg)
        bls.append(nucell)
    return [bls,regs]





#Exporting a XML file wih the materials. 
materials = openmc.Materials([urox, water, zirc, helium, graphite, 
                              carbonation_gas,alcohol, water_b, carbonation_solution,humulone, glucose, 
                              beer, calcium, iron, copper, zinc, chlorine, potassium, magnesium, manganese, 
                              sodium, silicon, phosphorus])
materials.cross_sections = "/home/mtorsvoll.linux/nndc_hdf5/cross_sections.xml"
materials.export_to_xml()


"""""Defining geometry"""
#Defining cylinder regions for the pincell
inneredge = "transmission"

base = 0.4 - 0.07

urcyl = openmc.ZCylinder(r = base, boundary_type = inneredge)
hecyl2 = openmc.ZCylinder(r = base + 0.02, boundary_type = inneredge)
zirccyl2 = openmc.ZCylinder(r = base + 0.07, boundary_type = inneredge)

#"""

#making box of 1.4 cm x 1.4cm x 2cm

#defining a simple variable for box edge properties
edge = "reflective"

#"""
xmin = openmc.XPlane(-0.7 , boundary_type = edge)
xmax = openmc.XPlane(0.7, boundary_type = edge)
ymin = openmc.YPlane(-0.7, boundary_type = edge)
ymax = openmc.YPlane(0.7, boundary_type = edge) 
#"""
zmin = openmc.ZPlane(-1, boundary_type = edge)
zmax = openmc.ZPlane(1, boundary_type = edge)

box = +xmin & - xmax & +ymin & - ymax & +zmin & - zmax


"""Making fuel cells and other material regions"""
fuel_reg = - urcyl & -zmax & +zmin
he_reg = + urcyl & - hecyl2 & -zmax & +zmin
zircclad = + hecyl2 & - zirccyl2 & -zmax & +zmin
water_reg = +zirccyl2 & box



fuel = openmc.Cell(name = "fuel")
fuel.fill = urox
fuel.region = fuel_reg

hegap = openmc.Cell(name = "hegap")
hegap.fill = helium
hegap.region = he_reg

cladding = openmc.Cell(name = "cladding")
cladding.fill = zirc
cladding.region = zircclad 



#Defining 

bubble = bubbles(int((1.4*1.4*2)/500*10**6),0.42)
#bubble = bubbles(10,0.42)

for i in range (0,len(bubble[1])):
    water_reg = ~ bubble[1][i] & water_reg

moderator = openmc.Cell(name = "moderator")
moderator.fill = beer
moderator.region = water_reg

cells = [fuel, hegap, cladding, moderator]

cells = cells + bubble[0]
root_universe = openmc.Universe(cells = cells)



#Export geometry
geometry = openmc.Geometry()
geometry.root_universe = root_universe
geometry.export_to_xml()


#"""
if __name__ == "__main__":
    root_universe.plot(basis = "xy", 
              colors = {fuel: "green", hegap: "brown", 
                        cladding: "gray", moderator: "yellow"})
    plt.show()

    root_universe.plot(basis = "xz", 
              colors = {fuel: "green", hegap: "brown", 
                        cladding: "gray", moderator: "yellow"})
    plt.show()

#"""

# OpenMC simulation parameters
settings = openmc.Settings()
settings.batches = 100
settings.inactive = 10
settings.particles = 10000


# Create an initial uniform spatial source distribution overfissionable zones
bounds = [-0.7, -0.7, -1, 0.7, 0.7, 1]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],only_fissionable=True)
settings.source = openmc.IndependentSource(space=uniform_dist)

# Export to "settings.xml"
settings.export_to_xml()



#extracting mesh coords
meshcor = box.bounding_box

#Constructing a mesh 
mesh = openmc.RegularMesh()
mesh.dimension = [100,100]
mesh.lower_left= [-0.7,-0.7]
mesh.upper_right= [0.7,0.7]



# Create mesh filter for tally
mesh_filter = openmc.MeshFilter(mesh)
#esh_filter.mesh = mesh

# Instantiate an empty Tallies object
tallies_file = openmc.Tallies()


erange = np.logspace(np.log10(10e-5), np.log10(20e6), 501)
energy_filter = openmc.EnergyFilter(erange)

# Create mesh tally to score flux and fission rate
tally = openmc.Tally(name = "Tally")
tally.filters = [mesh_filter]
tally.scores = ["flux", 'prompt-nu-fission']
tallies_file.append(tally)


tally2 = openmc.Tally(name = "Energy")
tally2.filters = [energy_filter]
tally2.scores = ["flux"]
tallies_file.append(tally2)




# Resonance Escape Probability tallies
therm_abs_rate = openmc.Tally(name='therm. abs. rate')
therm_abs_rate.scores = ['absorption']
therm_abs_rate.filters = [openmc.EnergyFilter([0.0,0.625])]
tallies_file.append(therm_abs_rate)

# Thermal Flux Utilization tallies
fuel_therm_abs_rate = openmc.Tally(name='fuel therm. abs. rate')
fuel_therm_abs_rate.scores = ['absorption']
fuel_therm_abs_rate.filters = [openmc.EnergyFilter([0.0, 0.625]),
                               openmc.CellFilter([fuel])]
tallies_file.append(fuel_therm_abs_rate)


tallies_file.export_to_xml()


if __name__ == "__main__":
    #Running the simulation
    openmc.run()

