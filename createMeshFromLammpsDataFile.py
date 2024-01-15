#data reader
import bpy, bmesh
import numpy as np
from mathutils import Vector
import os
# How to use this script to create lammps input from blender scenes:
# if you export a unit cell or crystals, make sure all models contain only one type of atoms
# Information for the structure with parameters can be passed by modifying the object name. 
# Here, each section is separated by "_", starting with the name and the default charge of each atom.

# global variables
filename = "\\molecule.data"
datapath = os.path.realpath(r"PATHWAYTOFILE"+filename)
f = open(datapath,"r")
data = f.readlines()
f.close()
atoms = []
atomtypes = []
bonds=[]
appendAtoms=False
appendBonds= False
nextline = 0
# iterate over the lines in the file. First, the section with the atom positions is extracted
# Then the Bond data is read
for l in data:
    if("Atoms" in l):
        appendAtoms=True
        nextline = 1
    if("Bonds" in l):
        appendBonds=True
        nextline = 1
        appendAtoms=False
        print(appendBonds,nextline)
    #print(nextline)
    if(l=="\n" and appendAtoms):
        nextline +=1
    if(appendAtoms and nextline==2):
        atoms.append(l)
    if(l=="\n" and appendBonds):
        nextline +=1
    if(appendBonds and nextline==2):
        bonds.append(l)
# From the lines with atom information, the atom index, charge and type is saved to arrays in order from the file
atomPos = []
atomTypes = []
charges = [] # use to assign size of default sphere
atomInd=[]
for a in atoms[1:]:
    print(a)
    line = a.split()
    print(line)
    atomPos.append((float(line[4]),float(line[5]),float(line[6])))
    atomTypes.append(int(float(line[2])))
    charges.append(float(line[3]))
    atomInd.append(int(float(line[0]))-1)
# The Bond information is read
bondarray=[]
indsort = np.argsort(atomInd)
possort= [atomPos[i] for i in indsort]
for b in bonds[1:]:
    b = b.split()
    bondarray.append([int(float(b[2]))-1,int(float(b[3]))-1])
# First, a mesh for all atoms according to the bond data is generated and saved as molecule_mesh

molecule_mesh = bpy.data.meshes.new(filename.split(".")[0][1::]+'bonds')
molecule_mesh.from_pydata(possort, bondarray, [])
molecule_mesh.update()
molecule_object = bpy.data.objects.new(filename.split(".")[0][1::], molecule_mesh)
new_collection = bpy.data.collections.new(filename.split(".")[0][1::]+'_collection')
bpy.context.scene.collection.children.link(new_collection)
new_collection.objects.link(molecule_object)

# The atoms in the bond mesh are assigned to vertex groups named according the scheme chosen for this project
# First, the group index which corresponds with the atom type ID in the Lammps data file and the charge are used as ID_charge_
# Also, surfspheres are added and set as a child to newly created meshes consisting of the atoms of the group as single vertices
# by setting the instance to VERTS, all vertices are displayed as the sphere, which was also scaled with a self-brewn function
# to let the displayed atom size corellate with the charge

for at in np.unique(atomTypes):
    typeAt = [atomPos[i] for i in atomInd if atomTypes[i]==at]
    indAt = [atomInd[i] for i in atomInd if atomTypes[i]==at]
    print(len(typeAt), typeAt)
    print(indAt)
    charge_here = [charges[i] for i in atomInd if atomTypes[i]==at][0]
    print(charge_here, at)
    atomspheres = bpy.data.meshes.new(str(at))
    atomspheres.from_pydata(typeAt,[],[])
    atomspheres.update()
    newAts = bpy.data.objects.new(str(at),atomspheres)
    new_collection.objects.link(newAts)
    
    new_vertex_group = molecule_object.vertex_groups.new(name=f"{at}_{charge_here}_")
    vertex_group_data = indAt
    new_vertex_group.add(vertex_group_data, 1.0, 'ADD')
    
    bpy.ops.surface.primitive_nurbs_surface_sphere_add(radius=0.6+np.arctan(float(abs(charge_here))*3)*0.3,enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1,1,1))
    a = bpy.context.selected_objects[-1]
    print(bpy.context.selected_objects)
    a.name = str(at)
    a.parent=newAts
    newAts.instance_type='VERTS'
    
