# Script to export data in a scene in blender to a data file for MD simulations in the LAMMPS format
import bpy, bmesh
import numpy as np
from mathutils import Vector
import os
# How to use this script to create lammps input from blender scenes:
# 1. Select all atoms and molecules which you want to export
#   all objects follow an object naming convention which decides the procedure for data export. 
#   For all types of atoms, information of the ID type in LAMMPS and a charge value is required. All parts are separated by _. 
#   This may vary with respect to molecules, crystals as explained respectively. Always include an additional _ after the name!
#   if you export a unit cell of a metal or ionic crystals, make sure all models contain only one type of atoms

# 2. Be careful to apply all transforms and modifiers before exporting. The global coordinates of all vertices are used directly representing the space in the MD units
# global variables
xlo=0
xhi= 0
ylo=0
yhi=0
zlo=0
zhi = 0
datapath = os.path.realpath(r"DATAPATH")
datafilename = "PEO3chains" # just teh base name without filetype suffix

# will return the default box dimensions as a bounding box of all atoms
def CheckBox(x,y,z):
    global xlo, xhi, ylo, yhi, zlo, zhi
    if(x<xlo):
        xlo=x
    if(x>xhi):
        xhi=x
    if(y<ylo):
        ylo=y
    if(y>yhi):
        yhi=y
    if(z<zlo):
        zlo=z
    if(z>zhi):
        zhi=z
    return xlo, xhi, ylo, yhi, zlo, zhi
  
# default function for writing full style atom data in lammps data format
# all verticies for the mesh are treated as atoms from one atom type
def MakeAtoms(inputPointMesh, index,out_moltypes):
    output= []
    
    global xlo, xhi, ylo, yhi, zlo, zhi
    mesh = inputPointMesh.data
    nam = inputPointMesh.name.split("_")
    print(nam)
    typeID = nam[2][4:]
    # get the bmesh data
    if mesh.is_editmode:
        bm = bmesh.from_edit_mesh(mesh)
    else:
        bm = bmesh.new()
        bm.from_mesh(mesh)

    # cycle through all vertices
    for vert in bm.verts:
        #+print(vert.co)
        #print(nam)
        xlo, xhi, ylo, yhi, zlo, zhi=CheckBox(vert.co[0],vert.co[1],vert.co[2])
        output.append(f"{index} 1 {typeID} {nam[1]} {vert.co[0]} {vert.co[1]} {vert.co[2]}\n")
        index +=1
        #print(vert.co[0],vert.co[2])
        
    out_moltypes[3].append(typeID)
    return index, output,out_moltypes

# polarizable two point atoms, which may overlap - require bond data for each pair
# They also require a CS-Info at the end, which is the Atom index and the core-shell ID
# The core shell ID is the molecule ID
# in each selection, there should be only one type of atoms or bonds
# typeID and bond type are extracted from the name
def MakePolarizable(inputPointMesh, index, moleculeID, lastbond,out_moltypes):
    atomOutput= []
    bond_out = []
    cs_out=[]
    mesh = inputPointMesh.data
    nam = inputPointMesh.name.split("_")
    #print(nam)
    
    typeID = int(nam[4][4:])
    bondType = int(nam[5][4:])
    # get the bmesh data
    if mesh.is_editmode:
        bm = bmesh.from_edit_mesh(mesh)
    else:
        bm = bmesh.new()
        bm.from_mesh(mesh)
    global xlo, xhi, ylo, yhi, zlo, zhi
    # cycle through all vertices
    for vert in bm.verts:
        lastbond+=1
        xlo, xhi, ylo, yhi, zlo, zhi=CheckBox(vert.co[0],vert.co[1],vert.co[2])
        #print(index, type(index), f"{1+1} a")
        atomOutput.append(f"{index} {moleculeID} {typeID} {nam[2]} {vert.co[0]} {vert.co[1]} {vert.co[2]}\n")
        atomOutput.append(f"{index+1} {moleculeID} {int(typeID)+1} {nam[3]} {vert.co[0]} {vert.co[1]} {vert.co[2]}\n")
        
        bond_out.append(f"{lastbond} {bondType} {index} {index+2} \n")
        cs_out.append(f"{index} {moleculeID}\n")
        cs_out.append(f"{index+1} {moleculeID}\n")
        index +=2
        moleculeID+=1
        #print(vert.co[0],vert.co[2])
    out_moltypes[0].append(bondType)
    out_moltypes[3].append(typeID)
    out_moltypes[3].append(typeID+1)
    return index, atomOutput, bond_out, cs_out, moleculeID,out_moltypes

# meshes with edge data use vertex groups to assign the vertices to atom groups
# Make sure to not have unassigned vertices!
# later, a default helper file including information about the created bond, angle and dihedral types is created, with which the force field parameters can be assigned
# This returns all the topological types, depending on the particle type determined by the vertex group
# the atom IDs and charges are stored in the name of the vertex group to which each vertex in the mesh representing a molecule needs to be assigned to.
# A molecule must only contain one connected mesh with all its covalent bonds
# Mesh name Name_Mol_, Vertex group names: type_charge_helpingInfo

# meshes including edge data in which atoms are marked by a material are used to create a molecule representation
# bonds, angles and dihedrals are assigned to default groups 
# To identify the nature of multiple type, an extra file is written in which an example is given for each bond/angle or dihedral
# this also serves as a template to write a PARM.lammps file
# The bond angle or dihedral types which if no extra force constant is added will then need to be written as zero.
# the atom IDs and charges are stored in the name of a vertex group to which each vertex in the mesh representing a molecule needs to be assigned to.
# A molecule must only contain one connected mesh with all its covalent bonds
# Mesh name: Name_Mol, Vertex group names: type_charge
def MakeMolecule(inputPointMesh, index, moleculeID, lastbond, lastangle, lastdihedral, out_moltypes):
    typearray = []
    atomtype_helperarray=[]
    out_atoms= []
    out_bonds=[]
    out_angles=[]
    out_dihedrals=[]
    chargearray=[]
    global xlo, xhi, ylo, yhi, zlo, zhi
    mesh = inputPointMesh.data
    nam = inputPointMesh.name.split("_")
    print(nam)
    typeID = nam[2][4:]
    # get the bmesh data
    if mesh.is_editmode:
        bm = bmesh.from_edit_mesh(mesh)
    else:
        bm = bmesh.new()
        bm.from_mesh(mesh)
    for i,group in enumerate(inputPointMesh.vertex_groups):
        nam=group.name.split("_")
        print(nam)
        atomtype_helperarray+= [ v.index for v in inputPointMesh.data.vertices if i in [ vg.group for vg in v.groups ] ]
        typearray += [ int(float(nam[0])) for v in inputPointMesh.data.vertices if i in [ vg.group for vg in v.groups ] ]
        chargearray += [ nam[1] for v in inputPointMesh.data.vertices if i in [ vg.group for vg in v.groups ] ]
        
    print([v.index for v in bm.verts])
    print("types: ", typearray)
    print("verts ", atomtype_helperarray)
    atomtypearray = [typearray[v] for v in np.argsort(atomtype_helperarray)]
    chargearray = [chargearray[v] for v in np.argsort(atomtype_helperarray)]
    bondtypearray=[]
    angletypearray=[]
    dihedraltypearray=[]
    for i,e in enumerate(inputPointMesh.data.edges):
        a = e.vertices[0]
        b = e.vertices[1]
        n = np.sort([a,b])
        out_bonds.append(f"{a} {b}")
        bondtypearray.append(f"{atomtypearray[n[0]]} {atomtypearray[n[1]]}")
        for j,ang in enumerate(inputPointMesh.data.edges):
            c = ang.vertices[0]
            d = ang.vertices[1]
            subchain=None
            if(i<j):
                if(b==c):
                    dtinter=f"{atomtypearray[a]} {atomtypearray[b]} {atomtypearray[d]}"
                    if(dtinter in out_moltypes[2] or dtinter in angletypearray):
                        subchain=[a,b,d]
                    elif(dtinter[::-1] in out_moltypes[2] or dtinter in angletypearray):
                        subchain=[d,b,a]
                    else: 
                        subchain=[a,b,d]
                    out_angles.append(f"{subchain[0]+1} {subchain[1]+1} {subchain[2]+1}")
                    angletypearray.append(dtinter)
                elif(b==d):
                    #out_angles.append(f"{a} {b} {c}")
                    dtinter=f"{atomtypearray[a]} {atomtypearray[b]} {atomtypearray[c]}"
                    if(dtinter in out_moltypes[2] or dtinter in angletypearray):
                        subchain=[a,b,c]
                    elif(dtinter[::-1] in out_moltypes[2] or dtinter in angletypearray):
                        subchain=[d,b,a]
                    else: 
                        subchain=[a,b,c]
                    angletypearray.append(dtinter)
                    out_angles.append(f"{subchain[0]+1} {subchain[1]+1} {subchain[2]+1}")
                elif(a==c):
                    #out_angles.append(f"{b} {c} {d}")
                    dtinter=f"{atomtypearray[b]} {atomtypearray[c]} {atomtypearray[d]}"
                    if(dtinter in out_moltypes[2] or dtinter in angletypearray):
                        subchain=[b,c,d]
                    elif(dtinter[::-1] in out_moltypes[2] or dtinter in angletypearray):
                        subchain=[d,c,b]
                    else: 
                        subchain=[b,c,d]
                    angletypearray.append(dtinter)
                    out_angles.append(f"{subchain[0]+1} {subchain[1]+1} {subchain[2]+1}")
                elif(a==d):
                    #out_angles.append(f"{b} {a} {c}")
                    dtinter=f"{atomtypearray[b]} {atomtypearray[a]} {atomtypearray[c]}"
                    if(dtinter in out_moltypes[2] or dtinter in angletypearray):
                        subchain=[b,a,c]
                    elif(dtinter[::-1] in out_moltypes[2] or dtinter in angletypearray):
                        subchain=[c,a,b]
                    else: 
                        subchain=[b,a,c]
                    angletypearray.append(dtinter)
                    out_angles.append(f"{subchain[0]+1} {subchain[1]+1} {subchain[2]+1}")
                for k,e_dihed in enumerate(inputPointMesh.data.edges):
                    e = e_dihed.vertices[0]
                    f = e_dihed.vertices[1]
                    if(j<k and subchain != None):
                        if(e==subchain[0]):
                            #print("matchdihed: ", i,j,k,": ", a, b,c,d,e,f)
                            out_dihedrals.append(f"{f+1} {subchain[0]+1} {subchain[1]+1} {subchain[2]+1}")
                            
                            dtinter = f"{atomtypearray[f]} {atomtypearray[subchain[0]]} {atomtypearray[subchain[1]]} {atomtypearray[subchain[2]]}"
                            if(dtinter in out_moltypes[3] or dtinter in dihedraltypearray):
                                dihedraltypearray.append(dtinter)
                            else: 
                                dihedraltypearray.append(dtinter[::-1])
                        elif(e==subchain[2]):
                            #print("matchdihed: ", i,j,k,": ", a, b,c,d,e,f)
                            out_dihedrals.append(f"{f+1} {subchain[0]+1} {subchain[1]+1} {subchain[2]+1}")
                            dtinter=f"{atomtypearray[f]} {atomtypearray[subchain[2]]} {atomtypearray[subchain[1]]} {atomtypearray[subchain[0]]}"
                            if(dtinter in out_moltypes[3] or dtinter in dihedraltypearray):
                                dihedraltypearray.append(dtinter)
                            else: 
                                dihedraltypearray.append(dtinter[::-1])

                        elif(f==subchain[0]):
                            #print("matchdihed: ", i,j,k,": ", a, b,c,d,e,f)
                            out_dihedrals.append(f"{e+1} {subchain[0]+1} {subchain[1]+1} {subchain[2]+1}")
                            dtinter=f"{atomtypearray[e]} {atomtypearray[subchain[0]]} {atomtypearray[subchain[1]]} {atomtypearray[subchain[2]]}"
                            if(dtinter in out_moltypes[3] or dtinter in dihedraltypearray):
                                dihedraltypearray.append(dtinter)
                            else: 
                                dihedraltypearray.append(dtinter[::-1])
                        elif(f==subchain[2]):
                            #print("matchdihed: ", i,j,k,": ", a, b,c,d,e,f)
                            out_dihedrals.append(f"{e+1} {subchain[0]+1} {subchain[1]+1} {subchain[2]+1}")
                            dtinter=f"{atomtypearray[e]} {atomtypearray[subchain[2]]} {atomtypearray[subchain[1]]} {atomtypearray[subchain[0]]}"
                            if(dtinter in out_moltypes[3] or dtinter in dihedraltypearray):
                                dihedraltypearray.append(dtinter)
                            else: 
                                dihedraltypearray.append(dtinter[::-1])
    for x,i in enumerate(np.unique(bondtypearray)):
        out_moltypes[1].append(f"bondID  {x+1}, atomtypes in this bond: {i}\n")  
    for x,i in enumerate(np.unique(angletypearray)):
        if(i not in out_moltypes[1]):
            out_moltypes[2].append(f"{i}")
    for x,i in enumerate(np.unique(dihedraltypearray)):
        if(i not in out_moltypes[2]):
            out_moltypes[3].append(f"{i}")
    for x,i in enumerate(np.unique(atomtypearray)):
        if(i not in out_moltypes[3]):
            out_moltypes[0].append(f"{i}")
    print("charges: ", chargearray)
    bondtypearray=[[x+1 for x,i in enumerate(np.unique(bondtypearray)) if (i==j)][0]  for j in bondtypearray]
    angletypearray=[[x+1 for x,i in enumerate(np.unique(angletypearray))if(i==j)][0] for j in angletypearray]
    dihedraltypearray=[[x+1 for x,i in enumerate(np.unique(dihedraltypearray))if(i==j)][0] for j in dihedraltypearray]
    for i,vert in enumerate(bm.verts):
        out_atoms.append(f"{index} {moleculeID} {atomtypearray[i]} {chargearray[i]} {vert.co[0]} {vert.co[1]} {vert.co[2]}\n")
        index+=1
    
    #WEITER HIER brauche eine Liste der Bond types über die atomsorte aus den vertex groups. Muss auch prüfen ob ab=ba
    for i,e in enumerate(inputPointMesh.data.edges):
        out_bonds[i] = f"{lastbond+i+1} {bondtypearray[i]} {e.vertices[0]+1} {e.vertices[1]+1}\n"
    print("uniqueangles: ",angletypearray, len(out_angles), ", ", len(angletypearray))
    for i,a in enumerate(out_angles):
        out_angles[i] = f"{lastangle+i+1} {angletypearray[i]} {a}\n"
    for i, d in enumerate(out_dihedrals):
        out_dihedrals[i] = f"{lastdihedral+i+1} {dihedraltypearray[i]} {d}\n"
    
    return index, out_atoms, out_bonds, out_angles, out_dihedrals, out_moltypes


    # cycle through all vertices
    
# system: as long as there is no GUI yet, every atom type contains its information with tags and numbers in the object name
# all selected objects are used to create a system
# the required information for each atom type is the charge 
# for polarizable atoms, if not created by geometry, the prefix Pol is used, with the name being double for the core-shell pair as used in the LAMMPS CORESHELL module
# Before running make sure to apply all relevant modifiers such as arrays. 
# If you represent the atoms with a skin modifier, delete it as the mesh vertices are treated as atoms
def WriteLammpsFile(inputSelection):
    atoms = [] # default array for all atom lines written as strings
    polarizable = []
    bonds = []
    angles = []
    dihedrals = []
    CSinfo = [] # default array for polarizable pairs written as strings
    #print(len(inputSelection))
    last_bond=0
    last_angle=0
    last_dihedral=0
    n_bondTypes=0
    # first iterate and write all atom to arrays
    last_atom = 1 # index for atom counter
    moleculeID = 1 # the molecule ID needs to be changed if polarizable two-point dipoles are used. Therefore, the ID must be counted up.
    # nonpolarizable atoms still use the molecule-ID 1
    
    inputList = []
    out_moltypes=[[],[],[],[]] # bonds, angles,dihedrals,atom types
    for obj in inputSelection:
        n = obj.name.split("_")
        if("Pol" in n):
            inputList.insert(0,obj)
        else:
            inputList.append(obj)
    #print(inputList)
    #print(inputSelection)
    for obj in inputList:
        #print(obj.name)
        nam = obj.name.split("_")
        if ("Pol" in nam):
            
            lastatom, o_at, o_bd, o_cs, moleculeID,out_moltypes =MakePolarizable(obj, moleculeID+1,  last_bond,out_moltypes)
            last_bond+=len(o_bd)
            for i in o_at:
                polarizable.append(i)
            for i in o_bd:
                bonds.append(i)
            for i in o_cs:
                CSinfo.append(i)
            
            n_bondTypes+=1
        elif ("Mol" in nam):
            last_atom, o_at, o_bd, o_angles, o_dihedral, out_moltypes_update= MakeMolecule(obj, moleculeID, last_atom, last_bond, last_angle, last_dihedral, out_moltypes)
            for i in out_moltypes_update[0]:
                if i not in out_moltypes[0]:
                    out_moltypes[0].append(i)
            for i in out_moltypes_update[1]:
                if i not in out_moltypes[1]:
                    out_moltypes[1].append(i)
            for i in out_moltypes_update[2]:
                if i not in out_moltypes[2]:
                    out_moltypes[2].append(i)
            for i in out_moltypes_update[3]:
                if i not in out_moltypes[3]:
                    out_moltypes[3].append(i)
            moleculeID+=1
            last_bond+=len(o_bd)
            last_angle+=len(o_angles)
            last_dihedral+=len(o_dihedral)
            for i in o_at:
                atoms.append(i)
            for i in o_bd:
                bonds.append(i)
            for i in o_angles:
                angles.append(i)
            for i in o_dihedral:
                dihedrals.append(i)
        else:
            print(nam)
            last_atom, o_at,out_moltypes = MakeAtoms(obj, last_atom,out_moltypes)
            #print(o_at) 
            for i in o_at:
                atoms.append(i)
    # second iteration: Open the file and write everything to add the header 
    #print(len(atoms)) 
    f = open(datapath+f"\{datafilename}_.data", "w")
    f.write(f"#This structure was created in Blender, plugin written by Rene Rekers\n\n")
    f.write(f"{len(atoms)+len(polarizable)} atoms\n{len(bonds)} bonds\n{len(angles)} angles\n{len(dihedrals)} dihedrals\n\n")
    f.write(f"{len(out_moltypes[0])} atom types\n{len(out_moltypes[1])} bond types\n{len(out_moltypes[2])} angle types\n{len(out_moltypes[3])} dihedral types\n\n")
    f.write(f"{-10} {60} xlo xhi\n")
    f.write(f"{-25} {25} ylo yhi\n")
    f.write(f"{-25} {25} zlo zhi\n\nAtoms # full\n\n")
    for i in polarizable:
        f.write(i)
    for i in atoms:
        f.write(i)
    if(len(bonds)>0):
        f.write("\nBonds\n\n")
        for i in bonds:
            f.write(i)
    if(len(angles)>0):
        f.write("\nAngles\n\n")
        for i in angles:
            f.write(i)
        
    if(len(dihedrals)>0):
        f.write("\nDihedrals\n\n")
        for i in dihedrals:
            f.write(i)
    if(len(CSinfo)>0):
        f.write("\nCS-Info\n\n")
        for i in CSinfo:
            f.write(i)
    f.close()
    
    #generate helper file
    # preprocess all info and generate the type file
    f = open(datapath+f"\{datafilename}_helper.txt","w")
    f.write("This file helps you with the force field parameters\n\n")    
    f.write("bond types\n")
    # sort the losts for uniques
    for x,w in enumerate(np.unique(out_moltypes[1])):
        f.write(f"bond{x+1}, atomtypes in this bond: "+w)
    f.write("\nangle types\n")
    # sort the losts for uniques
    for x,w in enumerate(np.unique(out_moltypes[2])):
        f.write(f"angleID {x+1}, atomtypes in this angle: "+w+"\n")
    f.write("\ndihedral types\n")
    for x,w in enumerate(np.unique(out_moltypes[3])):
        f.write(f"dihedralID {x+1}, atomtypes in this dihedral: "+w+"\n")
    f.close()
# 2*8.49051*cos(pi/6)
# 20.78957

# vs LTP: 8.5135 and 20.8605


# get the active object from the scene
active_object = bpy.context.view_layer.objects.active

# ensure that there is an active object
if (active_object != None):
    # call the method to select all vertices in bounding box vectors
    # object             the object to operate on
    # vector 1           the first coordinate of the bounding box
    # vector 2           the second coordinate of the bounding box
    #SelectVerticesInBound(active_object, Vector((10.63, -11.05,  -20.8)), Vector((-6.22853, 3.55585, -0.17 )))
    #print(bpy.context.selected_objects)
    WriteLammpsFile(bpy.context.selected_objects)
