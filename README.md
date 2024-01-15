# LAMMPS-blender-tools

In this project, python code is written to enable using blender as a molecular editor to prepare, inspect and modify LAMMPS data files.

# workflow
## object structure
Different types of materials require different types of simulations which require different types of procedures.
To specify the type of the object, the object name is used as an input. It follows the general scheme:
NAME(arbitrary)_TYPE_ 
_ is used as a separator and must not be used in the name
supported types so far:
-### Mol: Molecule
The bonds making up the molecular topology are simply represented by edges forming a mesh. All atom types are assigned using a vertex group
The vertex group also follows a naming convention:
AtomID_charge_OWNNAME
- the AtomID is used for the LAMMPS simulation
- the charge is written to the Atoms section
- the OWNNAME can be chosen freely to help assigning the atoms in the editor. Not yet reassigned when importing a file
Make sure to assign all atoms in the molecule to IDs
Also separate all loose parts in the molecule before exporting. If you still need to select all small molecules, just select all objects and go with the vertex group selection

# planned in the future

- [ ] generate an armature for polymer chains to properly reorient them
- [ ] create an UI for all scripts to create a proper blender addon
- [ ] automatically separate all single molecules in the importer
- [ ] write some checking code which ensures that there are no unassigned vertices in the mesh
