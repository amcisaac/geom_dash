import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
from openff.toolkit import Molecule, ForceField, Topology
import json
import tqdm
import sys


def add_mol_to_dict(qm_file,mm_file,ff,bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict):
    '''
    Function to enumerate FF parameters assigned to a molecule (e.g. bond length, angle, dihedral angle),
    calculate the actual value from both QM and MM, and compile a set of parameter dictionaries for plotting.

    Inputs:
        qm_file: SDF file with the optimized QM geometry
        mm_file: SDF file with the optimized MM geometry
        ff: ForceField object to assign parameters with
        bond_data_dict: dictionary where bond lengths are enumerated. Can be empty
        angle_data_dict: dictionary where angles are enumerated. Can be empty
        proper_data_dict: dictionary where dihedral angles corresponding to
                          proper torsion parameters are enumerated. Can be empty
        improper_data_dict: dictionary where dihedral angles corresponding to
                          improper torsion parameters are enumerated. Can be empty

    Returns:
        Nothing is returned, but data dictionaries are modified in place.
    '''
    qm_mol = Molecule(qm_file)
    qm_mol_rd = Chem.SDMolSupplier(qm_file,removeHs=False)[0]
    qm_conf = qm_mol_rd.GetConformer()

    mm_mol = Molecule(mm_file)
    mm_mol_rd = Chem.SDMolSupplier(mm_file,removeHs=False)[0]
    mm_conf = mm_mol_rd.GetConformer()

    recordid = qm_mol_rd.GetPropsAsDict()['Record QCArchive']
    smiles = qm_mol.to_smiles(mapped=True)

    topo = Topology.from_molecules([qm_mol])
    molecule_force_list = ff.label_molecules(topo) # dictionary of forces for the molecule, keys are type of force ('Bonds', 'Angles', etc)
    # these have the form {(atom1 idx, atom2 idx,...): FF parameter} for all bonds/angles/torsions in the molecule
    bond_dict = dict(molecule_force_list[0]['Bonds'])
    angle_dict = dict(molecule_force_list[0]['Angles'])
    proper_dict = dict(molecule_force_list[0]['ProperTorsions'])
    improper_dict = dict(molecule_force_list[0]['ImproperTorsions'])

    # qm_data_dict[recordid] = {'Bonds':{},'Angles':{},'ProperTorsions':{},'ImproperTorsions':{}} # record ID: {'Bonds': (idx):qm_bl} --> can be read in to avoid recalc
    # mm_data_dict[recordid] = {'Bonds':{},'Angles':{},'ProperTorsions':{},'ImproperTorsions':{}} # same but for MM

    # Bonds
    for idx in bond_dict:
        b = bond_dict[idx]

        # This part could potentially be eliminated by reading in QM geom data
        qm_bl = Chem.rdMolTransforms.GetBondLength(qm_conf,idx[0],idx[1])
        # qm_data_dict[recordid]['Bonds'][idx] = qm_bl

        mm_bl = Chem.rdMolTransforms.GetBondLength(mm_conf,idx[0],idx[1]) # A
        # mm_data_dict[recordid]['Bonds'][idx] = mm_bl

        if b.smirks in bond_data_dict:
            bond_data_dict[b.smirks]['molecules'].append(smiles)
            bond_data_dict[b.smirks]['envs'].append(idx)
            bond_data_dict[b.smirks]['qm_values'].append(qm_bl)
            bond_data_dict[b.smirks]['mm_values'].append(mm_bl)
        else:
            bond_data_dict[b.smirks] = {'ident':b.id,
                                        'sage_value': b.length.magnitude,
                                        'qm_values': [qm_bl],
                                        'mm_values': [mm_bl],
                                        'molecules': [smiles],
                                        'envs': [idx]}

    # Angles
    for idx in angle_dict:
        b = angle_dict[idx]

        # This part could potentially be eliminated by reading in QM geom data
        qm_ang = Chem.rdMolTransforms.GetAngleDeg(qm_conf,idx[0],idx[1],idx[2])
        # qm_data_dict[recordid]['Angles'][idx] = qm_ang
        mm_ang = Chem.rdMolTransforms.GetAngleDeg(mm_conf,idx[0],idx[1],idx[2])
        # mm_data_dict[recordid]['Angles'][idx] = mm_ang

        if b.smirks in angle_data_dict:
            angle_data_dict[b.smirks]['molecules'].append(smiles)
            angle_data_dict[b.smirks]['envs'].append(idx)
            angle_data_dict[b.smirks]['qm_values'].append(qm_ang)
            angle_data_dict[b.smirks]['mm_values'].append(mm_ang)
        else:
            angle_data_dict[b.smirks] = {'ident':b.id,
                                        'sage_value': b.angle.magnitude,
                                        'qm_values': [qm_ang],
                                        'mm_values': [mm_ang],
                                        'molecules': [smiles],
                                        'envs': [idx]}

    # Proper Torsions
    for idx in proper_dict:
        b = proper_dict[idx]

        # This part could potentially be eliminated by reading in QM geom data
        qm_tor = Chem.rdMolTransforms.GetDihedralDeg(qm_conf,idx[0],idx[1],idx[2],idx[3])
        # qm_data_dict[recordid]['ProperTorsions'][idx] = qm_tor
        mm_tor = Chem.rdMolTransforms.GetDihedralDeg(mm_conf,idx[0],idx[1],idx[2],idx[3])
        # mm_data_dict[recordid]['ProperTorsions'][idx] = mm_tor

        if b.smirks in proper_data_dict:
            proper_data_dict[b.smirks]['molecules'].append(smiles)
            proper_data_dict[b.smirks]['envs'].append(idx)
            proper_data_dict[b.smirks]['qm_values'].append(qm_tor)
            proper_data_dict[b.smirks]['mm_values'].append(mm_tor)
        else:
            proper_data_dict[b.smirks] = {'ident':b.id,
                                        # 'sage_value': b.angle.magnitude,
                                        'qm_values': [qm_tor],
                                        'mm_values': [mm_tor],
                                        'molecules': [smiles],
                                        'envs': [idx]}

    # Improper Torsions
    for idx in improper_dict:
        b = improper_dict[idx]

        # This part could potentially be eliminated by reading in QM geom data
        qm_tor = Chem.rdMolTransforms.GetDihedralDeg(qm_conf,idx[0],idx[1],idx[2],idx[3])
        # qm_data_dict[recordid]['ProperTorsions'][idx] = qm_tor
        mm_tor = Chem.rdMolTransforms.GetDihedralDeg(mm_conf,idx[0],idx[1],idx[2],idx[3])
        # mm_data_dict[recordid]['ProperTorsions'][idx] = mm_tor

        if b.smirks in improper_data_dict:
            improper_data_dict[b.smirks]['molecules'].append(smiles)
            improper_data_dict[b.smirks]['envs'].append(idx)
            improper_data_dict[b.smirks]['qm_values'].append(qm_tor)
            improper_data_dict[b.smirks]['mm_values'].append(mm_tor)
        else:
            improper_data_dict[b.smirks] = {'ident':b.id,
                                        # 'sage_value': b.angle.magnitude,
                                        'qm_values': [qm_tor],
                                        'mm_values': [mm_tor],
                                        'molecules': [smiles],
                                        'envs': [idx]}

        return

def main():
    # Collecting data for all molecules
    bond_data_dict = {}
    angle_data_dict = {}
    proper_data_dict = {}
    improper_data_dict = {}
    qm_data_dict = {}
    mm_data_dict = {}

    # TODO: make this modifiable easily
    dir = '/Users/lexiemcisaac/Documents/OpenFF/conformer_energy_ordering/swope_scripts/benchmarking/'
    compound_list = dir + 'all_molecules.txt'
    qm_dir = dir + 'b3lyp-d3bj_dzvp/'
    mm_dir = dir + 'openff-2.1.0/'
    sage = ForceField('openff-2.1.0.offxml')
    # jsonpath = 'test_bond2.json'

    all_mols = np.loadtxt(compound_list,dtype='str')[:100]
    # print(all_mols)
    for mol in tqdm.tqdm(all_mols,desc='Calculating geometric parameters'):
        qm_file = qm_dir + mol
        mm_file = mm_dir + mol
        add_mol_to_dict(qm_file,mm_file,sage,bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict)

    with open('bonds_test2.json','w') as jsonfile:
        json.dump(bond_data_dict,jsonfile,indent=4)

    with open('angles_test2.json','w') as jsonfile:
        json.dump(angle_data_dict,jsonfile,indent=4)

    with open('propers_test2.json','w') as jsonfile:
        json.dump(proper_data_dict,jsonfile,indent=4)

    with open('impropers_test2.json','w') as jsonfile:
        json.dump(improper_data_dict,jsonfile,indent=4)


if __name__ == '__main__':
    main()
