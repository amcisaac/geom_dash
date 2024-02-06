import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
from openff.toolkit import Molecule, ForceField, Topology
import json
import tqdm
import sys


def add_mol_to_dict(qm_file,mm_file,ff,bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict,filter=False):
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
    # filter_type = True
    # if filter and len()
    qm_mol = Molecule(qm_file,allow_undefined_stereo=True)
    qm_mol_rd = Chem.SDMolSupplier(qm_file,removeHs=False)[0]
    qm_conf = qm_mol_rd.GetConformer()

    mm_mol = Molecule(mm_file,allow_undefined_stereo=True)
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
    # print(angle_dict)

    # qm_data_dict[recordid] = {'Bonds':{},'Angles':{},'ProperTorsions':{},'ImproperTorsions':{}} # record ID: {'Bonds': (idx):qm_bl} --> can be read in to avoid recalc
    # mm_data_dict[recordid] = {'Bonds':{},'Angles':{},'ProperTorsions':{},'ImproperTorsions':{}} # same but for MM


    # Bonds
    # for idx in bond_dict:
    #     b = bond_dict[idx]
    #
    #     # This part could potentially be eliminated by reading in QM geom data
    #     qm_bl = Chem.rdMolTransforms.GetBondLength(qm_conf,idx[0],idx[1])
    #     # qm_data_dict[recordid]['Bonds'][idx] = qm_bl
    #
    #     mm_bl = Chem.rdMolTransforms.GetBondLength(mm_conf,idx[0],idx[1]) # A
    #     # mm_data_dict[recordid]['Bonds'][idx] = mm_bl
    #
    #     if b.smirks in bond_data_dict:
    #         bond_data_dict[b.smirks]['molecules'].append(smiles)
    #         bond_data_dict[b.smirks]['envs'].append(idx)
    #         bond_data_dict[b.smirks]['qm_values'].append(qm_bl)
    #         bond_data_dict[b.smirks]['mm_values'].append(mm_bl)
    #     else:
    #         bond_data_dict[b.smirks] = {'ident':b.id,
    #                                     'sage_value': b.length.magnitude,
    #                                     'qm_values': [qm_bl],
    #                                     'mm_values': [mm_bl],
    #                                     'molecules': [smiles],
    #                                     'envs': [idx]}

    # print(len(angle_dict.keys()))
    # Angles
    for idx in angle_dict.keys():
        b = angle_dict[idx]
        # print(idx)
        # print(qm_mol.chemical_environment_matches(filter))
        # print(idx in qm_mol.chemical_environment_matches(filter))
        # if filter:
        matches = qm_mol.chemical_environment_matches(filter)
        filter_condition = idx in qm_mol.chemical_environment_matches(filter)
        filter_condition = ((idx[0],) in matches or (idx[1],) in matches) or (idx[2],) in matches
        # else: filter_condition = True
        # print(filter_condition)
        if filter_condition:

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
    # for idx in proper_dict:
    #     b = proper_dict[idx]
    #
    #     # This part could potentially be eliminated by reading in QM geom data
    #     qm_tor = Chem.rdMolTransforms.GetDihedralDeg(qm_conf,idx[0],idx[1],idx[2],idx[3])
    #     # qm_data_dict[recordid]['ProperTorsions'][idx] = qm_tor
    #     mm_tor = Chem.rdMolTransforms.GetDihedralDeg(mm_conf,idx[0],idx[1],idx[2],idx[3])
    #     # mm_data_dict[recordid]['ProperTorsions'][idx] = mm_tor
    #
    #     if b.smirks in proper_data_dict:
    #         proper_data_dict[b.smirks]['molecules'].append(smiles)
    #         proper_data_dict[b.smirks]['envs'].append(idx)
    #         proper_data_dict[b.smirks]['qm_values'].append(qm_tor)
    #         proper_data_dict[b.smirks]['mm_values'].append(mm_tor)
    #     else:
    #         proper_data_dict[b.smirks] = {'ident':b.id,
    #                                     # 'sage_value': b.angle.magnitude,
    #                                     'qm_values': [qm_tor],
    #                                     'mm_values': [mm_tor],
    #                                     'molecules': [smiles],
    #                                     'envs': [idx]}
    #
    # # Improper Torsions
    # for idx in improper_dict:
    #     b = improper_dict[idx]
    #
    #     # This part could potentially be eliminated by reading in QM geom data
    #     qm_tor = Chem.rdMolTransforms.GetDihedralDeg(qm_conf,idx[0],idx[1],idx[2],idx[3])
    #     # qm_data_dict[recordid]['ProperTorsions'][idx] = qm_tor
    #     mm_tor = Chem.rdMolTransforms.GetDihedralDeg(mm_conf,idx[0],idx[1],idx[2],idx[3])
    #     # mm_data_dict[recordid]['ProperTorsions'][idx] = mm_tor
    #
    #     if b.smirks in improper_data_dict:
    #         improper_data_dict[b.smirks]['molecules'].append(smiles)
    #         improper_data_dict[b.smirks]['envs'].append(idx)
    #         improper_data_dict[b.smirks]['qm_values'].append(qm_tor)
    #         improper_data_dict[b.smirks]['mm_values'].append(mm_tor)
    #     else:
    #         improper_data_dict[b.smirks] = {'ident':b.id,
    #                                     # 'sage_value': b.angle.magnitude,
    #                                     'qm_values': [qm_tor],
    #                                     'mm_values': [mm_tor],
    #                                     'molecules': [smiles],
    #                                     'envs': [idx]}

    return

def main(mm_dir0, ff_file,qm_dir0,conformers=False,dir = '/Users/lexiemcisaac/Documents/OpenFF/conformer_energy_ordering/swope_scripts/benchmarking/'):
    # Collecting data for all molecules
    bond_data_dict = {}
    angle_data_dict = {}
    proper_data_dict = {}
    improper_data_dict = {}
    qm_data_dict = {}
    mm_data_dict = {}

    if conformers:
        compound_list = dir + 'all_molecules.txt'
    else:
        compound_list = dir + 'compound.list'
    qm_dir = dir + qm_dir0
    mm_dir = dir + mm_dir0
    sage = ForceField(ff_file,allow_cosmetic_attributes=True)

    all_mols = np.loadtxt(compound_list,dtype='str')#[:10]

    for mol in tqdm.tqdm(all_mols,desc='Calculating geometric parameters'):
        if not conformers:
            mol += '-00.sdf'
        qm_file = qm_dir +'/'+ mol
        mm_file = mm_dir +'/'+ mol
        try:
            add_mol_to_dict(qm_file,mm_file,sage,bond_data_dict,angle_data_dict,proper_data_dict,improper_data_dict,filter='[r4:1]')
        except OSError:
            pass

    # with open('bonds_qmv{}.json'.format(mm_dir0),'w') as jsonfile:
    #     json.dump(bond_data_dict,jsonfile,indent=4)

    with open('angles_qmv{}_r4.json'.format(mm_dir0),'w') as jsonfile:
        json.dump(angle_data_dict,jsonfile,indent=4)

    # with open('propers_qmv{}.json'.format(mm_dir0),'w') as jsonfile:
    #     json.dump(proper_data_dict,jsonfile,indent=4)
    #
    # with open('impropers_qmv{}.json'.format(mm_dir0),'w') as jsonfile:
    #     json.dump(improper_data_dict,jsonfile,indent=4)


if __name__ == '__main__':

    try: mm_dir_prefix = sys.argv[1]
    except IndexError:
        mm_dir_prefix = 'openff-2.1.0'
    try: qm_dir_prefix = sys.argv[3]

    except IndexError:
        qm_dir_prefix = 'b3lyp-d3bj_dzvp'
    try: ff = sys.argv[2]
    except IndexError:
        ff = 'openff-2.1.0.offxml'

    main(mm_dir_prefix,ff,qm_dir_prefix)
