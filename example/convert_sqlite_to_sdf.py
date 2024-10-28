from yammbs import MoleculeStore
from openff.toolkit.topology import Molecule
from openff.units import unit
from rdkit import Chem


store = MoleculeStore('OpenFF_Gen2_Coverage_sage220.sqlite')
ff_yammbs = 'openff_unconstrained-2.2.0.offxml'
mol_ids = store.get_molecule_ids()
qmdir = 'OpenFF_Gen2_Coverage_sage220_QM'
mmdir = 'OpenFF_Gen2_Coverage_sage220_MM'

all_confs = 0
for i in mol_ids:
    qca_ids = store.get_qcarchive_ids_by_molecule_id(i)
    all_confs += len(qca_ids)
    for j,qca_id in enumerate(qca_ids):

        mapped_smiles = store.get_smiles_by_molecule_id(i)
        molfile = 'mol-{:02}-conf-{:02}.sdf'.format(i,j)

        mm_conf = store.get_mm_conformers_by_molecule_id(i,force_field = ff_yammbs)[j]
        qm_conf = store.get_qm_conformers_by_molecule_id(i)[j]

        mm_mol = Molecule.from_mapped_smiles(mapped_smiles,allow_undefined_stereo=True)
        mm_mol.add_conformer(unit.Quantity(mm_conf,'angstrom'))
        mm_mol_rdkit = mm_mol.to_rdkit()

        writer = Chem.SDWriter('{}/{}'.format(mmdir,molfile))
        mm_mol_rdkit.SetProp('Record QCArchive',str(qca_id))
        mm_mol_rdkit.SetProp('Mapped SMILES',str(mapped_smiles)) # Just doing this for consistent ordering
        writer.write(mm_mol_rdkit)

        qm_mol = Molecule.from_mapped_smiles(mapped_smiles,allow_undefined_stereo=True)
        qm_mol.add_conformer(unit.Quantity(qm_conf,'angstrom'))
        qm_mol_rdkit = qm_mol.to_rdkit()

        writer = Chem.SDWriter('{}/{}'.format(qmdir,molfile))
        qm_mol_rdkit.SetProp('Record QCArchive',str(qca_id))
        qm_mol_rdkit.SetProp('Mapped SMILES',str(mapped_smiles))
        writer.write(qm_mol_rdkit)

print(all_confs)
