# Import packages
from dash import Dash, html, dash_table, dcc, callback, Output, Input, ctx,exceptions
import pandas as pd
import plotly.express as px
import json
import plotly.graph_objects as go
import numpy as np
import tqdm
# For drawing:
from openff.toolkit.topology import FrozenMolecule
from openff.toolkit import Molecule
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem.rdDepictor import Compute2DCoords
from IPython.display import SVG
import base64
import sys


DEFAULT_WIDTH = int(980/2)
DEFAULT_HEIGHT = int(400/2)
def draw_molecule_svg(
    molecule,
    image_width = DEFAULT_WIDTH,
    image_height= DEFAULT_HEIGHT,
    highlight_atoms= None,
    highlight_bonds = None,
    atom_notes = None,
    bond_notes = None,
    explicit_hydrogens = None,
    color_by_element = None,
):
    """Draw a molecule--stolen from OpenFF docs
    Parameters
    ==========
    molecule
        The molecule to draw.
    image_width
        The width of the resulting image in pixels.
    image_height
        The height of the resulting image in pixels.
    highlight_atoms
        A list of atom indices to highlight, or a map from indices to colors.
        Colors should be given as triplets of floats between 0.0 and 1.0.
    highlight_atoms
        A list of pairs of atom indices indicating bonds to highlight, or a map
        from index pairs to colors. Colors should be given as triplets of floats
        between 0.0 and 1.0.
    atom_notes
        A map from atom indices to a string that should be printed near the
        atom.
    bond_notes
        A map from atom index pairs to a string that should be printed near the
        bond.
    explicit_hydrogens
        If ``False``, allow uncharged monovalent hydrogens to be hidden. If
        ``True``, make all hydrogens explicit. If ``None``
    color_by_element
        If True, color heteroatoms according to their element; if False, the
        image will be black and white. By default, uses black and white when
        highlight_atoms or highlight_bonds is provided, and color otherwise.
    Raises
    ======
    KeyError
        When an atom or bond in highlight_atoms or highlight_bonds is missing
        from the image, including when it is present in the molecule but hidden.
    """
    if isinstance(molecule, FrozenMolecule):
        rdmol = molecule.to_rdkit()
    else:
        rdmol = molecule

    if color_by_element is None:
        color_by_element = highlight_atoms is None and highlight_bonds is None

    if explicit_hydrogens is None:
        idx_map = {i: i for i in range(rdmol.GetNumAtoms())}
    elif explicit_hydrogens:
        idx_map = {i: i for i in range(rdmol.GetNumAtoms())}
        rdmol = Chem.AddHs(rdmol, explicitOnly=True)
    else:
        idx_map = {
            old: new
            for new, old in enumerate(
                a.GetIdx()
                for a in rdmol.GetAtoms()
                if a.GetAtomicNum() != 1 and a.GetMass() != 1
            )
        }
        rdmol = Chem.RemoveHs(rdmol, updateExplicitCount=True)

    if highlight_atoms is None:
        highlight_atoms = []
        highlight_atom_colors = None
    elif isinstance(highlight_atoms, dict):
        highlight_atom_colors = {
            idx_map[i]: tuple(c) for i, c in highlight_atoms.items() if i in idx_map
        }
        highlight_atoms = list(highlight_atoms.keys())
    else:
        highlight_atoms = [idx_map[i] for i in highlight_atoms if i in idx_map]
        highlight_atom_colors = None

    if highlight_bonds is None:
        highlight_bonds = []
        highlight_bond_colors = None
    elif isinstance(highlight_bonds, dict):
        highlight_bond_colors = {
            rdmol.GetBondBetweenAtoms(idx_map[i_a], idx_map[i_b]).GetIdx(): tuple(v)
            for (i_a, i_b), v in highlight_bonds.items()
            if i_a in idx_map and i_b in idx_map
        }

        highlight_bonds = list(highlight_bond_colors.keys())
    else:
        highlight_bonds = [
            rdmol.GetBondBetweenAtoms(idx_map[i_a], idx_map[i_b])
            for i_a, i_b in highlight_bonds
            if i_a in idx_map and i_b in idx_map
        ]
        highlight_bond_colors = None

    if bond_notes is not None:
        for (i_a, i_b), note in bond_notes.items():
            if i_a not in idx_map or i_b not in idx_map:
                continue
            rdbond = rdmol.GetBondBetweenAtoms(idx_map[i_a], idx_map[i_b])
            rdbond.SetProp("bondNote", note)

    if atom_notes is not None:
        for i, note in atom_notes.items():
            if i not in idx_map:
                continue
            rdatom = rdmol.GetAtomWithIdx(idx_map[i])
            rdatom.SetProp("atomNote", note)

    Compute2DCoords(rdmol)

    drawer = Draw.MolDraw2DSVG(image_width, image_height)

    draw_options = drawer.drawOptions()
    if not color_by_element:
        draw_options.useBWAtomPalette()

    drawer.DrawMolecule(
        rdmol,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=highlight_atom_colors,
        highlightBonds=highlight_bonds,
        highlightBondColors=highlight_bond_colors,
    )
    drawer.FinishDrawing()

    svg_contents = drawer.GetDrawingText()

    return svg_contents

@callback(
   Output('click-output', 'children',allow_duplicate=True),
   Input('controls-and-graph','clickData'),
   Input('param-id',component_property='value'),
   Input('data_type',component_property='value'),
   prevent_initial_call=True
)
def draw_mols(clickData,fig_title,dtype):
    global JSONS
    json_file = JSONS[dtype]
    key = fig_title.split(': ')[1]
    try:
        idx = clickData['points'][0]['pointIndex']
        atom_idx = json_file[key]['envs'][idx]
        smiles = json_file[key]['molecules'][idx]

        svg = draw_molecule_svg(Molecule.from_mapped_smiles(smiles),highlight_atoms=atom_idx)

        # Took this part from Brent
        try:
            encoded = base64.b64encode(bytes(svg, "utf-8"))
        except Exception as e:
            print("error: ", e)
        pic = html.Img(src=f"data:image/svg+xml;base64,{encoded.decode()}")
        return pic
    except TypeError:
        pass

@callback(
    Output(component_id = 'controls-and-graph',component_property='figure'),
    Input(component_id = 'param-id',component_property='value'),
    Input(component_id = 'data_type',component_property='value')
)
def make_figure(fig_title,dtype,colors='blue'):
    global JSONS
    json = JSONS[dtype]
    key = fig_title.split(': ')[1]
    df = pd.DataFrame(json[key])
    if dtype == 'Bonds':
        unit = 'A'
    else:
        unit = 'deg'

    layout = go.Layout(
        title=fig_title,
        xaxis=dict(
            title="MM value ({})".format(unit)
        ),
        yaxis=dict(
            title="QM value ({})".format(unit)
    ) )

    fig = go.Figure(layout=layout)
    fig.add_trace(go.Scattergl(x = [df['qm_values'].min(),df['qm_values'].max()],y = [df['qm_values'].min(),df['qm_values'].max()],mode='lines',line=dict(color='black',dash='dash'),name='x=y'))
    fig.add_trace(go.Scattergl(x = df['mm_values'], y = df['qm_values'],mode='markers',marker=dict(color=colors),name='Data'))
    try:
        fig.add_trace(go.Scattergl(x = [df['sage_value'][0]], y = [df['sage_value'][0]],mode='markers', marker=dict(size=[15],color='green'),marker_symbol='star',name='Sage'))
    except KeyError:
        pass

    return fig

# Borrowed from Brent
@callback(
    Output("controls-and-graph",component_property='figure', allow_duplicate=True),
    [Input("submit", "n_clicks"), Input("smirks_input", "value")],
    Input('data_type','value'),
    Input('param-id','value'),
    prevent_initial_call=True,
)
def submit_smirks(_, smirks,dtype,fig_title):
    global JSONS
    if ctx.triggered_id == "submit":
        key = fig_title.split(': ')[1]
        rec = JSONS[dtype][key]
        colors = []
        for m, e in tqdm.tqdm(
            zip(rec['molecules'], rec['envs']),
            desc="Labeling molecules",
            total=len(rec['molecules']),
        ):
            mol = Molecule.from_mapped_smiles(m, allow_undefined_stereo=True)

            if (env := mol.chemical_environment_matches(smirks)) and (
                tuple(e) in env or tuple(e)[::-1] in env
            ):
                colors.append("red")
            else:
                colors.append("blue")
                # print(e,e[::-1],env)
        return make_figure(fig_title,dtype,colors=colors)
    raise exceptions.PreventUpdate()

######
#
# Defining the global inputs
#
######

suffix = sys.argv[1]
try: port = sys.argv[2]
except IndexError:
    port = 8050

# First load in all data so it's available easily for toggling
BOND_DATA_FILE='bonds_qmv{}.json'.format(suffix)
ANGLE_DATA_FILE='angles_qmv{}.json'.format(suffix)
PROPER_DATA_FILE='propers_qmv{}.json'.format(suffix)
IMPROPER_DATA_FILE='impropers_qmv{}.json'.format(suffix)

with open(BOND_DATA_FILE,'r') as jsonfile:
    BOND_JSON = dict(json.load(jsonfile))

with open(ANGLE_DATA_FILE,'r') as jsonfile:
    ANGLE_JSON = dict(json.load(jsonfile))

with open(PROPER_DATA_FILE,'r') as jsonfile:
    PROPER_JSON = dict(json.load(jsonfile))

with open(IMPROPER_DATA_FILE,'r') as jsonfile:
    IMPROPER_JSON = dict(json.load(jsonfile))

# Only JSONS winds up being needed as a global variable.
JSONS = {'Bonds':BOND_JSON,'Angles':ANGLE_JSON,'Proper Torsions':PROPER_JSON,'Improper Torsions':IMPROPER_JSON}
KEYS = {'Bonds':list(BOND_JSON.keys()),'Angles':list(ANGLE_JSON.keys()),'Proper Torsions':list(PROPER_JSON.keys()),'Improper Torsions':list(IMPROPER_JSON.keys())}

PARAM_IDS = {'Bonds': [JSONS['Bonds'][k]['ident'] for k in KEYS['Bonds']],
             'Angles': [JSONS['Angles'][k]['ident'] for k in KEYS['Angles']],
             'Proper Torsions': [JSONS['Proper Torsions'][k]['ident'] for k in KEYS['Proper Torsions']],
             'Improper Torsions': [JSONS['Improper Torsions'][k]['ident'] for k in KEYS['Improper Torsions']]}


#####
#
# Dashboard code begins
#
#####

# Set up connection between parameters available in dropdown menu and the type of data being displayed
all_label_options = {
    'Bonds': [PARAM_IDS['Bonds'][i] + ': '+KEYS['Bonds'][i] for i in range(0,len(PARAM_IDS['Bonds']))],
    'Angles': [PARAM_IDS['Angles'][i] + ': '+KEYS['Angles'][i] for i in range(0,len(PARAM_IDS['Angles']))],
    'Proper Torsions': [PARAM_IDS['Proper Torsions'][i] + ': '+KEYS['Proper Torsions'][i] for i in range(0,len(PARAM_IDS['Proper Torsions']))],
    'Improper Torsions': [PARAM_IDS['Improper Torsions'][i] + ': '+KEYS['Improper Torsions'][i] for i in range(0,len(PARAM_IDS['Improper Torsions']))]
}

for key in all_label_options.keys():
    all_label_options[key].sort()

@callback(
    Output('param-id','options'),
    Output('param-id','value'),
    Input('data_type','value')
)
def set_paramid(data_type):
    return all_label_options[data_type], all_label_options[data_type][0]

# Set up dashboard
app = Dash(__name__)
app.layout = html.Div([
    html.H1(children='QM vs MM geometry parameters'), # Title
    html.Div([  # Toggle buttons

        html.Div([
            dcc.RadioItems( # Button for data type
                    list(all_label_options.keys()),
                    value = 'Bonds',
                    id='data_type',
                    labelStyle={'display': 'inline-block', 'marginTop': '5px'}
            ),

            dcc.Dropdown( # Dropdown for parameter selection
                # options=[],
                # value=all_label_options['Bonds'][0],
                id='param-id'
            )
        ],style={'width': '49%', 'display': 'inline-block'})
    ]),
    dcc.Input(id="smirks_input", value='Enter tagged SMIRKs', style=dict(width="30vw")),
    html.Button("Submit", id="submit", n_clicks=0),

    html.Div([  # Graph/structures
        html.Div(  # Graph (will auto-update based on param ID chosen from dropdown)
            dcc.Graph(figure={}, id = 'controls-and-graph'),
            style={'width': '48%', 'padding': '0px 20px 20px 20px','display':'inline-block',}
        ),
        html.Div(  # Draw structures
            [],
            id="click-output",
            style={'width': '40%', 'padding': '0px 20px 20px 20px', 'display':'inline-block',
                "max-height": "90vh",
                "overflow": "hidden",
                "overflow-y": "scroll",
            },
        )
    ]),

])



if __name__ == '__main__':
    app.run(debug=True,port=port)
