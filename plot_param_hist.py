import pandas as pd
import json
import numpy as np
import tqdm
from openff.toolkit import ForceField
import sys
import matplotlib.pyplot as plt
import click
import os

def plot_hist(datasets,labels,filename,title='',xlim=None,ylim=None,legend=True,lw=1,xlab = ''):
    plt.figure()
    plt.axvline(x=0,linestyle='--',color='k',linewidth=0.5)
    for i,data in enumerate(datasets):
        plt.hist(data,histtype='step',bins=50,label=labels[i],linewidth=lw)

    if xlim != None: plt.xlim(xlim[0],xlim[1])
    plt.title(title)
    plt.ylabel('Count')
    plt.xlabel(xlab)
    if legend: plt.legend()
    if ylim != None: plt.ylim(ylim[0],ylim[1])
    plt.savefig(filename)
    plt.close()

def boxplot(datasets,labels,filename,title='',ylab = '',ylim = None):

    fig,ax = plt.subplots(nrows=1,ncols=1)
    ax.boxplot(datasets,tick_labels=labels)
    left,right=ax.get_xlim()
    ax.hlines(y=0,linestyle='--',color='k',linewidth=0.5,xmin=left,xmax=right)
    ax.set_title(title)
    ax.set_ylabel(ylab)
    if ylim != None: ax.set_ylim(ylim[0],ylim[1])
    plt.savefig(filename)
    plt.close()

def per_param_err(mm,qm,dtype):
    if dtype == 'bonds' or dtype == 'angles':
        return mm - qm
    else: # For torsions, need to make sure they're aligned correctly
        mm_min_qm = mm - qm

        # try shifting with minimum image convention
        x_size = 360

        mm_min_qm[mm_min_qm < -x_size*0.5] += x_size
        mm_min_qm[mm_min_qm >= x_size*0.5] -= x_size

        return mm_min_qm


def get_errs(param_list,suffixs,dtype,): # dtype = bonds, angles, propers, impropers
    # Calculate the error in each parameter of type dtype
    param_errs = {}
    for l,suffix in enumerate(suffixs):
        param_errs[suffix]={param : {'errors': []} for param in param_list}

        PARAM_DATA_FILE='{}/{}_qmv{}.json'.format(suffix,dtype,suffix.split('/')[-1])

        with open(PARAM_DATA_FILE,'r') as jsonfile:
            PARAM_JSON = dict(json.load(jsonfile))

        # Would be better to just loop these once, but the smirks changes
        # for each FF so data would get lost
        for smirks ,value in PARAM_JSON.items():
            param_id = value['ident']
            if param_id in param_list:
                param_errs[suffix][param_id]['errors']=per_param_err(np.array(value['mm_values'] ),np.array(value['qm_values']),dtype)
                param_errs[suffix][param_id]['smirks'] = smirks

    return param_errs

def plot_errs(param_errs,param_list,suffixs,dtype,save_dir,labels,ylim=None,xlim=None):
    # Plot histogram and boxplot of the error.

    fig_titles = {'bonds': 'Bond length error MM - QM (A)', 'angles': "Angle error MM - QM (deg)", 'propers': 'Proper torsion error MM - QM (deg)','impropers': 'Improper torsion error MM - QM (deg)'}

    for i,param in enumerate(param_list):
        # Need to have several try/excepts, in case the different FFs don't have the same SMIRKS. Could extend deeper
        try:
            param_errs_x =[param_errs[x][param]['errors'] for x in suffixs]

            try:
                smirks = param_errs[suffixs[-1]][param]['smirks']
            except KeyError:
                if len(suffixs) > 1:
                    try:
                        smirks = param_errs[suffixs[-2]][param]['smirks']
                    except KeyError:
                        smirks = ''
                else:
                    smirks = ''

            plot_hist(param_errs_x,labels,filename='{}/{}_hist.pdf'.format(save_dir,param),title=param + ' ' + smirks,xlab=fig_titles[dtype],ylim=ylim,xlim=xlim) # Generalize title + allow passing limits etc
            boxplot(param_errs_x,labels,filename='{}/{}_boxplot.pdf'.format(save_dir,param),title=param + ' ' + smirks ,ylab=fig_titles[dtype],ylim=ylim )

        except KeyError:
            print(param)
            pass


@click.command()
@click.option('--data_dir',multiple=True,help='Directory where data jsons are')
@click.option('--labels',multiple=True,default=[],help='Short label for each data directory to include in legend. If not provided, will use the data_dir')
@click.option('--ff_file',default='',help='Force field whose parameters should be used for grouping (optional)')
@click.option('--bonds',multiple=True,default=[],help="Parameter id of bonds to plot. If 'all', will plot all bonds defined by ff_file")
@click.option('--angles',multiple=True,default=[],help="Parameter id of angles to plot. If 'all', will plot all angles defined by ff_file")
@click.option('--propers',multiple=True,default=[],help="Parameter id of proper torsions to plot. If 'all', will plot all proper torsions defined by ff_file")
@click.option('--impropers',multiple=True,default=[],help="Parameter id of improper torsions to plot. If 'all', will plot all improper torsions defined by ff_file")
@click.option('--save_dir',help='Directory to save the data')
def main(data_dir,labels,ff_file,bonds,angles,propers,impropers,save_dir):
    suffixs = data_dir
    if len(labels) == 0:
        labels = data_dir

    try:
        os.mkdir(save_dir)
    except FileExistsError:
        pass

    # list all parameters from the ff file, if all parameters requested
    if 'all' in bonds:
        bonds = [param.id for param in ForceField(ff_file).get_parameter_handler("Bonds").parameters]
    if 'all' in angles:
        angles = [param.id for param in ForceField(ff_file).get_parameter_handler("Angles").parameters]
    if 'all' in propers:
        propers = [param.id for param in ForceField(ff_file).get_parameter_handler("ProperTorsions").parameters]
    if 'all' in impropers:
        impropers = [param.id for param in ForceField(ff_file).get_parameter_handler("ImproperTorsions").parameters]

    if len(bonds) > 0:
        bond_errs = get_errs(bonds,suffixs,'bonds')
        plot_errs(bond_errs,bonds,suffixs,'bonds',save_dir,labels) # These functions do accept y/x limits, but have to pass it manually for now

    if len(angles) > 0:
        angle_errs = get_errs(angles,suffixs,'angles')
        plot_errs(angle_errs,angles,suffixs,'angles',save_dir,labels)

    if len(propers) > 0:
        proper_errs = get_errs(propers,suffixs,'propers')
        plot_errs(proper_errs,propers,suffixs,'propers',save_dir,labels)

    if len(impropers) > 0:
        improper_errs = get_errs(impropers,suffixs,'impropers')
        plot_errs(improper_errs,impropers,suffixs,'impropers',save_dir,labels)


if __name__ == '__main__':
    main()
