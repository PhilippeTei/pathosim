import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt 
import sciris as sc 
import copy

# Constants
PCR_criteria_defs = {'cont_conf': 'Recent contact with PCR or RAT confirmed case',
                    'cont_vuln': 'Has household members aged 60+',
                    'cont_sx': 'Contact with symptomatic agent (COVID or ILI)',
                    'cont_cs_sx': 'Contact with symptomatic individual (fever AND cough or sore throat)',
                    'cont_ncs_sx': 'Contact with symptomatic individual (general ILI symptoms)',
                    'cs_sx': 'Has fever AND cough or sore throat', 
                    'ncs_sx': 'Has general ILI symptoms',
                    'pos_RAT': 'Recent positive RAT',
                    'pos_RAT_asx': 'Recent positive RAT without symptoms',
                    'neg_RAT_sx': 'Recent negative RAT with symptoms',
                    'other': 'Other reasons',
                    'work': 'Mandatory workplace testing'}

RAT_criteria_defs = {'cont_conf': 'Recent contact with PCR or RAT confirmed case',
                    'cont_vuln': 'Has household members aged 60+',
                    'cont_sx': 'Contact with symptomatic agent (COVID or ILI)',
                    'cont_cs_sx': 'Contact with symptomatic individual (fever AND cough or sore throat)',
                    'cont_ncs_sx': 'Contact with symptomatic individual (general ILI symptoms)',
                    'cs_sx': 'Has fever AND cough or sore throat', 
                    'ncs_sx': 'Has general ILI symptoms',
                    'pos_RAT': 'Recent positive RAT',
                    'neg_RAT': 'Recent negative RAT',
                    'pos_RAT_asx': 'Recent positive RAT without symptoms',
                    'neg_RAT_sx': 'Recent negative RAT with symptoms',
                    'other': 'Other reasons',
                    'work': 'Mandatory workplace testing'}

default_criteria = ['cont_conf', 'cont_vuln', 'cont_cs_sx', 'cont_ncs_sx', 'cs_sx', 'ncs_sx', 'pos_RAT', 'pos_RAT_asx', 'neg_RAT_sx', 'neg_RAT', 'other', 'work']

default_colors = [(0.25036271661125314, 0.9831770085338173, 0.3349247522917396),
 (0.2328287334209082, 0.30569880709893255, 0.7696900575827298),
 (0.3453211223438526, 0.8351921754410631, 0.9870821874532515),
 (0.9110926168081614, 0.8801156753945002, 0.24030899803089276),
 (0.6902410934441348, 0.2777056732325177, 0.26502620220400885),
 (0.4805119696395338, 0.6489481769846063, 0.5548401661670068),
 (0.6928470904849535, 0.24596124778484157, 0.836912769852372),
 (0.7214622676177648, 0.9756193515412254, 0.6465881227849036),
 (0.24587914941518646, 0.34272226584755033, 0.2553666082229418),
 (0.9927317230643592, 0.4889623746052707, 0.271020403451419),
 (0.621404442324434, 0.561305132515773, 0.9770518020000377)]

letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

# Functions for customizing the y-scale
def forward(x): 
    return x**(1/2)

def inverse(x): 
    return x**2


def plot_eligible(testobj, scale='root', style=None, dpi=400, legend='on',lw=2, xlim=None,
                ylim=None, fig=None, ax=None, show_criteria=default_criteria, colors=default_colors, fname=None): 

    # Set some basic parameters
    if 'PCR' in testobj.system: 
        criteria_defs = PCR_criteria_defs
    else: 
        criteria_defs = RAT_criteria_defs

    if style == "science": 
        plt.style.use(["science", "nature", "grid", "no-latex"])
    elif not(style is None): 
        raise RuntimeError("Invalid style provided.")
    
    # Generate plot canvas
    fig, ax = plt.subplots(figsize=(6, 4))

    fig.suptitle('Number of individuals meeting each criterion',y=0.94, fontsize=13)
    ax.set_xlabel('Days since beginning of outbreak', fontsize=11)
    ax.set_ylabel('Number of individuals', fontsize=11)

    ax.set_xlim((0, 200))
    if not(xlim is None):
        ax.set_xlim(xlim)
    if not(ylim is None): 
        ax.set_ylim(ylim)
    ax.set_xticks(np.arange(0, 201, 50))
    ax.tick_params(labelsize=11)

    if scale == 'root':
        ax.set_yscale('function', functions=(forward, inverse))
    elif scale == 'linear': 
        ax.set_yscale('linear')
    
    # Get valid keys
    valid = [key for key in show_criteria if key in testobj.alloc_groups.keys()]

    # Plot 
    l = 0
    for key, color in zip(valid, colors): 
        if key in testobj.crit_groups.keys(): 
            ax.plot(testobj.crit_groups[key], color=color, label=letters[l] + ": " + criteria_defs[key], alpha=0.7, lw=lw)
            l += 1 

    # Make the legend
    if legend != 'off': 
        lines = []
        labels = []
        for ax in fig.axes:
            axLine, axLabel = ax.get_legend_handles_labels()
            lines.extend(axLine)
            labels.extend(axLabel)
        
        fig.legend(lines, labels, loc='center left', fontsize=7, labelspacing=0.8, bbox_to_anchor=(0.91, 0.5), frameon=False)

    # Save 
    if not(fname is None): 
        fig.savefig(fname, dpi=fig.dpi, bbox_inches='tight')

    # Reset style 
    mpl.rcParams.update(mpl.rcParamsDefault)

    return fig, ax


def plot_seek(testobj, scale='root', style=None, dpi=400, legend='on',lw=2, xlim=None,
                ylim=None, fig=None, ax=None, show_criteria=default_criteria, colors=default_colors, fname=None): 

    # Set some basic parameters
    if 'PCR' in testobj.system: 
        criteria_defs = PCR_criteria_defs
    else: 
        criteria_defs = RAT_criteria_defs

    if style == "science": 
        plt.style.use(["science", "nature", "grid", "no-latex"])
    elif not(style is None): 
        raise RuntimeError("Invalid style provided.")
    
    # Generate plot canvas
    fig, ax = plt.subplots(figsize=(6, 4))

    fig.suptitle('Number of individuals meeting each criterion',y=0.94, fontsize=13)
    ax.set_xlabel('Days since beginning of outbreak', fontsize=11)
    ax.set_ylabel('Number of individuals', fontsize=11)

    ax.set_xlim((0, 200))
    if not(xlim is None):
        ax.set_xlim(xlim)
    if not(ylim is None): 
        ax.set_ylim(ylim)
    ax.set_xticks(np.arange(0, 201, 50))
    ax.tick_params(labelsize=11)

    if scale == 'root':
        ax.set_yscale('function', functions=(forward, inverse))
    elif scale == 'linear': 
        ax.set_yscale('linear')
    
    # Get valid keys
    valid = [key for key in show_criteria if key in testobj.alloc_groups.keys()]

    # Plot 
    l = 0
    for key, color in zip(valid, colors): 
        if key in testobj.seek_groups.keys(): 
            ax.plot(testobj.seek_groups[key], color=color, label=letters[l] + ": " + criteria_defs[key], alpha=0.7, lw=lw)
            l += 1 

    # Make the legend
    if legend != 'off': 
        lines = []
        labels = []
        for ax in fig.axes:
            axLine, axLabel = ax.get_legend_handles_labels()
            lines.extend(axLine)
            labels.extend(axLabel)
        
        fig.legend(lines, labels, loc='center left', fontsize=7, labelspacing=0.8, bbox_to_anchor=(0.91, 0.5), frameon=False)

    # Save 
    if not(fname is None): 
        fig.savefig(fname, dpi=fig.dpi, bbox_inches='tight')

    # Reset style 
    mpl.rcParams.update(mpl.rcParamsDefault)

    return fig, ax


def plot_allocate(testobj, scale='root', style=None, lw=2, dpi=400, legend='on', xlim=None, ylim=None, fig=None, ax=None,
                show_criteria=default_criteria, colors=default_colors, show_capacity=True, show_consumed=True, show_positives=True, fname=None):  


    if 'PCR' in testobj.system: 
        criteria_defs = PCR_criteria_defs
    else: 
        criteria_defs = RAT_criteria_defs

    if style == "science":  # Default to science
        plt.style.use(["science", "nature", "grid", "no-latex"])
    elif not(style is None): 
        raise RuntimeError("Invalid style provided. Please use 'science' or leave this argument empty to use default.")

    fig, ax = plt.subplots(figsize=(6, 4))

    # Set labels
    fig.suptitle('Number of test recipients meeting each criterion',y=0.94, fontsize=13)
    ax.set_xlabel('Days since beginning of outbreak', fontsize=11)
    ax.set_ylabel('Number of test recipients', fontsize=11)

    # Set axis limits 
    ax.set_xlim((0, 200))
    if not(xlim is None):
        ax.set_xlim(xlim)
    if not(ylim is None): 
        ax.set_ylim(ylim)
    ax.set_xticks(np.arange(0, 201, 50))
    ax.tick_params(labelsize=11)

    # Set y-scaling 
    if scale == 'root':
        ax.set_yscale('function', functions=(forward, inverse))
    elif scale == 'linear': 
        ax.set_yscale('linear')

    # Get valid keys (set implementation does not retain order)
    valid = [key for key in show_criteria if key in testobj.alloc_groups.keys()]

    # Plot 
    l = 0
    for key, color in zip(valid, colors): 
        if key in testobj.alloc_groups.keys(): 
            ax.plot(testobj.alloc_groups[key],color=color, label=letters[l] + ": " + criteria_defs[key], alpha=0.7, lw=lw)
        else: 
            print(key, "is not a valid criterion")
        l += 1 

    if show_capacity: 
        ax.plot(testobj.capacity, label=letters[l] + ': Test capacity', color='black', alpha=1.0, lw=lw, ls='dotted', dash_capstyle='round')
        l += 1

    if show_consumed: 
        ax.plot(testobj.tests_consumed, label=letters[l] + ': Tests consumed', color='green', alpha=1.0, lw=lw)
        l += 1
    
    if show_positives:
        ax.plot(testobj.pos_history, label=letters[l] + ': Positive PCR tests', color='red', alpha=1.0, lw=lw)
        l += 1

    # Make the legend
    if legend != 'off': 
        lines = []
        labels = []
        for ax in fig.axes:
            axLine, axLabel = ax.get_legend_handles_labels()
            lines.extend(axLine)
            labels.extend(axLabel)
        
        fig.legend(lines, labels, loc='center left', fontsize=7, labelspacing=0.8, bbox_to_anchor=(0.91, 0.5), frameon=False)

    # Save 
    if not(fname is None): 
        fig.savefig(fname, dpi=fig.dpi, bbox_inches='tight')

    # Reset style 
    mpl.rcParams.update(mpl.rcParamsDefault)

    return fig, ax 