import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt 
import sciris as sc 
import cmasher as cmr
import copy

# Constants
PCR_criteria_defs = {'cont_conf': 'Recent contact with PCR confirmed case',
                    'cont_vuln': 'Has household members aged 60+',
                    'cont_sx': 'Contact with symptomatic agent (COVID or ILI)',
                    'cont_cs_sx': 'Contact with symptomatic individual (fever AND cough or sore throat)',
                    'cont_ncs_sx': 'Contact with symptomatic individual (general ILI symptoms)',
                    'cs_sx': 'Has fever AND cough or sore throat', 
                    'ncs_sx': 'Has general ILI symptoms',
                    'pos_RAT': 'Recent positive RAT',
                    'pos_RAT_asx': 'Recent positive RAT without symptoms',
                    'neg_RAT_sx': 'Recent negative RAT with symptoms',
                    'travel': 'Mandatory travel testing',  # Not in use 
                    'work': 'Mandatory workplace testing'}

RAT_criteria_defs = {'cont_conf': 'Recent contact with PCR or RAT confirmed case',
                    'cont_vuln': 'Has household members aged 60+',
                    'cont_sx': 'Contact with symptomatic agent (COVID or ILI)',
                    'cont_cs_sx': 'Contact with symptomatic individual (fever AND cough or sore throat)',
                    'cont_ncs_sx': 'Contact with symptomatic individual (general ILI symptoms)',
                    'cs_sx': 'Has fever AND cough or sore throat', 
                    'ncs_sx': 'Has general ILI symptoms',
                    'pos_RAT': 'Recent positive RAT',
                    'pos_RAT_asx': 'Recent positive RAT without symptoms',
                    'neg_RAT_sx': 'Recent negative RAT with symptoms',
                    'travel': 'Mandatory travel testing',  # Not in use 
                    'work': 'Mandatory workplace testing'}

default_criteria = ['cont_conf', 'cont_vuln', 'cont_cs_sx', 'cont_ncs_sx', 'cs_sx', 'ncs_sx', 'pos_RAT', 'pos_RAT_asx', 'neg_RAT_sx', 'work']

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

# Functions for customizing the y-scale
def forward(x): 
    return x**(1/2)

def inverse(x): 
    return x**2


def compute_incidence_multisim_avg(sims): 
    '''
    Takes as input a list of completed simulations and returns the average incidence
    '''
    incidences = np.array([sim.results['new_infections'].values for sim in sims])
    return np.mean(incidences, axis=0)


def compute_pos_multisim_avg(slist):
    '''
    Takes as input a list of completed surveillance simulations and returns average number of positive cases each day
    Works for antigen and PCR surveillance. 
    ''' 
    pos = np.array([p.pos_history for p in slist])
    return np.mean(pos, axis=0)


# General utility functions 
def compute_PCR_multisim_avg(slist):
    '''
    Given a list of PCR surveillance systems, compute average signals
    '''
    dummy = copy.deepcopy(slist[0])

    PCR_cg = dict()
    PCR_sg = dict()
    PCR_ag = dict()

    for key in dummy.crit_groups.keys(): 
        cg = np.array([trial.crit_groups[key] for trial in slist])
        cg_avg = np.mean(cg, axis=0)
        PCR_cg[key] = cg_avg

    for key in dummy.seek_groups.keys(): 
        sg = np.array([trial.seek_groups[key] for trial in slist])
        sg_avg = np.mean(sg, axis=0)
        PCR_sg[key] = sg_avg

    for key in dummy.alloc_groups.keys(): 
        ag = np.array([trial.alloc_groups[key] for trial in slist])
        ag_avg = np.mean(ag, axis=0)
        PCR_ag[key] = ag_avg

    # Main signals
    dummy.crit_groups = PCR_cg
    dummy.seek_groups = PCR_sg
    dummy.alloc_groups = PCR_ag

    # Other signals
    dummy.tests_consumed = np.mean(np.array([trial.tests_consumed for trial in slist]), axis=0)

    return dummy



## Basic plots ## 
def plot_incidence_positives(sim, surv, y1_lim=(0, 3500), y2_lim=(0, 175), strategy='PCR', save=False):
    '''
    Plot incidence of infection against the number of positive cases. This should work for antigen and PCR surveillance objects. 
    This has not been tested yet. 
    '''
    fig, ax1 = plt.subplots(figsize=(3, 2))

    ax1.plot(sim.results['new_infections'], label='True incidence', color='black', alpha=0.7, linewidth=1.0)

    se = surv.sensitivity
    sp = surv.specificity

    ax1.set_xlabel('Days since beginning of outbreak', fontsize=7)
    ax1.set_ylabel('Number of individuals', fontsize=7)
    if strategy == 'PCR': 
        ax1.set_title('PCR strategy (SE: ' + str(se) + ', SP: ' + str(sp) + ')', fontsize=7)
    elif strategy == 'antigen': 
        ax1.set_title('Antigen strategy (SE: ' + str(se) + ', SP: ' + str(sp) + ')', fontsize=7)
    ax1.set_ylim((0, 3500))
    ax1.set_xlim(y1_lim)
    ax1.tick_params(axis='y', colors='black')
    ax1.set_xticks(np.arange(0, 201, 50))

    ax2 = ax1.twinx()
    if strategy == 'PCR': 
        ax2.plot(surv.pos_history, label='Positive PCR tests', color='red', alpha=0.8, linewidth=1.0)
    elif strategy == 'antigen':
        ax2.plot(surv.pos_history, label='Positive rapid antigen tests', color='red', alpha=0.8, linewidth=1.0)
    ax2.set_ylim(y2_lim)
    ax2.tick_params(axis='y', colors='red')

    lines = []
    labels = []

    for ax in fig.axes:
        axLine, axLabel = ax.get_legend_handles_labels()
        lines.extend(axLine)
        labels.extend(axLabel)
        
    fig.legend(lines, labels, loc=(0.605, 0.760), fontsize=5, frameon=False)

    if save:
        if strategy == 'PCR': 
            fig.savefig('PCR strategy.jpg', dpi=fig.dpi, bbox_inches='tight')
        elif strategy == 'antigen':
            fig.savefig('Antigen strategy.jpg', dpi=fig.dpi, bbox_inches='tight')
    return 

def plot_incidence_positives_AB(inc, pos_s, pos_ns, fill=False, save=True, se=0.85, sp=1.00): 
    '''
    Plot PCR incidence against true incidence of infection for broad and stringent eligibility. Not intended for antigen surveillance. 
    '''
    plt.style.use("default")
    mpl.rcParams['figure.dpi'] = 1000

    fig, ax = plt.subplots(1, 2, figsize=(6, 2.5))

    fig.supxlabel('Days since beginning of outbreak', fontsize=7, y=0.1)
    fig.suptitle('Passive surveillance based on PCR testing (SE: ' + str(se) + ', SP: ' + str(sp) + ')', fontsize=8, y=0.95)


    ax[0].plot(inc, label='Incidence of infection', color='black', alpha=0.7, lw=1.0)
    ax[0].set_ylim((0, 4000))
    ax[0].set_xlim((0, 200))
    ax[0].set_xticks(np.arange(0, 201, 50))
    ax[0].set_title('A) Broad PCR Testing Eligibility', fontsize=7)
    ax[0].set_ylabel('Number of individuals', labelpad=4, fontsize=7)

    ax01 = ax[0].twinx()
    ax01.plot(pos_ns, label='PCR-identified cases', color='red', alpha=0.7, lw=1.0)
    ax01.get_yaxis().set_visible(False)
    ax01.set_ylim((0, 200))
    # ax01.tick_params(axis='y', colors='red')

    ax[1].plot(inc, label='Incidence of infection', color='black', alpha=0.7, lw=1.0)
    ax[1].set_ylim((0, 4000))
    ax[1].set_xlim((0, 200))
    ax[1].set_xticks(np.arange(0, 201, 50))
    ax[1].set_title('B) Stringent PCR Testing Eligibility', fontsize=7)
    ax[1].get_yaxis().set_visible(False)

    ax11 = ax[1].twinx()
    ax11.plot(pos_s, label='PCR-identified cases', color='red', alpha=0.7, lw=1.0)
    ax11.set_ylim((0, 200))
    ax11.tick_params(axis='y', colors='red')

    lines, labels = [], []
    for axis in [ax[0], ax01]: 
        axLine, axLabel = axis.get_legend_handles_labels()
        lines.extend(axLine)
        labels.extend(axLabel)
        
    fig.axes[1].legend(lines, labels, loc="upper right", fontsize=5, frameon=False)

    if fill: 
        x = np.arange(0, 200)
        ax[0].fill_between(x, inc[0:200], color='black', alpha=0.1)
        ax01.fill_between(x, pos_ns, color='red', alpha=0.1)
        ax[1].fill_between(x, inc[0:200], color='black', alpha=0.1)
        ax11.fill_between(x, pos_s, color='red', alpha=0.1)

    fig.tight_layout()

    if save: 
        fig.savefig('PCR_incidence_positives.jpg', dpi=fig.dpi, bbox_inches='tight')
    return


def plot_incidence_positives_RAT(inc, pos, fill=True, se=0.70, sp=0.99, f=3, p=0.5, save=True):
    
    mpl.rcParams['figure.dpi'] = 400

    fig, ax1 = plt.subplots(1, 1, figsize=(5, 3))

    ax1.plot(inc, label='Incidence of infection', c='black', alpha=0.7)
    ax1.plot(pos, label=f'RAT-identified cases\n(f={f}, p={p})', c='red',alpha=0.7)

    ax1.set_title('Active surveillance based on RAT (SE: ' + str(se) + ', SP: ' + str(sp) + ')', fontsize=11, y=1.03)
    ax1.set_ylabel('Number of individuals')
    ax1.set_xlabel('Days since beginning of outbreak')
    ax1.set_xlim((0, 200))
    ax1.set_ylim((0, 6000))
    ax1.set_xticks(np.arange(0, 201, 50))

    if fill: 
        x = np.arange(0, 200)
        ax1.fill_between(x, inc[0:200], color='black', alpha=0.1)
        ax1.fill_between(x, pos, color='red', alpha=0.1)

    lines, labels = ax1.get_legend_handles_labels()
    ax1.legend(lines, labels, loc='upper right', fontsize=7, labelspacing=0.8, frameon=False)  # Setting legend on axis instead of figure is better - avoids most manual tuning

    fig.tight_layout()

    if save: 
        plt.savefig('Antigen surveillance.jpg', dpi=1000)



## Plotting PCR allocation framework ##
def plot_meets_criteria(surv, 
                        scale='root', 
                        style=None, 
                        dpi=400, 
                        legend='on',
                        lw=2,
                        xlim=None,
                        ylim=None,
                        fig=None,
                        ax=None,
                        show_criteria=default_criteria,
                        colors=default_colors,
                        cmap=None,
                        fname=None): 
    '''
    Plot who meets each criteria

    Args: 
        surv            (surveillance): Surveillance system. Currently these plots are geared for PCR. 
        scale           (str)         : Scale y-axis linearly or according to sqrt function
        style           (str)         : Plotting style to use
        dpi             (int)         : Figure quality
        lw              (float)       : Linewidth 
        xlim            (tuple)       : x-axis limits, if provided 
        ylim            (tuple)       : y-axis limits, if provided
        show_criteria   (list)        : Testing criteria variable names 
        fname           (str)         : Location to save image. No saving if not provided. 
    '''
    mpl.rcParams['figure.dpi'] = dpi

    if surv.system == 'PassivePCR': 
        criteria_defs = PCR_criteria_defs
    else: 
        criteria_defs = RAT_criteria_defs
    
    if not(cmap is None): 
        colors = sc.vectocolor(len(surv.alloc_groups), cmap='jet')
    
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    if style == "science":  # Default to science
        plt.style.use(["science", "nature", "grid", "no-latex"])
    elif not(style is None): 
        raise RuntimeError("Invalid style provided. Please use 'science' or leave this argument empty to use default.")
    
    if (fig is None) and (ax is None): 
        fig, ax = plt.subplots(figsize=(6, 4))
        fig.suptitle('Number of individuals meeting each criterion',y=0.94, fontsize=13)

    # Set labels
    ax.set_xlabel('Days since beginning of outbreak', fontsize=11)
    ax.set_ylabel('Number of individuals', fontsize=11)
    ax.set_xlim((0, 200))

    # Set axis limits 
    if not(xlim is None):
        ax.set_xlim(xlim)
    if not(ylim is None): 
        ax.set_ylim(ylim)
    ax.set_xticks(np.arange(0, 201, 50))

    # Adjust y-scaling 
    if scale == 'root':
        ax.set_yscale('function', functions=(forward, inverse))
    elif scale == 'linear': 
        ax.set_yscale('linear')

    ax.tick_params(labelsize=11)

    # Get valid keys (set implementation does not retain order)
    valid = [key for key in show_criteria if key in surv.alloc_groups.keys()]

    # Plot 
    letter = 0
    for key, color in zip(valid, colors): 
        if key in surv.crit_groups.keys(): 
            ax.plot(surv.crit_groups[key], color=color, label=letters[letter] + ": " + criteria_defs[key], alpha=0.7, lw=lw)
            letter += 1 

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


def plot_criteria_seekers(surv,
                        scale='root',
                        style=None,
                        dpi=400,
                        legend='on',
                        lw=2,
                        xlim=None,
                        ylim=None,
                        fig=None,
                        ax=None,
                        show_criteria=default_criteria,
                        colors=default_colors,
                        cmap=None,
                        fname=None):  
    '''
    Plot test seekers for each criteria

    Args: 
        surv            (surveillance): Surveillance system. Currently these plots are geared for PCR. 
        scale           (str)         : Scale y-axis linearly or according to sqrt function
        style           (str)         : Plotting style to use
        dpi             (int)         : Figure quality
        lw              (float)       : Linewidth 
        xlim            (tuple)       : x-axis limits, if provided 
        ylim            (tuple)       : y-axis limits, if provided
        show_criteria   (list)        : Testing criteria variable names 
        fname           (str)         : Location to save image. No saving if not provided.
    '''
    mpl.rcParams['figure.dpi'] = dpi

    if surv.system == 'PassivePCR': 
        criteria_defs = PCR_criteria_defs
    else: 
        criteria_defs = RAT_criteria_defs

    if not(cmap is None): 
        colors = sc.vectocolor(len(surv.alloc_groups), cmap='jet')

    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    if style == "science":  # Default to science
        plt.style.use(["science", "nature", "grid", "no-latex"])
    elif not(style is None): 
        raise RuntimeError("Invalid style provided. Please use 'science' or leave this argument empty to use default.")

    if (fig is None) and (ax is None): 
        fig, ax = plt.subplots(figsize=(6, 4))
        fig.suptitle('Number of test seekers meeting each testing criterion',y=0.94, fontsize=13)

    # Set labels
    ax.set_xlabel('Days since beginning of outbreak', fontsize=11)
    ax.set_ylabel('Number of test seekers', fontsize=11)
    ax.set_xlim((0, 200))

    # Set axis limits 
    if not(xlim is None):
        ax.set_xlim(xlim)
    if not(ylim is None): 
        ax.set_ylim(ylim)
    ax.set_xticks(np.arange(0, 201, 50))

    # Adjust y-scaling 
    if scale == 'root':
        ax.set_yscale('function', functions=(forward, inverse))
    elif scale == 'linear': 
        ax.set_yscale('linear')

    ax.tick_params(labelsize=11)

    # Get valid keys (set implementation does not retain order)
    valid = [key for key in show_criteria if key in surv.alloc_groups.keys()]
    # set1 = set(show_criteria)  # User defined
    # set2 = set(surv.seek_groups.keys())  # All keys available
    # valid = set1.intersection(set2)

    # Plot 
    letter = 0
    for key, color in zip(valid, colors): 
        if key in surv.seek_groups.keys(): 
            ax.plot(surv.seek_groups[key], color=color, label=letters[letter] + ": " + criteria_defs[key], alpha=0.7, lw=lw)
        else: 
            print(key, "is not a valid criterion")
        letter += 1 

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


    
def plot_criteria_testers(surv,
                        scale='root',
                        style=None,
                        lw=2,
                        dpi=400,
                        legend='on',
                        xlim=None,
                        ylim=None,
                        fig=None,
                        ax=None,
                        show_criteria=default_criteria,
                        colors=default_colors,
                        cmap=None,
                        show_capacity=True,
                        show_consumed=True,
                        show_positives=True,
                        fname=None):  
    '''
    Plot testers for each criteria

    Args: 
        surv                (surveillance) : Surveillance system. Currently these plots are geared for PCR. 
        scale               (str)          : Scale y-axis linearly or according to sqrt function
        style               (str)          : Plotting style to use
        dpi                 (int)          : Figure quality
        lw                  (float)        : Linewidth 
        xlim                (tuple)        : x-axis limits, if provided 
        ylim                (tuple)        : y-axis limits, if provided
        show_criteria       (list)         : Testing criteria variable names
        show_capacity       (bool)         : Whether to plot test capacity 
        show_consumed       (bool)         : Whether to plot tests consumed
        show_positives      (bool)         : Whether to plot positive tests
        fname               (str)          : Location to save plot and file type. No saving if not provided.
    '''
    mpl.rcParams['figure.dpi'] = dpi

    if surv.system == 'PassivePCR': 
        criteria_defs = PCR_criteria_defs
    else: 
        criteria_defs = RAT_criteria_defs

    if not(cmap is None): 
        colors = sc.vectocolor(len(surv.alloc_groups), cmap='jet')

    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    if style == "science":  # Default to science
        plt.style.use(["science", "nature", "grid", "no-latex"])
    elif not(style is None): 
        raise RuntimeError("Invalid style provided. Please use 'science' or leave this argument empty to use default.")
    
    #TODO: Finish implementing this case where we want to plot figures in a multiplot
    if (fig is None) and (ax is None): 
        fig, ax = plt.subplots(figsize=(6, 4))
        fig.suptitle('Number of test recipients meeting each criterion',y=0.94, fontsize=13)


    # Set labels
    ax.set_xlabel('Days since beginning of outbreak', fontsize=11)
    ax.set_ylabel('Number of test recipients', fontsize=11)
    ax.set_xlim((0, 200))

    # Set axis limits 
    if not(xlim is None):
        ax.set_xlim(xlim)
    if not(ylim is None): 
        ax.set_ylim(ylim)
    ax.set_xticks(np.arange(0, 201, 50))

    # Adjust y-scaling 
    if scale == 'root':
        ax.set_yscale('function', functions=(forward, inverse))
    elif scale == 'linear': 
        ax.set_yscale('linear')

    ax.tick_params(labelsize=11)

    # Get valid keys (set implementation does not retain order)
    valid = [key for key in show_criteria if key in surv.alloc_groups.keys()]
    # set1 = set(show_criteria)  # User defined
    # set2 = set(surv.alloc_groups.keys())  # All keys available
    # valid = set1.intersection(set2)

    # Plot 
    letter = 0
    for key, color in zip(valid, colors): 
        if key in surv.alloc_groups.keys(): 
            ax.plot(surv.alloc_groups[key],color=color, label=letters[letter] + ": " + criteria_defs[key], alpha=0.7, lw=lw)
        else: 
            print(key, "is not a valid criterion")
        letter += 1 

    if show_capacity: 
        ax.plot(surv.capacity, label=letters[letter] + ': Test capacity', color='black', alpha=1.0, lw=lw, ls='dotted', dash_capstyle='round')
        letter += 1

    if show_consumed: 
        ax.plot(surv.tests_consumed, label=letters[letter] + ': Tests consumed', color='green', alpha=1.0, lw=lw)
        letter += 1
    
    if show_positives:
        ax.plot(surv.pos_history, label=letters[letter] + ': Positive PCR tests', color='red', alpha=1.0, lw=lw)
        letter += 1

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



def plot_meets_criteria_AB(dummy_ns, dummy_s, save=True, scale='linear', ylim=(0, 100000)): 
    plt.style.use("default")
    mpl.rcParams['figure.dpi'] = 400

    fig, ax = plt.subplots(1, 2, figsize=(9, 4))

    plot_meets_criteria(dummy_ns, 
                        scale=scale, 
                        legend='off',
                        ylim=ylim,
                        fig = fig,
                        ax = ax[0])

    plot_meets_criteria(dummy_s,  
                        scale=scale, 
                        legend='off',
                        ylim=ylim,
                        fig=fig,
                        ax=ax[1])

    lines_ns, labels_ns = fig.axes[0].get_legend_handles_labels()
    lines_s, labels_s = fig.axes[1].get_legend_handles_labels()

    lines_s[1].set_color(lines_ns[4].get_color())

    ax[1].legend(lines_ns, labels_ns, fontsize=5.5, loc='upper center', labelspacing=0.8, frameon=False)

    ax[0].set_xlabel('')
    ax[1].set_xlabel('')
    ax[1].set_ylabel('')
    ax[0].yaxis.labelpad=4
    ax[1].get_yaxis().set_ticklabels([])

    ax[0].set_title('A) Broad PCR Testing Eligibility', fontsize=10)
    ax[1].set_title('B) Stringent PCR Testing Eligibility', fontsize=10)

    ax[0].minorticks_off()
    ax[1].minorticks_off()

    fig.supxlabel('Days since beginning of outbreak', fontsize=11, y=0.03)
    fig.suptitle('Number of individuals meeting each testing criterion')

    fig.tight_layout()

    if save: 
        plt.savefig('eligibility.jpg', dpi=1000, bbox_inches='tight')
    


def plot_criteria_seekers_AB(dummy_ns, dummy_s, save=True, scale='linear', ylim=(0, 30000)): 
    plt.style.use('default')
    mpl.rcParams['figure.dpi'] = 400

    fig, ax = plt.subplots(1, 2, figsize=(9, 4))
    plt.minorticks_off()

    plot_criteria_seekers(dummy_ns, 
                                scale=scale, 
                                legend='off',
                                ylim=ylim,
                                fig = fig,
                                ax = ax[0])

    plot_criteria_seekers(dummy_s, 
                                scale=scale, 
                                legend='off',
                                ylim=ylim,
                                fig = fig,
                                ax = ax[1])

    lines_ns, labels_ns = fig.axes[0].get_legend_handles_labels()
    lines_s, labels_s = fig.axes[1].get_legend_handles_labels()

    lines_s[1].set_color(lines_ns[4].get_color())

    ax[1].legend(lines_ns, labels_ns, fontsize=5.5, loc='upper center', labelspacing=0.8, frameon=False)

    ax[0].set_xlabel('')
    ax[0].yaxis.labelpad=6

    ax[1].set_xlabel('')
    ax[1].set_ylabel('')
    ax[1].get_yaxis().set_ticklabels([])

    ax[0].set_title('A) Broad PCR Testing Eligibility', fontsize=10)
    ax[1].set_title('B) Stringent PCR Testing Eligibility', fontsize=10)

    ax[0].minorticks_off()
    ax[1].minorticks_off()

    fig.supxlabel('Days since beginning of outbreak', fontsize=11, y=0.03)
    fig.suptitle('Number of test seekers meeting each criterion')

    fig.tight_layout()

    if save: 
        plt.savefig('seeking.jpg', dpi=1000, bbox_inches='tight')


def plot_criteria_testers_AB(dummy_ns, 
                            dummy_s, 
                            save=True, 
                            scale='linear', 
                            ylim=(0, 800), 
                            show_capacity=True,
                            show_consumed=False,
                            show_positives=False): 

    plt.style.use("default")
    mpl.rcParams['figure.dpi'] = 400

    fig, ax = plt.subplots(1, 2, figsize=(9, 4))

    plot_criteria_testers(dummy_ns, 
                                scale=scale, 
                                ylim=ylim,
                                legend='off',
                                fig = fig,
                                ax = ax[0],
                                show_capacity=show_capacity, 
                                show_consumed=show_consumed,
                                show_positives=show_positives)

    plot_criteria_testers(dummy_s, 
                                scale=scale, 
                                ylim=ylim,
                                legend='off',
                                fig = fig,
                                ax = ax[1],
                                show_capacity=show_capacity, 
                                show_consumed=show_consumed,
                                show_positives=show_positives)

    lines_ns, labels_ns = fig.axes[0].get_legend_handles_labels()
    lines_s, labels_s = fig.axes[1].get_legend_handles_labels()

    lines_s[1].set_color(lines_ns[4].get_color())


    ax[1].legend(lines_ns, labels_ns, fontsize=5.5, loc='upper center', labelspacing=0.8, frameon=False)

    ax[0].set_xlabel('')
    ax[1].set_xlabel('')
    ax[1].set_ylabel('')
    ax[0].yaxis.labelpad=6
    ax[1].get_yaxis().set_ticklabels([])

    ax[0].set_title('A) Broad PCR Testing Eligibility', fontsize=10)
    ax[1].set_title('B) Stringent PCR Testing Eligibility', fontsize=10)

    ax[0].minorticks_off()
    ax[1].minorticks_off()

    fig.supxlabel('Days since beginning of outbreak', fontsize=11, y=0.03)
    fig.suptitle('Number of test recipients meeting each criterion')

    fig.tight_layout()
    
    if save: 
        plt.savefig('allocation.jpg', dpi=1000, bbox_inches='tight')


def plot_multi_fp_RAT(results, multisim=True, save=True, se=0.70, sp=0.99): 
    fig, ax = plt.subplots(1, 1, figsize=(5, 3))

    ax.set_title('Different frequencies (f) and participation (p) in population-wide RAT', fontsize=7)
    ax.set_xlabel('Days since beginning of outbreak')
    ax.set_ylabel('Number of positive cases')

    ax.set_xlim((0, 200))
    ax.set_ylim((0, 9000))
    ax.set_xticks(np.arange(0, 201, 50))

    if multisim: 
        for slist in results: 
            f = str(slist[0].frequency)
            p = str(int(slist[0].participation['all']*100))+'%'
            ax.plot(compute_pos_multisim_avg(slist), label='f: '+f+' days'+' '*(6-2*len(f))+'p: '+p, alpha=0.7)
    else: 
        for r in results: 
            f = str(r.frequency)
            p = str(int(r.participation['all']*100))+'%'
            ax.plot(r.pos_history, label='f: '+f+' days'+' '*(6-2*len(f))+'p: '+p, alpha=0.7)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.3)
    ax.text(5, 7800, 'SE: '+str(se*100)+'%\n'+'SP: '+str(sp*100)+'%', bbox=props)

    lines, labels = ax.get_legend_handles_labels()
    ax.legend(lines, labels, frameon=False)

    if save: 
        plt.savefig('antigen_manipulations.jpg', dpi=1000)

