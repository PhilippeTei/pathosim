
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from calibration.paths import BASE_PATH
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go

def plot_multiregion_results(sim):
    """
    Plots new diagnoses and severe cases for multiple regions. 
    """
    # parse
    regions = sim.rnames

    # Make a dataframe. 
    mr_results = pd.DataFrame()
    mr_results['date'] = sim.datevec

    for metric in ['new_diagnoses', 'new_severe']:
        for region in regions:
            col_name = region + '_' + metric
            mr_results[col_name] = sim.results[col_name].values
            
    # 2 subplots
    fig = make_subplots(rows=2, cols=1)


    for i, metric in enumerate(['new_diagnoses', 'new_severe']):
        # Supress legend for all but first subplot
        show_legend = True
        if i > 0:
            show_legend = False

        for j, region in enumerate(regions):
            
            col_name = region + '_' + metric
            fig.add_trace(go.Scatter
                                (x=mr_results['date'], y=mr_results[col_name], 
                                name=region, showlegend=show_legend, 
                                marker=dict(color=px.colors.qualitative.Plotly[j])
                                ), 
                            row=i+1, col=1
                            )
        fig.update_yaxes(title_text=metric, row=i+1, col=1)
        fig.update_xaxes(title_text='Date', row=i+1, col=1)

    fig.update_layout(title_text="Multi-Region Results", showlegend=True)
    import sys
    normal_stdout = sys.stdout
    sys.stdout = open('trash.txt', 'w')
    fig.show(renderer = 'png') # This function for some reason spits out a lot of text. Redirect it to a trash file. 
    sys.stdout = normal_stdout

def plot_eval_variant_prevelances(sim, variants, start_day, end_day, fpath):
    """
    1. Load the variant prevalences from the Sanger data.
    2. Extract the simulation data.
    3. Plot overlay. 
    """

    df_vars = pd.read_csv(fpath)
    # Convert to datetime
    df_vars['date'] = pd.to_datetime(df_vars['date'])

    sim_data = sim.results['variant']

    # Generate dates from start_date to end_date. 
    x_axis = pd.date_range(start=start_day, end=end_day)

    # Two figures. Top subplot for percentages. Bottom for counts. Counts are for sim only. 
    fig, ax = plt.subplots(2, 1, figsize=(10,10))
    # Top subplot
    ax[0].set_title('New infections by variant')
    ax[0].set_ylabel('Proportion of new infections')
    ax[0].set_xlabel('Date')

    for i, var in enumerate(variants):
        # Plot the normalized data. 
        totals = sim_data['new_infections_by_variant'].values.sum(axis=0)
        # Add a small number to avoid divide by zero.
        totals += 1e-10
        ax[0].plot(x_axis, sim_data['new_infections_by_variant'].values[i,:]/totals, label=var + ' sim')

        # Don't plot beta from the sanger data. 
        if var != 'beta':
            ax[0].plot(df_vars['date'], df_vars[var]/100, label=var + ' actual')

    ax[0].legend()

    # Bottom subplot
    ax[1].set_ylabel('Number of new infections')
    ax[1].set_xlabel('Date')
    for i, var in enumerate(variants):
        ax[1].plot(x_axis, sim_data['new_infections_by_variant'].values[i,:], label=var + 'sim')
        # plt.plot(df_vars['date'], df_vars[var], label=var)
    ax[1].legend()

    # Space out dates
    plt.xticks(rotation=45)

    plt.show()