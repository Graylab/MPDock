import argparse, os
import random
import numpy as np
import pandas as pd
from tqdm import tqdm


def get_scorefile_list(direc, filename):
    '''The function to obtain a list of score files for the PDB targets in a
    directory. Filename specifies the specific energy function to deal with.
    '''
    scorefile_list = []
    directory = os.path.abspath(direc)
    for pdbid in os.listdir(directory):
        if os.path.exists( os.path.join(directory, pdbid, filename) ):
            scorefile_list.append( os.path.join(directory, pdbid, filename) )
    
    return scorefile_list


def count_success(sample, n):
    '''
    This functions measures the metric to calculate success. It takes in the list and
    the randomly selected decoys that need to be chosen and returns the success count
    
    **Parameters**
        
        sample : *list,list,float,float,str*
        Nested list containing floats and strings .i.e. scores, rms and description
        of the PDB decoys
        
        n : *int*
        The integer of decoys that need to be randomly selected.
        
    **Returns**
        
        success_count : *int*
        The successful/acceptble decoys present in the random selection from 1000 decoys.
        
    '''
    success_count = 0

    for i in range(n):
        if sample[i][1] >= 1.0: # CAPRI Accepable or better
            success_count += 1

    return success_count


def estimate_Irms(sample, n):

    Irms = []
    for i in range(n):
        Irms.append(sample[i][2])
    
    avg_irms = np.mean(Irms)
    std_irms = np.std(Irms)

    return avg_irms, std_irms

def best_irms(sample):
    # I_sc, CAPRI_rank, Irms, fnat, desc
    sample.sort(key=lambda x: x[2])
    irms = []
    for i in range(5):
        irms.append(sample[i][2])
    avg_irms = np.mean(irms)
    best_irms = sample[0][2]
    return best_irms, avg_irms

def estimate_fnat(sample, n):

    fnat = []
    for i in range(n):
        fnat.append(sample[i][3])
    
    avg_fnat = np.mean(fnat)
    std_fnat = np.std(fnat)
    
    return avg_fnat, std_fnat

def estimate_fnat_irms(sample):
    fnat = []
    # sort so that you get the lowest RMSD metric
    sample.sort(key=lambda x: x[2])

    
    avg_fnat_n10 = np.mean(fnat[:10])
    # avg_fnat_n100 = np.std(fnat[:100])

    return avg_fnat_n10

def randomly_sample(scorefile, scoretype):
    '''
    The function assesses a population of 1000 randomly picked sample from the output,
    and return the N5, N10, N20, N100 and N1000 metric data.
    
    **Parameters**
        
        scorefile : *str*
        The filename of the scorefile.
        
    **Returns**
        
        n5,n10,n100,n200,n1000 : *int,int,int,int,int*
        The possibility of having a successful CAPRI acceptable solution within randomly 
        selected <N> decoys.
        
    '''

    #row_dict = { 'ref15' : 41, 'mp15' : 44, 'fa19' : 43, 'default' : 0 }    # for rigid 
    # row_dict = { 'ref15' : 47, 'mp15' : 50, 'fa19' : 49, 'default' : 0 }        # for ensemble
    row_dict = {  'mp15' : 50, 'fa19' : 39, 'fa23' : 41, 'mp15fa19': 49, 'mdsmp15': 41, 'mp15fa23': 43, 'default' : 0 }     # for updated ensemble (39 for updated, 47 for mp)
    row_len = row_dict[scoretype]

    #for one sample, to be repeated 1000 times
    sample = []

    for i in range(len(scorefile)-2):
        entry = random.randint(3, len(scorefile)-2)
        scores = scorefile[entry].split()
        
        # total_score, cen_rms, desc
        if len(scores) == row_len :
            if( scores[-1] == 'description'):
                continue
            else:
                # I_sc, CAPRI_rank, Irms, fnat, desc
                sample.append([float(scores[5]), float(scores[3]), float(scores[6]), float(scores[4]), scores[-1]])
        
    # sort for high res
    sample.sort(key=lambda x: x[0])

    return sample

def best_irms_all_sample(scorefile, scoretype):
    '''
    The function assesses a population of 1000 randomly picked sample from the output,
    and return the N5, N10, N20, N100 and N1000 metric data.
    
    **Parameters**
        
        scorefile : *str*
        The filename of the scorefile.
        
    **Returns**
        
        n5,n10,n100,n200,n1000 : *int,int,int,int,int*
        The possibility of having a successful CAPRI acceptable solution within randomly 
        selected <N> decoys.
        
    '''

    #row_dict = { 'ref15' : 41, 'mp15' : 44, 'fa19' : 43, 'default' : 0 }    # for rigid 
    # row_dict = { 'ref15' : 47, 'mp15' : 50, 'fa19' : 49, 'default' : 0 }        # for ensemble
    row_dict = {  'mp15' : 50, 'fa19' : 39, 'fa23' : 41, 'mp15fa19': 49, 'mdsmp15': 41, 'mp15fa23': 43, 'default' : 0 }     # for updated ensemble (39 for updated, 47 for mp)
    row_len = row_dict[scoretype]

    #for one sample, to be repeated 1000 times
    sample = []

    for i in range(len(scorefile)-2):
        
        scores = scorefile[i].split()
        
        # total_score, cen_rms, desc
        if len(scores) == row_len :
            if( scores[-1] == 'description'):
                continue
            else:
                # I_sc, CAPRI_rank, Irms, fnat, desc
                sample.append([float(scores[5]), float(scores[3]), float(scores[6]), float(scores[4]), scores[-1]])

    # sort for high res
    sample.sort(key=lambda x: x[2])
    return sample[0][2]

def bootstrap(source, filename, scoretype, destination):
    '''
    The function bootstraps all the files in the directory and generates
    necessary Enrichment data metrics which is essential to tabulate performance.
    
    **Parameters**
    
        source : *dir,str*
            The name of the source directory.
    
        destination : *dir,str*
            The name of the destination directory.
    
    '''

    scorefile_list = get_scorefile_list(source, filename)

    stats_file = open(destination, 'w')
    stats_file.write("PDBID\tAvg_N5\tStdev_N5\tAvg_N10\tStdev_N10\t"\
                     "Avg_N100\tStdev_N100\tAvg_Irms_N10\tStdev_Irms_N10\t"\
                     "Avg_Irms_N100\tStdev_Irms_N100\tAvg_fnat_N10\tStdev_fnat_N10\t"\
                     "Avg_fnat_N100\tStdev_fnat_N100\tE1\tE5\tBest_IRMS\tBest_IRMS_N1000\tAvg_IRMS_N5\tAvg_fnat_N10_irms\n")
    stats_file.flush()

    for scorefile_name in scorefile_list:

        print("Calculating for PDB: ", scorefile_name.split('/')[-2])
        if(scorefile_name.split('/')[-2] == '2qjy' and scoretype == 'fa19'):
            continue
        data = { 'N5' : [], 'N10' : [], 'N100' : [], 'Irms_N10' : [],
                 'std_Irms_N10' : [], 'std_Irms_N100' : [],
                 'Irms_N100' : [], 'fnat_N10' : [], 'fnat_N100' : [],
                 'std_fnat_N10' : [], 'std_fnat_N100' : [],
                 'E1' : [], 'E5' : [], 'best_irms' : [], 'best_irms_N1000' : [], 'avg_irms_N5': [],
                 'avg_fnat_N10_irms' : []}

        scorefile = open(scorefile_name, 'r').readlines()
        best_irms_all = best_irms_all_sample(scorefile, scoretype)
        # sample 1000 times
        for i in tqdm(range(1000)):
            sample = randomly_sample(scorefile, scoretype)
            data['N5'].append( count_success(sample, 5) )
            data['N10'].append( count_success(sample, 10) )
            data['N100'].append( count_success(sample, 100) )
            I_10, _ = estimate_Irms(sample, 10)
            I_100, _ =  estimate_Irms(sample, 100)
            f_10, _ = estimate_fnat(sample, 10)
            f_100, _ =  estimate_fnat(sample, 100)
            best_irms_1000, avg_irms_1000 = best_irms(sample)
            data['Irms_N10'].append( I_10 )
            data['Irms_N100'].append( I_100 )
            data['fnat_N10'].append( f_10 )
            data['fnat_N100'].append( f_100 )
            data['E1'].append( float(count_success(sample, 10) / 10) )
            data['E5'].append( float(count_success(sample, 50) / 50) )
            data['best_irms'].append(best_irms_all)
            data['best_irms_N1000'].append(best_irms_1000)
            data['avg_irms_N5'].append(avg_irms_1000)
            data['avg_fnat_N10_irms'].append(estimate_fnat_irms(sample))

    
        stats_file.write(scorefile_name.split('/')[-2] + "\t" + \
                     str(round(np.mean(data['N5']),3)) + "\t" + \
                     str(round(np.std(data['N5']),3)) + "\t" + str(round(np.mean(data['N10']),3)) + \
                     "\t" + str(round(np.std(data['N10']),3)) + "\t" + \
                     str(round(np.mean(data['N100']),3)) + "\t" + str(round(np.std(data['N100']),3)) +  "\t" + \
                     str(round(np.mean(data['Irms_N10']),3)) + "\t" + str(round(np.std(data['Irms_N10']),3)) +  "\t" + \
                     str(round(np.mean(data['Irms_N100']),3)) + "\t" + str(round(np.std(data['Irms_N100']),3)) +  "\t" + \
                     str(round(np.mean(data['fnat_N10']),3)) + "\t" + str(round(np.std(data['fnat_N10']),3)) +  "\t" + \
                     str(round(np.mean(data['fnat_N100']),3)) + "\t" + str(round(np.std(data['fnat_N100']),3)) +  "\t" + \
                     str(round(np.mean(data['E1']),3)) + "\t" + str(round(np.mean(data['E5']),3)) + "\t" +\
                     str(round(np.mean(data['best_irms']),3)) + "\t" + str(round(np.mean(data['best_irms_N1000']),3)) + "\t" + str(round(np.mean(data['avg_irms_N5']),3)) + "\n")
        stats_file.flush()
    
    stats_file.close()

def compare_irms():
    
    import plotly.graph_objects as go
    
    filepath_jabber = "/home/rsamant2/scratch16-jgray21/rsamant2/ensemble_docking/MPDockTest/scripts/Jabberdock.csv"
    df_jabber_data = pd.read_csv(filepath_jabber, delimiter=",")
    df_final = pd.DataFrame()
    for score_fxn in ['fa19','fa23','mp15']:
        file_in = '/home/rsamant2/scratch16-jgray21/rsamant2/ensemble_docking/Results/Ensembles/bootstrapping_results_'+score_fxn+'.csv'
        df = pd.read_csv(file_in, delimiter="\t")
        # print(df_extracted)
        
        if(df_final.empty):
            df_final = pd.merge(df_jabber_data, df, on="PDBID", suffixes=('_jd', '_'+score_fxn))
        else:
            df_final = pd.merge(df_final, df, on="PDBID", suffixes=('', '_'+score_fxn))
        
    df_final_medium = df_final[ df_final['Flexibility'] == 'Medium' ].sort_values( by=['Nat-Irms']).reset_index(drop=True)
    df_final_difficult = df_final[ df_final['Flexibility'] == 'Difficult' ].sort_values( by=['Nat-Irms']).reset_index(drop=True)
    
    fig = go.Figure()
    fig.add_trace( go.Bar(x=df_final_medium['PDBID'], y=df_final_medium['Best_IRMS_mp15'], error_y=dict(type='data', array=abs(df_final_medium['Avg_IRMS_N5_mp15']-df_final_medium['Best_IRMS_mp15'])), name="MP15", marker_color="gray") )
    fig.add_trace( go.Bar(x=df_final_medium['PDBID'], y=df_final_medium['Best_IRMS'], error_y=dict(type='data', array=abs(df_final_medium['Avg_IRMS_N5']-df_final_medium['Best_IRMS'])), name="Fa19", marker_color="skyblue") )
    fig.add_trace( go.Bar(x=df_final_medium['PDBID'], y=df_final_medium['Best_IRMS_fa23'], error_y=dict(type='data', array=abs(df_final_medium['Avg_IRMS_N5_fa23']-df_final_medium['Best_IRMS_fa23'])), name="Fa23", marker_color="teal") )
    fig.add_trace( go.Bar(x=df_final_medium['PDBID'], y=df_final_medium['Irms'], name="Jabberdock", marker_color="rgb(179,226,205)", marker_pattern_shape='x' ) )
 
    
    fig.update_xaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
    fig.update_yaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
    fig.update_xaxes(mirror=True, linewidth=2, linecolor='black',
                        tickfont=dict(family="Arial", color='black',size=24),ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
    fig.update_yaxes(title_text = 'IRMS', range=[0,10],tickfont=dict(family="Arial", color='black',size=24),
                        mirror=True, linewidth=2, linecolor='black',showgrid=True,ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
    fig.update_traces(marker_line_color='black',marker_line_width=0.7)
    fig.update_layout(showlegend=True)
    fig.update_layout(autosize=True, height=450, width=50,legend=dict(orientation="h",yanchor="bottom",y=1.02,xanchor="left",x=0,font=dict(size=18),itemwidth=30))

    fig.update(layout_coloraxis_showscale=True)
    fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
    fig.write_image( "comparing_irms_allscore_medium.png", format="png", width=600, height=450, scale=2)
    # fig.show()
    
    fig = go.Figure()
    fig.add_trace( go.Bar(x=df_final_difficult['PDBID'], y=df_final_difficult['Best_IRMS_mp15'], error_y=dict(type='data', array=abs(df_final_difficult['Avg_IRMS_N5_mp15']-df_final_difficult['Best_IRMS_mp15'])), name="MP15", marker_color="gray") )
    fig.add_trace( go.Bar(x=df_final_difficult['PDBID'], y=df_final_difficult['Best_IRMS'], error_y=dict(type='data', array=abs(df_final_difficult['Avg_IRMS_N5']-df_final_difficult['Best_IRMS'])), name="Fa19", marker_color="skyblue") )
    fig.add_trace( go.Bar(x=df_final_difficult['PDBID'], y=df_final_difficult['Best_IRMS_fa23'], error_y=dict(type='data', array=abs(df_final_difficult['Avg_IRMS_N5_fa23']-df_final_difficult['Best_IRMS_fa23'])), name="Fa23", marker_color="teal") )
    fig.add_trace( go.Bar(x=df_final_difficult['PDBID'], y=df_final_difficult['Irms'], name="Jabberdock", marker_color="rgb(179,226,205)", marker_pattern_shape='x' ) )
 
    
    fig.update_xaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
    fig.update_yaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
    fig.update_xaxes(mirror=True, linewidth=2, linecolor='black',
                        tickfont=dict(family="Arial", color='black',size=24), ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
    fig.update_yaxes(title_text = 'IRMS', range=[0,10],tickfont=dict(family="Arial", color='black',size=24),
                        mirror=True, linewidth=2, linecolor='black',showgrid=True,ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
    fig.update_traces(marker_line_color='black',marker_line_width=0.7)
    fig.update_layout(showlegend=True)
    fig.update_layout(autosize=True, height=450, width=50,legend=dict(orientation="h",yanchor="bottom",y=1.02,xanchor="left",x=0,font=dict(size=18),itemwidth=30))

    fig.update(layout_coloraxis_showscale=True)
    fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
    fig.write_image( "comparing_irms_allscore_difficult.png", format="png", width=600, height=450, scale=2)
    # fig.show()
    
def compare_irms_score():
    
    import plotly.graph_objects as go
    
    filepath_jabber = "/home/rsamant2/scratch16-jgray21/rsamant2/ensemble_docking/MPDockTest/scripts/Jabberdock.csv"
    df_jabber_data = pd.read_csv(filepath_jabber, delimiter=",")
    df_final = pd.DataFrame()
    for score_fxn in ['fa19','fa23','mp15']:
        file_in = '/home/rsamant2/scratch16-jgray21/rsamant2/ensemble_docking/Results/Ensembles/bootstrapping_results_'+score_fxn+'.csv'
        df = pd.read_csv(file_in, delimiter="\t")
        # print(df_extracted)
        
        if(df_final.empty):
            df_final = pd.merge(df_jabber_data, df, on="PDBID", suffixes=('_jd', '_'+score_fxn))
        else:
            df_final = pd.merge(df_final, df, on="PDBID", suffixes=('', '_'+score_fxn))
        
    df_final_medium = df_final[ df_final['Flexibility'] == 'Medium' ].sort_values( by=['Nat-Irms']).reset_index(drop=True)
    df_final_difficult = df_final[ df_final['Flexibility'] == 'Difficult' ].sort_values( by=['Nat-Irms']).reset_index(drop=True)
    
    fig = go.Figure()
    fig.add_trace( go.Bar(x=df_final_medium['PDBID'], y=df_final_medium['Avg_Irms_N10_mp15'], name="MP15", marker_color="gray") )
    fig.add_trace( go.Bar(x=df_final_medium['PDBID'], y=df_final_medium['Avg_Irms_N10'], name="Fa19", marker_color="skyblue") )
    fig.add_trace( go.Bar(x=df_final_medium['PDBID'], y=df_final_medium['Avg_Irms_N10_fa23'], name="Fa23", marker_color="teal") )
    fig.add_trace( go.Bar(x=df_final_medium['PDBID'], y=df_final_medium['Irms'], name="Jabberdock", marker_color="rgb(179,226,205)", marker_pattern_shape='x' ) )
 
    
    fig.update_xaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
    fig.update_yaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
    fig.update_xaxes(mirror=True, linewidth=2, linecolor='black',
                        tickfont=dict(family="Arial", color='black',size=24),ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
    fig.update_yaxes(title_text = 'IRMS', range=[0,10],tickfont=dict(family="Arial", color='black',size=24),
                        mirror=True, linewidth=2, linecolor='black',showgrid=True,ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
    fig.update_traces(marker_line_color='black',marker_line_width=0.7)
    fig.update_layout(showlegend=True)
    fig.update_layout(autosize=True, height=450, width=50,legend=dict(orientation="h",yanchor="bottom",y=1.02,xanchor="left",x=0,font=dict(size=18),itemwidth=30))

    fig.update(layout_coloraxis_showscale=True)
    fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
    fig.write_image( "comparing_irms_scoresort_medium.png", format="png", width=600, height=450, scale=2)
    # fig.show()
    
    fig = go.Figure()
    fig.add_trace( go.Bar(x=df_final_difficult['PDBID'], y=df_final_difficult['Avg_Irms_N10_mp15'], name="MP15", marker_color="gray") )
    fig.add_trace( go.Bar(x=df_final_difficult['PDBID'], y=df_final_difficult['Avg_Irms_N10'], name="Fa19", marker_color="skyblue") )
    fig.add_trace( go.Bar(x=df_final_difficult['PDBID'], y=df_final_difficult['Avg_Irms_N10_fa23'], name="Fa23", marker_color="teal") )
    fig.add_trace( go.Bar(x=df_final_difficult['PDBID'], y=df_final_difficult['Irms'], name="Jabberdock", marker_color="rgb(179,226,205)", marker_pattern_shape='x' ) )
 
    
    fig.update_xaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
    fig.update_yaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
    fig.update_xaxes(mirror=True, linewidth=2, linecolor='black',
                        tickfont=dict(family="Arial", color='black',size=24), ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
    fig.update_yaxes(title_text = 'IRMS', range=[0,10],tickfont=dict(family="Arial", color='black',size=24),
                        mirror=True, linewidth=2, linecolor='black',showgrid=True,ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
    fig.update_traces(marker_line_color='black',marker_line_width=0.7)
    fig.update_layout(showlegend=True)
    fig.update_layout(autosize=True, height=450, width=50,legend=dict(orientation="h",yanchor="bottom",y=1.02,xanchor="left",x=0,font=dict(size=18),itemwidth=30))

    fig.update(layout_coloraxis_showscale=True)
    fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
    fig.write_image( "comparing_irms_scoresort_difficult.png", format="png", width=600, height=450, scale=2)
    # fig.show()
    
    
    
def plotly_scatter(destination, compare_scorefile, scoretype, compare_scoretype):
    import plotly.graph_objects as go
    
    filepath_jabber = "/home/rsamant2/scratch16-jgray21/rsamant2/ensemble_docking/MPDockTest/scripts/Jabberdock.csv"
    df_jabber_data = pd.read_csv(filepath_jabber, delimiter=",")
    df_data = pd.read_csv(destination, delimiter="\t")
    
    df_comparedata = pd.read_csv(compare_scorefile, delimiter="\t")
    df_comparedata.columns = [str(col) + '_' + compare_scoretype if col!='PDBID' else str(col) for col in df_comparedata.columns ]
    print(df_data)
    print(df_comparedata)
    
    df_final = pd.merge(pd.merge(df_data,df_comparedata, on='PDBID'),df_jabber_data, on='PDBID')
    print(df_final)
    
    rigid_data = df_final[df_final['Flexibility']=='Rigid']
    medium_data = df_final[df_final['Flexibility']=='Medium']
    difficult_data = df_final[df_final['Flexibility']=='Difficult']
    
    
    fig = go.Figure()

    max_lim = float(max(max(df_final['Avg_N5']),max(df_final['Avg_N5_'+compare_scoretype])))+0.2
    min_lim = float(min(min(df_final['Avg_N5']),min(df_final['Avg_N5_'+compare_scoretype])))-0.2
    
    x = np.linspace(min_lim, max_lim, num=5)
    x_1 = np.linspace(min_lim, max_lim-1, num=5)
    x_m1 = np.linspace(min_lim+1, max_lim, num=5)
    
    y = x
    y_1 = x_1+1
    y_m1 = x_m1-1 
    
    # fig.add_trace(go.Scatter(x=df_final['Avg_N5_'+compare_scoretype], y=df_final['Avg_N5'], text=list(df_final['PDBID'][:]),                   
    #                     mode='markers+text',marker_symbol='circle',marker_color="rgb(240,96,96)",marker_line_color="rgb(191,191,191)", 
    #                     marker_line_width=1.0,marker_size=10,textposition='top left',textfont=dict(family="Helvetica",size=15,color="rgb(240,96,96)"),
    #                     showlegend=False))  
    fig.add_trace(go.Scatter(x=rigid_data['Avg_N5_'+compare_scoretype], y=rigid_data['Avg_N5'], text=list(df_final['PDBID'][:]),                   
                        mode='markers',marker_symbol='circle',marker_color="rgb(128,128,128)",marker_line_color="rgb(0,0,0)", 
                        marker_line_width=1.0,marker_size=10,textposition='top left',textfont=dict(family="Helvetica",size=15,color="rgb(240,96,96)"),
                        showlegend=True, name='Rigid'))
    fig.add_trace(go.Scatter(x=medium_data['Avg_N5_'+compare_scoretype], y=medium_data['Avg_N5'], text=list(df_final['PDBID'][:]),                   
                        mode='markers',marker_symbol='triangle-up',marker_color="rgb(0,0,205)",marker_line_color="rgb(0,0,0)", 
                        marker_line_width=1.0,marker_size=10,textposition='top left',textfont=dict(family="Helvetica",size=15,color="rgb(240,96,96)"),
                        showlegend=True, name='Medium'))  
    fig.add_trace(go.Scatter(x=difficult_data['Avg_N5_'+compare_scoretype], y=difficult_data['Avg_N5'], text=list(df_final['PDBID'][:]),                   
                        mode='markers',marker_symbol='diamond',marker_color="rgb(240,96,96)",marker_line_color="rgb(0,0,0)", 
                        marker_line_width=1.0,marker_size=10,textposition='top left',textfont=dict(family="Helvetica",size=15,color="rgb(240,96,96)"),
                        showlegend=True, name='Flexible')) 
     
    fig.add_trace(go.Scatter(x=x, y=y,mode='lines',
                        line = dict(color='black', width=1.5), showlegend=False))
    fig.add_trace(go.Scatter(x=x_1, y=y_1,mode='lines',
                        line = dict(color='gray', width=1.5), showlegend=False))
    fig.add_trace(go.Scatter(x=x_m1, y=y_m1,mode='lines',
                        line = dict(color='gray', width=1.5), showlegend=False))
    
    fig.update_xaxes(title_text = compare_scoretype + ' <N5>', range=[min_lim, max_lim], mirror=True, linewidth=2, linecolor='black',
                     tickfont=dict(family="Arial", color='black',size=20), ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
        
    fig.update_yaxes(title_text = scoretype + ' <N5>',range=[min_lim, max_lim],tickfont=dict(family="Arial", color='black',size=20),
                    mirror=True, linewidth=2, linecolor='black',showgrid=True, ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
    
    fig.update_layout(font=dict(family="Arial", color='black',size=20))
    fig.update_layout(xaxis=dict(dtick=1.0, tickvals = [0,1,2,3,4,5]))
    fig.update_layout(yaxis=dict(dtick=1.0,tickvals = [0,1,2,3,4,5]))
    
    fig.update_xaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
    fig.update_yaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
    fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
    fig.update_layout(autosize=True, height=450, width=50,legend=dict(orientation="h",yanchor="bottom",y=1.02,xanchor="left",x=0,font=dict(size=18),itemwidth=30))

    filename = 'comparing_fnat_'+scoretype+'_'+compare_scoretype+'.png'
    fig.write_image( os.path.join('./',filename ), format="png", width=550, height=550, scale=2)
    
#########  COMMAND LINE COMPATIBILITY   ###############################
def compare_all_score():
    
    import seaborn as sns
    import matplotlib
    import matplotlib.pyplot as plt
    
    
    pdbid=['3kly','3chx','1zoy','1e12','1bl8']
    df_final = pd.DataFrame()
    for score_fxn in ['fa19','fa23','mdsmp15','mp15','mp15fa19']:
        file_in = '/home/rsamant2/scratch16-jgray21/rsamant2/ensemble_docking/Results/Ensembles/bootstrapping_results_'+score_fxn+'.csv'
        df = pd.read_csv(file_in, delimiter="\t")
        df_extracted = df[df['PDBID'].isin(pdbid)]
        # print(df_extracted)
        df_extracted['method'] = score_fxn
        # print(df_extracted)
        df_final = pd.concat([df_final,df_extracted])
    print(df_final)    
    
    theme = {'axes.grid': True,
        'grid.linestyle': '',
        'xtick.labelsize': 18,
        'ytick.labelsize': 18,
        "font.weight": 'regular',
        "font.family": 'Arial',
        'xtick.color': 'black',
        'ytick.color': 'black',
        "axes.titlesize": 20,
        "axes.labelsize": 18,
        'axes.linewidth':2
    }

    matplotlib.rcParams.update(theme)
    my_pal = {'mdsmp15':'#333a89','fa19':'#636bc5','fa23':'#b1b5e2','mp15':'#113015','mp15fa19':'#21612a'}#,'mp15fa23':'#719075'}
    ax = sns.boxplot(x = df_final['method'],
                y = df_final['Avg_N5'],
                palette = my_pal,
                order=['mdsmp15','fa19','fa23','mp15','mp15fa19'])#,'mp15fa23'])
    
            
    ax.set_xlabel('Score Types')
    ax.set_ylabel('<N5>')
    ax.set_xticklabels(['Mp15','Fa19','Fa23','Mp15','Fa19'])#,'Fa23'])
    ax.xaxis.set_tick_params(width=2)
    ax.yaxis.set_tick_params(width=2)
    plt.tight_layout(rect=[0.0, 0.07, 1, 1.0])
    plt.xticks(rotation=90)
    # plt.yticks(np.arange(-4.0, 4.0, 0.5))
    outfile = 'comparing_fnat_allscore.png'
    plt.savefig(outfile, dpi=600, transparent=True)
    plt.close()
        

def main():

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-s', '--src', help='Directory with the score files', required = True)
    parser.add_argument('-d', '--destination', help='Directory with the bound score files', required = True)
    parser.add_argument('-f', '--scorefile', help = 'number of decoys to consider', required=True)
    parser.add_argument('-t', '--scoretype', help = 'directory to write plot files', required = True)
    parser.add_argument('-c', '--compare_scorefile', help = 'Directory with the bound compare-score files', default='./Ensembles/bootstrapping_results_fa19.csv')
    parser.add_argument('-p', '--compare_scoretype', help = 'Compare scoretype, i.e. fa19, mp15, fa23, etc', default = 'fa19')
    
    
    args = parser.parse_args()
    src = args.src
    destination = args.destination
    scorefile = args.scorefile
    scoretype = args.scoretype

    # bootstrap( src, scorefile, scoretype, destination )
    # compare_scorefile = args.compare_scorefile
    # compare_scoretype = args.compare_scoretype

    # if(compare_scoretype != None):
    #     plotly_scatter(destination, compare_scorefile, scoretype, compare_scoretype)
        
    compare_all_score()
    # compare_irms()
    # compare_irms_score()
if __name__ == "__main__":
    main()







