#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import argparse, os, sys, math, optparse
import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go


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


def get_native_list(direc, filename):
    '''The function to obtain a list of the native files of the PDB targets
    from a directory with defined target type
    '''
    scorefile_list = []
    directory = os.path.abspath(direc)
    for pdbid in os.listdir(directory):
        if os.path.exists( os.path.join(directory, pdbid, filename) ):
            scorefile_list.append( os.path.join(directory, pdbid, filename) )
    
    return scorefile_list


def parse_score_file(filename, scoretype, n):
    '''This function parses the high resolution score file to obtain necessary data
    '''

    # row_dict = {  'mp15' : 44, 'fa19' : 43, 'fa23': 46, 'default' : 0 }    # for rigid 
    #rigid
    #row_dict = {  'mp15' : 44, 'fa19' : 43, 'fa23':41, 'default' : 0 }        # for ensemble
    row_dict = {  'mp15' : 50, 'fa19' : 39, 'fa23': 41, 'default' : 0 }     # for updated ensemble (39 for updated, 47 for mp)
    row_len = row_dict[scoretype]

    print("Parsing file : ", filename)
    scorefile = open(filename, 'r').readlines()
    incorrect = []
    acceptable = []
    medium = []
    high = []
    low_cen_sc = 0
    high_cen_sc = -float("inf")
    
    if n == -1:
        n = len(scorefile) - 2
    elif n > len(scorefile) - 2:
        print("The score file only has " + str (len(scorefile) - 2) + "entries.")
        n = len(scorefile) - 2 
    
    for i in range(1, n+2):
        score_split = scorefile[i].split()

        if len(score_split) == row_len and score_split[-1] != 'description':

            if score_split[-1] != 'description':
                
                Isc = float(score_split[5])
            
                # To determine the Y-axes limits later on
                if Isc < low_cen_sc:
                    low_cen_sc = Isc
            
                if Isc > high_cen_sc:
                    high_cen_sc = Isc
                
                score_terms = [ score_split[-1], float(score_split[6]), Isc ]

                capri_rank = float(score_split[3])
                if( capri_rank.is_integer() ):
                    if score_split[3] == "0.000":
                        incorrect.append( score_terms )
                    elif score_split[3] == "1.000":
                        acceptable.append( score_terms )
                    elif score_split[3] == "2.000":
                        medium.append( score_terms )
                    elif score_split[3] == "3.000":
                        high.append( score_terms )
                    else:
                        sys.exit("Unreachable CAPRI rank")
                else:
                    print("CAPRI rank is not an integer")
            
            else:
                if( score_split[-1] == 'description' ):
                    continue
                else:
                    print( "WARNING: There is an issue with line number: " + str(i+1) )
            
    incorrect = list(map( list, zip(*incorrect) ))
    if len(acceptable) > 0:
        acceptable = list( map( list, zip(*acceptable) ))
    if len(medium) > 0:
        medium = list(map( list, zip(*medium) ))
    if len(high) > 0:
        high = list(map( list, zip(*high) ))
    
    scores = [incorrect, acceptable, medium, high, low_cen_sc, high_cen_sc]
    
    return scores


def parse_native_file(filename, scoretype, n=100):
    '''This function parses the score file to obtain necessary data
    '''
    
    print("Parsing file : ", filename)
    scorefile = open(filename, 'r').readlines()
    bound_scores = []
    low_cen_sc = 0

    row_dict = { 'mp15' : 50, 'fa19' : 39, 'fa23': 41, 'default' : 0 }
    row_len = row_dict[scoretype]

    if n == -1:
        n = len(scorefile) - 2
    elif n > len(scorefile) - 2:
        print("The score file only has " + str (len(scorefile) - 2) + "entries.")
        n = len(scorefile) - 2 

    for i in range(2, len(scorefile)):
        score_split = scorefile[i].split()

        #if len(score_split) == row_len and score_split[-1] != 'description':
        if score_split[-1] != 'description':
            if float(score_split[5]) < low_cen_sc :
                low_cen_sc = float(score_split[5])
            
            # Each list has decoy_name, Irms, I_sc
            Irms = float(score_split[6])
            Isc = float(score_split[5])
            if( Irms < 2.6 ) and (Isc < low_cen_sc+15 ) :
                bound_scores.append( [score_split[-1], Irms, Isc ] )
        
        else:
             print( "WARNING: There is an issue with line number: " + str(i+1) )
    
    bound_scores = list( map( list, zip(*bound_scores) ))
    scores = [bound_scores, low_cen_sc]
    
    return scores


def get_scorefile_quadruplets(mp_list, fa_list, fa23_list, native_mp_list, native_fa_list, native_fa23_list):
    
    if len(mp_list) != len(fa_list) or len(mp_list) != len(native_mp_list) or len(mp_list) != len(native_fa_list) :
        if len(fa_list) != len(fa23_list) or len(fa_list) != len(native_fa_list) or len(fa23_list) != len(native_fa23_list) :
            print( "\n"+"Lists are not the same: mp15 : " + str( len(mp_list) ) 
                + " v/s fa19 " + str( len(fa_list) ) + " v/s fa23 " + str( len(fa23_list) ) +" v/s Native : " + str( len(native_mp_list) ) + "\n" )

    
    file_quadruplets = []
    
    
    for mp_filename in mp_list:
        pdb = mp_filename.split('/')[-2]
        partners_found = False
        # print(pdb)
        for fa_filename in fa_list:
            if fa_filename.split('/')[-2] == pdb:
                # print(fa_filename.split('/')[-2])
                for fa23_filename in fa23_list:
                    if fa23_filename.split('/')[-2] == pdb:
                        
                        for native_mp in native_mp_list:
                            if native_mp.split('/')[-2] == pdb:

                                for native_fa in native_fa_list:
                                    if native_fa.split('/')[-2] == pdb:
                        
                                        for native_fa23 in native_fa23_list:
                                            if native_fa23.split('/')[-2] == pdb:
                                                # print(pdb)
                                                file_quadruplets.append( [ mp_filename, fa_filename, fa23_filename, native_mp, native_fa, native_fa23] )
                                                partners_found = True
                                                break
                                                
                                    if partners_found:
                                        break

                            if partners_found:
                                break

                    if partners_found:
                        break
            if partners_found:
                break
        
        if not partners_found:
            print("\n" + "Partners not found for '" + str(mp_filename) + "\n")

    return file_quadruplets
    

def make_plots(pdbid, mp_scores, fa_scores, fa23_scores, bound_mp, bound_fa, bound_fa23, plot_dir):
    
    y1_min = min( float(mp_scores[4]), float(bound_mp[1]) ) - 0.5
    y2_min = min( float(fa_scores[4]), float(bound_fa[1]) ) - 0.5
    y3_min = min( float(fa23_scores[4]), float(bound_fa23[1]) ) - 0.5
    y_min = min( y1_min, y2_min, y3_min)
    
    # fig = make_subplots( rows=1, cols=1, shared_yaxes=True, horizontal_spacing=0.1, vertical_spacing=0.1 )
    fig = go.Figure()
    # if( len(mp_scores[0]) > 0 ):
    #     fig.add_trace( go.Scatter(x = mp_scores[0][1], y = mp_scores[0][2], mode="markers", marker_color="rgb(191,191,191)", marker_line_width=0.5 ), row=1, col=1 )
    # if( len(mp_scores[1]) > 0 ):
    #     fig.add_trace( go.Scatter(x = mp_scores[1][1], y = mp_scores[1][2], mode="markers", marker_color="rgb(255,204,51)",  marker_line_width=0.5 ), row=1, col=1 )
    # if( len(mp_scores[2]) > 0 ):
    #     fig.add_trace( go.Scatter(x = mp_scores[2][1], y = mp_scores[2][2], mode="markers", marker_color="rgb(242,64,0)", marker_line_width=0.5 ), row=1, col=1 )
    # if( len(mp_scores[3]) > 0 ):
    #     fig.add_trace( go.Scatter(x = mp_scores[3][1], y = mp_scores[3][2], mode="markers", marker_color="rgb(0,204,0)", marker_line_width=0.5 ), row=1, col=1 )
    # if( len(bound_mp[0]) > 0 ):
    #     fig.add_trace( go.Scatter(x = bound_mp[0][1], y = bound_mp[0][2], mode="markers", marker_color="rgb(0,0,191)", marker_symbol=17, marker_line_width=0.3 ),
    #                   row=1, col=1 )
    print(fa_scores)
    
    if( len(fa_scores[0]) > 0 ):
        fig.add_trace( go.Scatter(x = fa_scores[0][1], y = fa_scores[0][2], mode="markers", marker_color="rgb(191,191,191)", marker_line_width=0.5 ) )
    if( len(fa_scores[1]) > 0 ):
        fig.add_trace( go.Scatter(x = fa_scores[1][1], y = fa_scores[1][2], mode="markers", marker_color="rgb(255,204,51)",  marker_line_width=0.5 ) )
    if( len(fa_scores[2]) > 0 ):
        fig.add_trace( go.Scatter(x = fa_scores[2][1], y = fa_scores[2][2], mode="markers", marker_color="rgb(242,64,0)", marker_line_width=0.5 ) )
    if( len(fa_scores[3]) > 0 ):
        fig.add_trace( go.Scatter(x = fa_scores[3][1], y = fa_scores[3][2], mode="markers", marker_color="rgb(0,204,0)", marker_line_width=0.5 ) )
    if( len(bound_fa[0]) > 0 ):
        fig.add_trace( go.Scatter(x = bound_fa[0][1], y = bound_fa[0][2], mode="markers", marker_color="rgb(0,0,191)", marker_symbol=17, marker_line_width=0.3 ) )

    # if( len(fa23_scores[0]) > 0 ):
    #     fig.add_trace( go.Scatter(x = fa23_scores[0][1], y = fa23_scores[0][2], mode="markers", marker_color="rgb(191,191,191)", marker_line_width=0.5 ), row=1, col=3 )
    # if( len(fa23_scores[1]) > 0 ):
    #     fig.add_trace( go.Scatter(x = fa23_scores[1][1], y = fa23_scores[1][2], mode="markers", marker_color="rgb(255,204,51)",  marker_line_width=0.5 ), row=1, col=3 )
    # if( len(fa23_scores[2]) > 0 ):
    #     fig.add_trace( go.Scatter(x = fa23_scores[2][1], y = fa23_scores[2][2], mode="markers", marker_color="rgb(242,64,0)", marker_line_width=0.5 ), row=1, col=3 )
    # if( len(fa23_scores[3]) > 0 ):
    #     fig.add_trace( go.Scatter(x = fa23_scores[3][1], y = fa23_scores[3][2], mode="markers", marker_color="rgb(0,204,0)", marker_line_width=0.5 ), row=1, col=3 )
    # if( len(bound_fa23[0]) > 0 ):
    #     fig.add_trace( go.Scatter(x = bound_fa23[0][1], y = bound_fa23[0][2], mode="markers", marker_color="rgb(0,0,191)", marker_symbol=17, marker_line_width=0.3 ),
    #                   row=1, col=3 )
        
    y1_max = min(0, float(mp_scores[5])) - 0.5
    y2_max = min(0, float(fa_scores[5])) - 0.5
    y3_max = min(0, float(fa23_scores[5])) - 0.5
    y_max = min( y1_max, y2_max, y3_max)
    
    fig.update_xaxes(range=[0, 15], mirror=True, linewidth=2, linecolor='black',
                     tickfont=dict(family="Arial", color='black',size=18))
    # fig.update_xaxes(range=[0, 15], mirror=True, linewidth=2, linecolor='black',
    #                  tickfont=dict(family="Arial", color='black',size=18), row=1, col=2)  
    # fig.update_xaxes(range=[0, 15], mirror=True, linewidth=2, linecolor='black',
    #                  tickfont=dict(family="Arial", color='black',size=18), row=1, col=3)              

    fig.update_yaxes(range=[y_min, y_max],tickfont=dict(family="Arial", color='black',size=18),
                     mirror=True, linewidth=2, linecolor='black',showgrid=True)
    # fig.update_yaxes(range=[y_min, y_max],tickfont=dict(family="Arial", color='black',size=18),
    #                  mirror=True, linewidth=2, linecolor='black',showgrid=True, row=1, col=2)
    # fig.update_yaxes(range=[y_min, y_max],tickfont=dict(family="Arial", color='black',size=18),
    #                 mirror=True, linewidth=2, linecolor='black',showgrid=True, row=1, col=3)


    fig.update_layout( xaxis=dict(tick0 = 0, dtick=3))#, xaxis2= dict(tick0 = 0, dtick=3), xaxis3= dict(tick0 = 0, dtick=3) )
    fig.update_xaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
    fig.update_yaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
    fig.update_yaxes(tickmode="auto", nticks=7)
    fig.update_traces(marker=dict(size=5,line=dict(width=0.4,color='Black')))
    fig.update_layout(showlegend=False)
    fig.update(layout_coloraxis_showscale=False)
    fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
    # fig.update_layout( yaxis3=dict( title = pdbid.upper(), titlefont=dict( family="Arial", size=24), side="right" ) )
    
    print("write_the_image")
    fig.write_image( os.path.join(plot_dir, pdbid +"_ensemble_ppy_v2.png"), format="png", width=350, height=350, scale=2)



def main():
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-s', '--src', help='Directory with the score files', required = True)
    parser.add_argument('-b', '--native', help='Directory with the bound score files', required = True)
    parser.add_argument('-n', '--num_decoys', help = 'number of decoys to consider', type=int, default=-1)
    parser.add_argument('-d', '--plot_dir', help = 'directory to write plot files', required = True)
    parser.add_argument('-p', '--pdbid', help = 'plot for particular pdb ID')
    
    args = parser.parse_args()
    src = args.src
    native = args.native
    n = args.num_decoys
    plot_dir = args.plot_dir
    pdb = args.pdbid
    print(pdb)
    mp_list = get_scorefile_list(src, 'docking_mp15.sc')
    fa_list = get_scorefile_list(src, 'docking_fa19.sc')
    fa23_list = get_scorefile_list(src, 'updated_docking_fa23.sc')
    # print(mp_list)
    # print(fa_list)
    # print(fa_list)
    # print(fa23_list)
    # break

    native_mp_list = get_native_list(native, 'native_mp15.sc')
    native_fa_list = get_native_list(native, 'native_fa19.sc')
    native_fa23_list = get_native_list(native, 'native_fa23.sc')
    # print(native_fa_list)
    # print(native_fa23_list)
    
    scorefile_quadruplets = get_scorefile_quadruplets( mp_list, fa_list, fa23_list, native_mp_list, native_fa_list, native_fa23_list)
    # print(scorefile_quadruplets)
    
    for scorefiles in scorefile_quadruplets:
        mp_scores = parse_score_file(scorefiles[0], 'mp15', n)
        fa_scores = parse_score_file(scorefiles[1], 'fa19', n)
        fa23_scores = parse_score_file(scorefiles[2], 'fa23', n)

        bound_mp = parse_native_file(scorefiles[3], 'mp15')
        bound_fa = parse_native_file(scorefiles[4], 'fa19')
        bound_fa23 = parse_native_file(scorefiles[5], 'fa23')
        
        pdbid = scorefiles[0].split('/')[-2]
        print(pdbid)
        # if(pdbid=='2ks1'):
        #     print('mp')
        #     print(bound_mp)
        #     print('fa19')
        #     print(bound_fa)
        #     print('fa23')
        #     print(bound_fa23)
        if pdb != None:
            if pdbid == pdb:
                make_plots(pdbid,  mp_scores, fa_scores, fa23_scores, bound_mp, bound_fa, bound_fa23, plot_dir)
        else:
            # print(fa_scores)
            make_plots(pdbid, mp_scores, fa_scores, fa23_scores, bound_mp, bound_fa, bound_fa23, plot_dir)
    

if __name__ == "__main__":
    main()