from __future__ import annotations
import Rapa_Tool as rt
import os 
import sys
from Rapa_Tool import readInputFiles
from sklearn.feature_selection import GenericUnivariateSelect, chi2,mutual_info_regression 
from sklearn.linear_model import LogisticRegression
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
import matplotlib.patches as mpatches
import numpy as np
import logging
import pickle

protien_iso= ["HA_IgG","HA_IgM","NA_IgG","NA_IgM"]
import_Lib_path = os.path.join(rt.ROOT,rt.LIB_PATH)
sys.path.append(import_Lib_path)
from LielTools import PlotTools
from LielTools import StatsTools
from LielTools import DataTools
import seaborn as sns 
import matplotlib.pyplot as plt
import pandas as pd 
output_path = os.path.join(rt.ROOT,rt.PROJECT_PATH,rt.project_name,rt.OUTPUT_PARENT_FOLDER)
input_path = os.path.join(rt.ROOT,rt.PROJECT_PATH,rt.project_name,rt.DATA_SOURCE_FOLDER)
barplot_color = ['#66cccc',
                 # '#2D7DD2',
                 '#009999',
                 '#cc0000',
                 # '#F45D01',
                 '#800000']

def auc_bar_plot(df,ycol,fig_size = (8, 6), color_pallette = None, plot_save_path=None, 
                 value_label_fontsize = 19, value_labels_rotation = 0, value_label_color = "#4d4d4d"
                 , xtick_rotation = 0, xtick_fontsize = 22, ytick_fontsize = 20
                 , xlabel_title = "", ylabel_title = "", xylabel_fontsize = 25
                 ,title = "",title_fonsize = 20, title_font = "bold",ylim=(0.5, 1)
                 , y_error=None , show_bar_labels = False,
                 color_specific_labels=None,
                 add_legend = False,Legend_title = "",Legend_label =[]):
    plt.figure(figsize=fig_size)
    plt.rcParams['font.family'] = 'Arial'
# sns.set_palette(['#472255', '#2f3c5e', '#175467', '#006f71'])
# sns.set_palette(['#6C5B7B', '#C06C84', '#F67280', '#F8B195'])
    if color_pallette is None:
        sns.set_palette(sns.color_palette("Set2"))
    else:
        sns.set_palette(color_pallette)

    if y_error is not None:
        ax = sns.barplot(y=df[ycol], x=df.index,yerr=df[y_error].values)
    else:
        ax = sns.barplot(y=df[ycol], x=df.index)

    if show_bar_labels:
        PlotTools.bar_plot_add_value_labels(ax, spacing=2, float_num_digits=2,
                              fontsize=value_label_fontsize, value_labels_rotation=value_labels_rotation,
                              color= value_label_color)

    if color_specific_labels is not None:
        for set_operation in color_specific_labels.keys():
                for color,pep_list in color_specific_labels[set_operation].items():
                     for xticklabel in ax.get_xticklabels():
                        xticklabel_text = str(xticklabel.get_text())
                        for label_from_list in pep_list:
                            if xticklabel_text == str(label_from_list):
                                 xticklabel.set_color(color)
                                 xticklabel.set_weight("bold")
    
    if add_legend:
        LR = mpatches.Patch(color='lime', label=Legend_label[0])
        LR_and_FET = mpatches.Patch(color='darkgreen', label=Legend_label[1])
        plt.legend(handles=[LR,LR_and_FET], title = Legend_title, prop={'size':20,"weight":"bold"}, loc='upper left')
    plt.xticks(rotation=xtick_rotation,fontsize = 15)
    ax.set_xticks(np.arange(len(df.index)))
    plt.yticks(fontsize = 15,weight="bold")
    plt.ylim(ylim)
    plt.rc('xtick', labelsize=xtick_fontsize)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=ytick_fontsize)  # fontsize of the tick labels
    plt.ylabel(ylabel_title, fontdict={'size': xylabel_fontsize})
    plt.xlabel(xlabel_title, fontdict={'size': xylabel_fontsize})
# plt.grid(axis='y', color='#808080', linestyle=':', linewidth=2.5)
    plt.grid(axis='y', color='#cccccc', linestyle=':', linewidth=1.5)
    ax = plt.gca()
    ax.tick_params(axis=u'both', which=u'both', length=0)
    plt.title(title, fontweight=title_font, fontsize=title_fonsize)
   
    plt.tight_layout()
    if plot_save_path is None:
        plt.show()
    else:
        plt.savefig(plot_save_path, dpi=300)
    plt.close("all")


def Roc_Auc_plot(df,fig_size = (8,6),xlabel="",ylabel="",xylabel_fontsize=22,
                 title = "",title_fonsize = 20, title_font = "bold",
                 plot_save_path = None):
    fig = plt.figure(figsize=fig_size)
    plt.rcParams['font.family'] = 'Arial'
    for index in df.index:
        x = np.fromstring(df.loc[index]['FPR'][1:-1], dtype=float, sep=' ')
        y = np.fromstring( df.loc[index]['TPR'][1:-1], dtype=float, sep=' ')
        plt.plot(x,y, label="{}, AUC={:.3f}".format(index, df.loc[index]['Auc']))
    plt.plot([0, 1], [0, 1], linestyle='--',color="orange")
    plt.xticks(np.arange(0.0,1,0.1),fontsize = 15)
    plt.yticks(np.arange(0.0, 1.1, step=0.1),fontsize = 15,weight="bold")
    plt.xlabel(xlabel,fontdict={"size":xylabel_fontsize})
    plt.ylabel(ylabel,fontdict={"size":xylabel_fontsize})
    plt.title(title, fontweight=title_font, fontsize=title_fonsize)
    plt.legend(prop={'size':18,"weight":"bold"}, loc='lower right')
    plt.tight_layout()
    if plot_save_path is None:
        plt.show()
    else:
        plt.savefig(plot_save_path, dpi=300)
    plt.close("all")



def get_data(filename,file_path,pro_iso,index=None):
     dir = os.path.join(file_path,pro_iso)
     df = readInputFiles(dir,filename,index)
     return df
#results_df = rt.readInputFiles(os.path.join(output_path,"HA_NA_IgG_IgM"),"RoC_AUC_for_different_features_by_importanc.xlsx",0)
results_df = rt.readInputFiles(output_path,"RoC_AUC_protiens_summaryforall.xlsx",0)
def plot_summary_stats_auc():

    auc_bar_plot(results_df, ycol = "Auc",plot_save_path= os.path.join(output_path,"Final Figures","Bar plot ROC AUC summary 1.jpg"), xtick_rotation = 45
             , xlabel_title = "Model", ylabel_title = "AUC",
             title="AUC Summary", show_bar_labels = False,fig_size=(8,6))
    Roc_Auc_plot(results_df,plot_save_path= os.path.join(output_path,"Final Figures","ROC AUC Summary 6.jpg"),
             xlabel = "False Positive Rate", ylabel="True Positive Rate",title= "Roc AUC summary", xylabel_fontsize = 25,fig_size=(8,6))
    




def plot_peptide_params():
    for pro_iso in protien_iso:
        eft_df=readInputFiles(os.path.join(input_path,"Array_"+pro_iso),"all_sig_peptides.xlsx",0)
        lr_df = get_data("All_sig_peptide_ShaprioTest.xlsx",output_path,pro_iso)
        lr_df,eft_df = rt.get_common_peptides_between_LR_EFT(lr_df,eft_df,"Live_Dead")
        intersection_pep = list((lr_df&eft_df))
        diff_pep = list((lr_df-eft_df))
    
        color_specfic_label_dict = {}
        color_specfic_label_dict.update({"LR_int_EFT":{"darkgreen":intersection_pep},"LR_sub_EFT":{"lime":diff_pep}})
        
        glm_result_df = rt.readInputFiles(os.path.join(output_path,pro_iso),"params_stats____GLM_ShaprioTest_{}.xlsx".format(pro_iso),0)
        glm_res_less = glm_result_df.loc[(glm_result_df["Mean"]<-0.15)]
        glm_res_more = glm_result_df.loc[(glm_result_df["Mean"]>0.15)]
        glm_result_df = pd.concat([glm_res_less,glm_res_more])
        print(glm_result_df.shape[0])
        pro_iso = pro_iso.replace("_","-")
        pro = pro_iso.split("-")
        title_text = 'y=Survival, Protien={}, Iso-type = {}, Auc={}, N =38.'.format(pro[0],pro[1],np.round(results_df.at[pro_iso,"Auc"], 3))
        auc_bar_plot(glm_result_df, ycol = "Mean",plot_save_path= os.path.join(output_path,"Final Figures","stats bar plot {} 1.jpg".format(pro_iso)), xtick_rotation = 90
             , xlabel_title = "Features (Peptides)", ylabel_title = "Mean Feature (peptide) Weight in CV Models", ylim = (-1.2,1.2), fig_size=(30,20),
             title=title_text,y_error="Std",color_specific_labels=color_specfic_label_dict
             ,add_legend = True,Legend_label=["Logistic Regression"," LR intersection Fisher's"])



def plot_LR_vs_FET(compareGroup):

    venn_dict={}
    for pro_iso in protien_iso:
        eft_df=readInputFiles(os.path.join(input_path,"Array_"+pro_iso),"all_sig_peptides.xlsx",0)
        lr_df = get_data("All_sig_peptide_ShaprioTest.xlsx",output_path,pro_iso)
        lr_df,eft_df = rt.get_common_peptides_between_LR_EFT(lr_df,eft_df,compareGroup)
        venn_dict.update({pro_iso:[lr_df,eft_df]})
    #plot venn diagram 
    group_name = rt.get_comparisionname_from_group(compareGroup)
    title = "Significant peptides {}".format(group_name)
    rt.draw_venndiagram(venn_dict,label=("Logistic Regression","Fisher's exact Test"),fig_rows=2,fig_cols=2,
                         title=title,fig_save_path=os.path.join(output_path,"Final Figures","Sig_peptide_venn_{} 1.jpg".format(group_name)),item_fontsize=10,
                         isannotate=True,
                         title_fontsize=25,colormap=("pink","lightblue"),figsize=(30, 20))

def spearman():
    
        for pro_iso in protien_iso:
            df = get_data("GLM_{}_ShaprioTest_Live_Dead.xlsx".format(pro_iso),output_path,pro_iso=pro_iso)
            original_df = get_data("array_data_RAPA3_{}.xlsx".format(pro_iso),input_path,pro_iso="Array_"+pro_iso)
            original_df = DataTools.get_df_with_cols(original_df,df.columns)
            #original_df = original_df.drop("survival_group",axis=1)
            PlotTools.plot_boxplot_subplots(original_df,"survival_group",original_df.columns[30:60],output_file_path = output_path+"\{}_survival_box_2.jpg".format(pro_iso,pro_iso),
                              fig_rows=6, fig_cols=5,xy_title_fontsize=10,xRotation=0)
        

        #corrmat = StatsTools.getCorrelationMat(original_df,"spearman")
      


def get_short_name_for_strain(protien,strain):
    if protien[0:2]=="HA":
        if strain == "X31":
            return "H3"
        elif strain == "Cal09":
            return "H1"
        else:
            return "H5"
    else:
        if strain == "X31":
          return "N2"
        elif strain == "Cal09":
            return "Cal_N1"
        else:
            return "Vie_N1"


def get_protien_length(protien,strain):
    a_file = open("C:\\Users\Vikas jain\\Downloads\\align.pkl", "rb") 
    output = pickle.load(a_file)
    pro_len = 0
    if protien[0:2]=="HA":
        if strain == "X31":
            pro_len = max(output["H3"].keys())
        elif strain == "Cal09":
            pro_len = max(output["H1"].keys())
        else:
            pro_len = max(output["H5"].keys())
    else:
        if strain == "X31":
            pro_len = max(output["N2"].keys())
        elif strain == "Cal09":
            pro_len = max(output["Cal_N1"].keys())
        else:
            pro_len = max(output["Vie_N1"].keys())
    a_file.close()
    return pro_len



def get_scores_protien_structure():
    """     for pro_iso in protien_iso:
        ax = PlotTools.plot_heatmap(original_df.corr(),figsize=(40,30),xRotation = 90,yRotation=0,annotate_text = True,
                                    font_scale = 1.2)
        plt.savefig(output_path+"\{}_heatmap.jpg".format(pro_iso),dpi=500)
        
        plt.close("all")
        #path = os.path.join(input_path,"Array_"+pro_iso)
        #df2 = readInputFiles(path,"all_sig_peptides.xlsx",0)
        title = "Significant peptide {}".format(pro_iso)
        '''
        #bar plot of the params
        #get comman peptides in EFT and LR
        eft_df=readInputFiles(os.path.join(input_path,"Array_"+pro_iso),"all_sig_peptides.xlsx",0)
        lr_df = get_data("All_sig_peptide_ShaprioTest.xlsx",output_path,pro_iso)
        #lr_df,eft_df = rt.get_common_peptides_between_LR_EFT(lr_df,eft_df,"Live_Dead")
        #intersection_pep = list((lr_df&eft_df))
        #diff_pep = list((lr_df-eft_df))
        ''''
        color_specfic_label_dict = {}
        color_specfic_label_dict.update({"LR_int_EFT":{"darkgreen":intersection_pep},"LR_sub_EFT":{"lime":diff_pep}})
        '''
        glm_result_df = rt.readInputFiles(os.path.join(output_path,pro_iso),"params_stats____GLM_ShaprioTest_{}.xlsx".format(pro_iso),0)
        glm_res_less = glm_result_df.loc[(glm_result_df["Mean"]<-0.15)]
        glm_res_more = glm_result_df.loc[(glm_result_df["Mean"]>0.15)]
        glm_result_df = pd.concat([glm_res_less,glm_res_more])
        print(glm_result_df.shape[0])
        '''
        title_text = 'y=Survival, Auc={}.'.format(np.round(results_df.at[pro_iso,"Auc"], 3))
        auc_bar_plot(glm_result_df, ycol = "Mean",plot_save_path= os.path.join(output_path,"Final Figures","stats bar plot {}.jpg".format(pro_iso)), xtick_rotation = 90
             , xlabel_title = "Variable", ylabel_title = "Mean Variable Weight in CV Models", ylim = (-1.2,1.2), fig_size=(30,20),
             title=title_text,y_error="Std",color_specific_labels=color_specfic_label_dict)
        
        title = "Significant peptide {}".format(pro_iso)
        rt.draw_venndiagram(lr_df,eft_df,label=("Logistic Regression","Fisher's exact Test"),fig_rows=3,fig_cols=2,
                         title=title,fig_save_path=os.path.join(output_path,"Final Figures","Sig_peptide_venn_{}.jpg".format(pro_iso)),item_fontsize=10,isannotate=True,
                         title_fontsize=25,colormap=("pink","lightblue"))
        '''
        """
    for pro_iso in protien_iso:
        glm_result_df = rt.readInputFiles(os.path.join(output_path,pro_iso),"params_stats____GLM_ShaprioTest_{}.xlsx".format(pro_iso),0)
        glm_res_less = glm_result_df.loc[(glm_result_df["Mean"]<-0.15)]
        glm_res_more = glm_result_df.loc[(glm_result_df["Mean"]>0.15)]
        glm_result_df = pd.concat([glm_res_less,glm_res_more])
        overlap_df = pd.DataFrame(columns=["scores"],index=rt.strains)
        if pro_iso[0:2]=="HA":
            resdict = dict.fromkeys(["H3","H5","H1"])
        else:
            resdict = dict.fromkeys(["N2","Vie_N1","Cal_N1"])
        for strain in rt.strains:
            df =  glm_result_df[glm_result_df.index.str.contains(strain)].copy()
            df.set_index(df.index.str.removeprefix(pro_iso[0:2]+"_"+strain+"_"),inplace=True)
            df.index = df.index.astype(int)
            prolen = get_protien_length(pro_iso,strain)
            score_dict=dict.fromkeys(range(1,prolen+1))
            for k, v in score_dict.items():
                if v is None:
                    score_dict[k] = 0.0
            for pos in df.index:
                lastpos = pos+20
                overlap_list = sorted(df.index[(df.index >pos) &(df.index<=lastpos)].tolist())
                templist = []
                for overlap in overlap_list:
                    templist.append(df.at[overlap,"Mean"])
                if len(overlap_list)==0:
                    score = df.at[pos,"Mean"]
                    for index in range(pos,lastpos):
                        score_dict.update({index:score})
                else:
                #all of the overlap is positive 
                    templist.append(df.at[pos,"Mean"])
                    if all(val > 0 for val in templist):
                        score = max(templist)
                        for index in range(pos,lastpos):
                            score_dict.update({index:score})
                    elif all(val < 0 for val in templist):
                        score = max(templist)
                        for index in range(pos,lastpos):
                            score_dict.update({index:score})
                    else:
                        for index in range(pos,overlap_list[0]):
                            score_dict.update({index:df.at[pos,"Mean"]})
        
        resdict.update({get_short_name_for_strain(pro_iso,strain):score_dict})
     

    a_file = open(os.path.join(output_path,pro_iso,"overlap_{}.pkl".format(pro_iso)), "wb")
    pickle.dump(resdict, a_file)
    a_file.close()
    #rt.writeOutputFile(os.path.join(output_path,pro_iso),"overlap_{}_1.xlsx".format(pro_iso),overlap_df,True) 
        


#for g in rt.compared_groups:
    #plot_LR_vs_FET(g)
#plot_peptide_params()
plot_summary_stats_auc()
#spearman()

