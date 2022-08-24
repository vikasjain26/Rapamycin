
import Logger
import Rapa_Tool as rt
import os 
import sys
from Rapa_Tool import readInputFiles
from sklearn.feature_selection import GenericUnivariateSelect, chi2,mutual_info_regression 
from sklearn.linear_model import LogisticRegression
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
import numpy as np
import logging

import_Lib_path = os.path.join(rt.ROOT,rt.LIB_PATH)
sys.path.append(import_Lib_path)
from LielTools import PlotTools
import seaborn as sns 
import matplotlib.pyplot as plt
import pandas as pd


JUST_RUN_ONE_PROTIEN = False
draw_boxplot = False
draw_distplot = False
#plotpath = r"D:\PhD-Tomer\Projects\Rapamycin\output\{}\16_06_2022".format(protien)

if JUST_RUN_ONE_PROTIEN: 
    protien_iso= ["NA_IgM"]

else:
    protien_iso= ["HA_IgG","HA_IgM","NA_IgG","NA_IgM"]



loggername = "Rapa_Logger"
logObj = Logger.get_Logger(loggername)
format = "%(asctime)s- %(levelname)s - %(message)s"
logObj.setLevel(logging.DEBUG)
output_path = os.path.join(rt.ROOT,rt.PROJECT_PATH,rt.project_name,rt.OUTPUT_PARENT_FOLDER)
input_path = os.path.join(rt.ROOT,rt.PROJECT_PATH,rt.project_name,rt.DATA_SOURCE_FOLDER)
def rapamain(isUnivariate = False):

    logginFormat = "%(asctime)s- %(levelname)s - %(message)s"
    #logger = logging.getLogger("my_logger")
   
    groupSplit = False


    for protien in protien_iso:
        datestr = rt.get_str_of_date(True)
        timeStr = rt.get_str_of_date()
        #create folder with data inside the protien iso_top
        datedir = rt.createOutputDirectories(output_path,protien,datestr,"")
        #create a folder withgroups inside the date folder
        filename = protien+"_"+datestr+"_"+timeStr+".txt"
        loggerfilename = os.path.join(datedir,filename)
        #Logger.setup_logger(logObj,loggername,loggerfilename,logging.DEBUG)
        File_handler = logging.FileHandler(loggerfilename)
        formatter = logging.Formatter(format)
        File_handler.setFormatter(formatter)
        logObj.addHandler(File_handler)
        logObj.debug("*"*20)
        
        # do expirment for groups 
        #figfolder = rt.FileTools.create_folder(os.path.join(datedir,timeStr))
        #rt.startAnalysis(protien,input_path,figfolder,"",arrayFileName,"",selectedPep,demoFile)
        if isUnivariate:
            logObj.debug("Begin univariate analysis for protien iso_type %s",protien)
            univariateAnalysis(protien,datedir,timeStr,True)
        else:
            logObj.debug("Begin multivariate analysis for protien iso_type %s",protien)
            start_multivariate_analysis(datedir,timeStr,protien)
        """
       for expType in rt.expDict:
            for s in rt.standardizemethod:
                for g in rt.compared_groups:
                    logObj.debug("start analysis for st method %s and exp %s ",s,expType)
                    cexpdir = rt.FileTools.create_folder(os.path.join(datedir,timeStr,rt.expDict[expType]))
                    logObj.debug("start analysis for experiment %s",expType)
                    rt.startAnalysis(protien,input_path,cexpdir,expType,stdmethod = s,groupType=g)
                
            if isplot:
                
                path = rt.get_last_createdFolder(datedir)
                cexpdir = path+r"\{}".format(rt.expDict[expType])
                rt.plot_venn_diagram(cexpdir,expType,protien,s)
                #rt.plot_lollipop(cexpdir,expType,protien,s) """
            
    logObj.debug("Analysis End for iso_type %s",protien)
    logObj.debug("*"*20)
    File_handler.close()
    logObj.removeHandler(File_handler)
    res_df = pd.DataFrame(columns=["Auc","FPR","TPR"],index = protien_iso)
    res_df.index.set_names('Protien')
    for protien_name in rt.GLM_res_dict.keys():
        res_df.at[protien_name,"Auc"] = rt.GLM_res_dict[protien_name]['roc_auc_info']['auc']
        res_df.at[protien_name,"FPR"] = rt.GLM_res_dict[protien_name]['roc_auc_info']['rates']['false positive rate']
        res_df.at[protien_name,"TPR"] = rt.GLM_res_dict[protien_name]['roc_auc_info']['rates']['true positive rate']
        

    rt.writeOutputFile(output_path,"RoC_AUC_protiens.xlsx",res_df,write_index = True)



def univariateAnalysis(protien,DateDirName,timeStr,drawplots=False):

    for expType in rt.expDict:
        for s in rt.standardizemethod:
            modelResDict = {}
            for g in rt.compared_groups:
                logObj.debug("start analysis for st method %s and exp %s and for group %s ",s,expType,g)
                cexpdir = rt.FileTools.create_folder(os.path.join(DateDirName,timeStr,rt.expDict[expType]))
                modelres = rt.startuniAnalysis(protien,input_path,cexpdir,expType,stdmethod = s,groupType=g)
                modelResDict.update({g:modelres})
            rt.get_fisher_test_sig_peptide(expPath=cexpdir, expType=expType, protien_isotype= protien,stMethod=s)
            if drawplots:
                intersection_dict = drawing_plots_venn_lollipop(modelres,os.path.join(output_path,protien,DateDirName,timeStr),expType,
                                                                pro_iso=protien,standardmethod =s,group=g)
               # for g in rt.compared_groups:
                    #draw_plots_bar_dist(modelResDict[g],os.path.join(output_path,protien,DateDirName,timeStr),expType,standardmethod =s,
                                       # pro_iso=protien,group=g,common_pep_list =intersection_dict[g])

                

            



def start_multivariate_analysis(dateDirectory,TimeStr,protien_iso):
    """
    """
    for s in rt.standardizemethod:
        for g in rt.compared_groups:
            logObj.debug("start analysis for st method %s ",s)
            cexpdir = rt.FileTools.create_folder(os.path.join(dateDirectory,TimeStr))
            processedFile = output_path+r"\{}".format(protien_iso)
            rt.multivaraite_analysis(protien_iso,inputFilePath=input_path,outputpath=cexpdir,stdmethod = s,processedFilepath=processedFile,group_type=g,standarise=True)


def drawing_plots_venn_lollipop(modeldf,fig_save_folder,experiment_type,standardmethod,pro_iso,group):
    """
    """
    figsave_folder = rt.get_last_createdFolder(fig_save_folder)
    #figsave_folder = figsave_folder+r"\{}".format(rt.expDict[experiment_type])
    fig_path = rt.FileTools.create_folder(os.path.join(figsave_folder,"Figures"))
    lr_eft_intersection_dict = rt.plot_venn_diagram(figsave_folder,experiment_type,pro_iso,standardmethod)
    rt.plot_lollipop(figsave_folder,experiment_type,pro_iso,standardmethod)
    
    '''
    array_df = rt.get_rapa_array_df(input_path,pro_iso)
    y_cols = modeldf[modeldf["pvalue"]<=0.05].index.to_list()
    logObj.debug("significant peptide for {} for {} with tranformation method {} are {}".format(pro_iso,group,get_log_method_str(standardmethod),y_cols))
    plot_title = "Logistic Regression for Significant peptide with {} for {}".format(get_log_method_str(standardmethod),group)
   
    for g in rt.compared_groups:
        if draw_boxplot:
            fig_filename = "sig_peptide_box_plot_logregression{}_{}_{}.jpg".format(standardmethod,pro_iso,g)
            fig_box_path = os.path.join(fig_path,fig_filename)
            rt.rapa_draw_df_box_plot("survival_group",y_cols,modeldf,df = array_df,output_file_path=fig_box_path, fig_rows=5, fig_cols=7, figsize=(25, 20),
                          title=plot_title, title_fontsize=18, title_y=1.03, font_scale=1, sns_style='ticks',hue_col=None,plot_with_boder=lr_eft_intersection_dict[g])
        if draw_distplot:
            fig_filename = "sig_peptide_dist_logregression{}_{}_{}.jpg".format(standardmethod,pro_iso,g)
            fig_dist_path = os.path.join(fig_path,fig_filename)
            fig = rt.pT.plot_columns_dist(array_df[y_cols], 
                                  fig_rows=5, fig_cols=7, figsize=(30, 20), 
                                  bins=20, font_scale=1.4,
                                  kde_color='black', hist_color='purple', hist_alpha=0.3,
                                  title=plot_title, title_fontsize=18, title_y=1.03,
                                  rug=True, rug_color='black', rug_alpha=0.3, rug_linewidth=1, rug_height=0.03, 
                                  output_file_path=fig_dist_path)
    '''
    return lr_eft_intersection_dict


def draw_plots_bar_dist(modeldf,fig_save_folder,experiment_type,standardmethod,pro_iso,group,common_pep_list =[]):

    figsave_folder = rt.get_last_createdFolder(fig_save_folder)
    #figsave_folder = figsave_folder+r"\{}".format(rt.expDict[experiment_type])
    fig_path = rt.FileTools.create_folder(os.path.join(figsave_folder,"Figures"))
    array_df = rt.get_rapa_array_df(input_path,pro_iso)
    y_cols = modeldf[modeldf["pvalue"]<=0.05].index.to_list()
    logObj.debug("significant peptide for {} for {} with tranformation method {} are {}".format(pro_iso,group,get_log_method_str(standardmethod),y_cols))
    plot_title = "Logistic Regression for Significant peptide with {} for {}".format(get_log_method_str(standardmethod),group)
   
    if draw_boxplot:
        fig_filename = "sig_peptide_box_plot_logregression{}_{}_{}.jpg".format(standardmethod,pro_iso,group)
        fig_box_path = os.path.join(fig_path,fig_filename)
        rt.rapa_draw_df_box_plot("survival_group",y_cols,modeldf,df = array_df,output_file_path=fig_box_path, fig_rows=5, fig_cols=7, figsize=(25, 20),
                          title=plot_title, title_fontsize=18, title_y=1.03, font_scale=1, sns_style='ticks',hue_col=None,plot_with_boder=common_pep_list)
    if draw_distplot:
        fig_filename = "sig_peptide_dist_logregression{}_{}_{}.jpg".format(standardmethod,pro_iso,group)
        fig_dist_path = os.path.join(fig_path,fig_filename)
        fig = rt.pT.plot_columns_dist(array_df[y_cols], 
                                  fig_rows=5, fig_cols=7, figsize=(30, 20), 
                                  bins=20, font_scale=1.4,
                                  kde_color='black', hist_color='purple', hist_alpha=0.3,
                                  title=plot_title, title_fontsize=18, title_y=1.03,
                                  rug=True, rug_color='black', rug_alpha=0.3, rug_linewidth=1, rug_height=0.03, 
                                  output_file_path=fig_dist_path,plot_with_boder=common_pep_list)
    

def get_log_method_str(standardmethod):
    if standardmethod == "ShaprioTest":
        return "log tranform based on  shaprio-wilk test "
    elif standardmethod == "TransformAll":
        return "log tranform all"
    elif standardmethod == "No_logTranform":
        return "No log transformation"
    else:
        return ""



rapamain(isUnivariate=True)

