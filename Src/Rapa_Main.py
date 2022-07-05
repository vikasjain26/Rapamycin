
import Logger
import Rapa_Tool as rt
import os 
import sys
from sklearn.feature_selection import GenericUnivariateSelect, chi2,mutual_info_regression 
from sklearn.linear_model import LogisticRegression
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
import numpy as np
import logging

import_Lib_path = os.path.join(rt.ROOT,rt.LIB_PATH)
sys.path.append(import_Lib_path)
from LielTools import PlotTools as pT

isplot = True
JUST_RUN_ONE_PROTIEN = True
draw_boxplot = True
draw_distplot = True
#plotpath = r"D:\PhD-Tomer\Projects\Rapamycin\output\{}\16_06_2022".format(protien)
if JUST_RUN_ONE_PROTIEN: 
    protien_iso= ["NA_IgG"]

else:
    protien_iso= ["HA_IgM","HA_IgG","NA_IgM","NA_IgG"]


isUnivariate = False
loggername = "Rapa_Logger"
logObj = Logger.get_Logger(loggername)
format = "%(asctime)s- %(levelname)s - %(message)s"
logObj.setLevel(logging.DEBUG)
output_path = os.path.join(rt.ROOT,rt.PROJECT_PATH,rt.project_name,rt.OUTPUT_PARENT_FOLDER)
input_path = os.path.join(rt.ROOT,rt.PROJECT_PATH,rt.project_name,rt.DATA_SOURCE_FOLDER)
def rapamain():

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
            univariateAnalysis(protien,datedir,timeStr)
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



def univariateAnalysis(protien,DateDirName,timeStr):

    for expType in rt.expDict:
        for s in rt.standardizemethod:
            for g in rt.GROUP:
                logObj.debug("start analysis for st method %s and exp %s and for group %s ",s,expType,g)
                cexpdir = rt.FileTools.create_folder(os.path.join(DateDirName,timeStr,rt.expDict[expType]))
                modelres = rt.startuniAnalysis(protien,input_path,cexpdir,expType,stdmethod = s,groupType=g)
                if isplot:
                    drawing_plots(modelres,cexpdir,expType,pro_iso=protien,standardmethod =s,group=g)



def start_multivariate_analysis(dateDirectory,TimeStr,protien_iso):
    """
    """
    for s in rt.standardizemethod:
        for g in rt.GROUP:
            logObj.debug("start analysis for st method %s ",s)
            cexpdir = rt.FileTools.create_folder(os.path.join(dateDirectory,TimeStr))
            processedFile = output_path+r"\{}".format(protien_iso)
            rt.multivaraite_analysis(protien_iso,inputFilePath=input_path,outputpath=cexpdir,stdmethod = s,processedFilepath=processedFile,group_type=g)


def drawing_plots(modeldf,fig_save_folder,experiment_type,standardmethod,pro_iso,group):
    """
    """
    #figsave_folder = rt.get_last_createdFolder(fig_save_folder)
    #figsave_folder = figsave_folder+r"\{}".format(rt.expDict[experiment_type])
    #rt.plot_venn_diagram(cexpdir,expType,protien,s)
            #rt.plot_lollipop(cexpdir,expType,protien,s)
    
    array_df = rt.get_rapa_array_df(input_path,pro_iso)
    y_cols = modeldf[modeldf["pvalue"]<=0.05].index.to_list()
    logObj.debug("significant peptide for {} for {} with tranformation method {} are {}".format(pro_iso,group,get_log_method_str(standardmethod),y_cols))
    plot_title = "Logistic Regression for Significant peptide with {} for {}".format(get_log_method_str(standardmethod),group)
    fig_path = rt.FileTools.create_folder(os.path.join(fig_save_folder,"Figures"))
    if draw_boxplot:
        fig_filename = "sig_peptide_box_plot_logregression{}_{}_{}.jpg".format(standardmethod,pro_iso,group)
        fig_box_path = os.path.join(fig_path,fig_filename)
        rt.rapa_draw_df_box_plot("survival_group",y_cols,modeldf,df = array_df,output_file_path=fig_box_path, fig_rows=5, fig_cols=7, figsize=(25, 20),
                          title=plot_title, title_fontsize=18, title_y=1.03, font_scale=1, sns_style='ticks',hue_col=None)
    if draw_distplot:
        fig_filename = "sig_peptide_dist_logregression{}_{}_{}.jpg".format(standardmethod,pro_iso,group)
        fig_dist_path = os.path.join(fig_path,fig_filename)
        fig = rt.pT.plot_columns_dist(array_df[y_cols], 
                                  fig_rows=5, fig_cols=7, figsize=(30, 20), 
                                  bins=20, font_scale=1.4,
                                  kde_color='black', hist_color='purple', hist_alpha=0.3,
                                  title=plot_title, title_fontsize=18, title_y=1.03,
                                  rug=True, rug_color='black', rug_alpha=0.3, rug_linewidth=1, rug_height=0.03, 
                                  output_file_path=fig_dist_path)



def get_log_method_str(standardmethod):
    if standardmethod == "ShaprioTest":
        return "log tranform based on  shaprio-wilk test "
    elif standardmethod == "TransformAll":
        return "log tranform all"
    elif standardmethod == "No_logTranform":
        return "No log transformation"
    else:
        return ""
rapamain()
'''
def checkinput_arraydata_with_significantPepList():

    for proiso in protien_iso:
        
        path = input_path+r"\Array_{}".format(proiso)
        filename = "array_data_RAPA3_{}.xlsx".format(proiso)
        for cgroups in rt.compared_groups:
            if cgroups == "Live_Dead":
                f = "all_sig_peptides_Live_Dead_peptide_{}.xlsx".format(proiso)
            else:
                f = "all_sig_4_comp_peptides_{}.xlsx".format(proiso)
            df = rt.readInputFiles(path,filename,None)
            comparedf = rt.readInputFiles(path,f,None)

            match cgroups:
                case "Rapa_Live_Rapa_Dead":
                    df = df[(df["group"]=="Rapa_Dead") | (df["group"]=="Rapa_Live")]
                    col_name = comparedf.loc[(comparedf["group_names"]=="Rapa_Live vs. Rapa_Dead"),"Antigen"].to_numpy()
                case "PBS_Live_PBS_Dead":
                    df = df[(df["group"]=="PBS_Dead") | (df["group"]=="PBS_Live")]
                    col_name = comparedf.loc[(comparedf["group_names"]=="PBS_Live vs. PBS_Dead"),"Antigen"].to_numpy()
                case "Rapa_Dead_PBS_Dead":
                    df = df[(df["group"]=="Rapa_Dead")| (df["group"]=="PBS_Dead")]
                    col_name = comparedf.loc[(comparedf["group_names"]=="Rapa_Dead vs. PBS_Dead"),"Antigen"].to_numpy()
                case "Rapa_Live_PBS_Live":
                    df = df[(df["group"]=="Rapa_Live") | (df["group"]=="PBS_Live")]
                    col_name = comparedf.loc[(comparedf["group_names"]=="Rapa_Live vs. PBS_Live"),"Antigen"].to_numpy()
                case "Live_Dead":
                    col_name = comparedf.loc[(comparedf["group_names"]=="Live vs. Dead"),"Antigen"].to_numpy()
            print(col_name)
            for i in range(len(col_name)):
                x= df["group"]
                print(col_name[i])
                y = df[col_name[i]]
                fig_file= output_path+r"\{}\\".format(proiso)+r"\box_plot_{}".format(cgroups)
                pT.plot_boxplot(x,y,saveFolder=fig_file)
            col_name = np.append(col_name,"group")
            df = df[col_name]
            fig_file= output_path+r"\{}\\".format(proiso)+"{}_{}.jpg".format(proiso,cgroups)
            pT.plot_columns_dist_hue(df,"group",output_file_path=fig_file)
  

#checkinput_arraydata_with_significantPepList()
'''
