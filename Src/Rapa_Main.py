
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
from LielTools import DataTools
import seaborn as sns 
import matplotlib.pyplot as plt
import pandas as pd

 
JUST_RUN_ONE_PROTIEN = False
draw_boxplot = False
draw_distplot = False
#plotpath = r"D:\PhD-Tomer\Projects\Rapamycin\output\{}\16_06_2022".format(protien)

if JUST_RUN_ONE_PROTIEN: 
    #protien_iso= ["NA_IgM"]
    protien_isotype= ["HA_NA_IgG_IgM"]

else:
    protien_iso= ["HA_IgG","HA_IgM","NA_IgG","NA_IgM"]
    protien_isotype= ["HA_NA_IgG","HA_NA_IgM","HA_NA_IgG_IgM"]



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

checkwithsingelmodel = True


def combine_IgG_IgM(cv):
    for pro in protien_isotype:
        prot = pro.split("_")
        datestr = rt.get_str_of_date(True)
        timeStr = rt.get_str_of_date()
       
        if pro != "HA_NA_IgG_IgM":
        #create folder with data inside the protien iso_top
            iso = prot[2]
            glmoutputpath = rt.FileTools.create_folder(os.path.join(os.path.join(output_path,"{}_{}_{}".format(prot[0],prot[1],iso)),datestr,timeStr))
            if not checkwithsingelmodel:
                ha_path = os.path.join(input_path,"Array_{}_{}".format(prot[0],iso))
                na_path = os.path.join(input_path,"Array_{}_{}".format(prot[1],iso))
                df_ha = readInputFiles(ha_path,"array_data_RAPA3_{}_{}.xlsx".format(prot[0],iso))
                df_na = readInputFiles(na_path,"array_data_RAPA3_{}_{}.xlsx".format(prot[1],iso))
            #df_IgG_IgM = df_igG.set_index("ptid").join(df_IgM.set_index("ptid"))
                df_HA_NA = df_ha.join(df_na,how="left",rsuffix='_right')
                rt.writeOutputFile(os.path.join(output_path,"{}_{}_{}".format(prot[0],prot[1],iso)),"array_data_RAPA3_{}_{}_IgG.xlsx".format(prot[0],prot[1]),df_HA_NA)
                filtered_peptide_ha_df = readInputFiles(ha_path,"filtered_Live_Dead_peptide_lists_{}_{}.xlsx".format(prot[0],iso))
                filtered_peptide_na_df = readInputFiles(na_path,"filtered_Live_Dead_peptide_lists_{}_{}.xlsx".format(prot[1],iso))
                demographycolnames = readInputFiles(input_path,"Demographic.xlsx")
                col_names = demographycolnames.iloc[:,1].values.tolist()
                filteredpeptide_ha = filtered_peptide_ha_df.iloc[:,1].values.tolist() 
                filteredpeptide_na = filtered_peptide_na_df.iloc[:,1].values.tolist()
                cols_include = filteredpeptide_ha+filteredpeptide_na+col_names
                df  =  DataTools.get_df_with_cols(df_HA_NA,cols_include)
                df = rt.get_log_transformation(df,rt.standardizemethod[0])
                df["survival_group"] = df["survival_group"].replace(["Live","Dead"],[1,0])
                df["Group"] = df["Group"].replace(["PBS","RAP"],[0,1])
                outpath = os.path.join(output_path,"{}_{}_{}".format(prot[0],prot[1],iso))
                rt.writeOutputFile(outpath,"log_tranformed_data_RAPA3_{}_{}_{}.xlsx".format(prot[0],prot[1],iso),df)
                rt.multivaraite_analysis("{}_{}_{}".format(prot[0],prot[1],iso),glmoutputpath,stdmethod=rt.standardizemethod[0],group_type=rt.compared_groups[0],standarise=True
                                    ,cv=True,data=df)
            else:
                ha_na_path = os.path.join(output_path,"{}_{}_{}".format(prot[0],prot[1],iso))
                df_HA_NA  = readInputFiles(ha_na_path,"log_tranformed_data_RAPA3_{}_{}_{}.xlsx".format(prot[0],prot[1],iso))
                df_l1_l2_weights = readInputFiles(ha_na_path,"L1_L2_weights_ShaprioTest_{}_{}_{}.xlsx".format(prot[0],prot[1],iso))
                
                '''
                df_ha_mean_weights = readInputFiles(ha_na_path,"params_stats____GLM_ShaprioTest_{}_{}.xlsx".format(prot[0],iso),read_index=0)
                df_na_mean_weights = readInputFiles(ha_na_path,"params_stats____GLM_ShaprioTest_{}_{}.xlsx".format(prot[1],iso),read_index=0)
                glm_res_less = df_ha_mean_weights.loc[(df_ha_mean_weights["Mean"]<-0.15)]
                glm_res_more = df_ha_mean_weights.loc[(df_ha_mean_weights["Mean"]>0.15)]
                df_ha_mean_weights = pd.concat([glm_res_less,glm_res_more])
                glm_res_less = df_na_mean_weights.loc[(df_na_mean_weights["Mean"]<-0.15)]
                glm_res_more = df_na_mean_weights.loc[(df_na_mean_weights["Mean"]>0.15)]
                df_na_mean_weights = pd.concat([glm_res_less,glm_res_more])
                df_features = df_ha_mean_weights.index.to_list()+df_na_mean_weights.index.to_list()+["survival_group","Group","weight"]
                df_HA_NA  =  DataTools.get_df_with_cols(df_HA_NA,df_features)
                '''
                rt.multivaraite_analysis("{}_{}_{}".format(prot[0],prot[1],iso),glmoutputpath,stdmethod=rt.standardizemethod[0],group_type=rt.compared_groups[0],standarise=True
                                    ,cv=cv,data=df_HA_NA,l1_l2_weight_df = df_l1_l2_weights )
        else:
            if not checkwithsingelmodel:
                glmoutputpath = rt.FileTools.create_folder(os.path.join(os.path.join(output_path,"{}_{}_{}_{}".format(prot[0],prot[1],prot[2],prot[3])),datestr,timeStr))
                ha_na_path_igg = os.path.join(output_path,"{}_{}_{}".format(prot[0],prot[1],prot[2]))
                ha_na_path_igm = os.path.join(output_path,"{}_{}_{}".format(prot[0],prot[1],prot[3]))
                df_ha_na_igg= readInputFiles(ha_na_path_igg,"log_tranformed_data_RAPA3_{}_{}_{}.xlsx".format(prot[0],prot[1],prot[2]))
                df_ha_na_igm= readInputFiles(ha_na_path_igm,"log_tranformed_data_RAPA3_{}_{}_{}.xlsx".format(prot[0],prot[1],prot[3]))
                df_HA_NA_igg_igm = df_ha_na_igg.join(df_ha_na_igm,how="left",rsuffix='_IgM')
                df_HA_NA_igg_igm = df_HA_NA_igg_igm.drop(["survival_group_IgM","Group_IgM","weight_IgM"],axis=1)
                rt.writeOutputFile(os.path.join(output_path,"{}_{}_{}_{}".format(prot[0],prot[1],prot[2],prot[3])),"log_tranformed_data_RAPA3_{}_{}_IgG_IgM.xlsx".format(prot[0],prot[1])
                                ,df_HA_NA_igg_igm)
                rt.multivaraite_analysis("{}_{}_{}_{}".format(prot[0],prot[1],prot[2],prot[3]),glmoutputpath,stdmethod=rt.standardizemethod[0],group_type=rt.compared_groups[0],standarise=True
                                    ,cv=True,data=df_HA_NA_igg_igm)
            else:
                ha_na_igg_igm_path = os.path.join(os.path.join(output_path,"{}_{}_{}_{}".format(prot[0],prot[1],prot[2],prot[3])))
                df_HA_NA_igg_igm = readInputFiles(ha_na_igg_igm_path,"log_tranformed_data_RAPA3_{}_{}_IgG_IgM.xlsx".format(prot[0],prot[1]))
                df_l1_l2_weights = readInputFiles(ha_na_igg_igm_path,"L1_L2_weights_ShaprioTest_{}_{}_IgG_IgM.xlsx".format(prot[0],prot[1]))
                '''
                df_ha_mean_weights_igg = readInputFiles(ha_na_igg_igm_path,"params_stats____GLM_ShaprioTest_{}_IgG.xlsx".format(prot[0]),read_index=0)
                df_na_mean_weights_igg = readInputFiles(ha_na_igg_igm_path,"params_stats____GLM_ShaprioTest_{}_IgG.xlsx".format(prot[1]),read_index=0)
                df_ha_mean_weights_igm = readInputFiles(ha_na_igg_igm_path,"params_stats____GLM_ShaprioTest_{}_IgM.xlsx".format(prot[0]),read_index=0)
                df_na_mean_weights_igm = readInputFiles(ha_na_igg_igm_path,"params_stats____GLM_ShaprioTest_{}_IgM.xlsx".format(prot[1]),read_index=0)
                glm_res_less = df_ha_mean_weights_igg.loc[(df_ha_mean_weights_igg["Mean"]<-0.15)]
                glm_res_more = df_ha_mean_weights_igg.loc[(df_ha_mean_weights_igg["Mean"]>0.15)]
                df_ha_mean_weights_igg = pd.concat([glm_res_less,glm_res_more])
                glm_res_less = df_na_mean_weights_igg.loc[(df_na_mean_weights_igg["Mean"]<-0.15)]
                glm_res_more = df_na_mean_weights_igg.loc[(df_na_mean_weights_igg["Mean"]>0.15)]
                df_na_mean_weights_igg = pd.concat([glm_res_less,glm_res_more])
                glm_res_less = df_ha_mean_weights_igm.loc[(df_ha_mean_weights_igm["Mean"]<-0.15)]
                glm_res_more = df_ha_mean_weights_igm.loc[(df_ha_mean_weights_igm["Mean"]>0.15)]
                df_ha_mean_weights_igm = pd.concat([glm_res_less,glm_res_more])
                glm_res_less = df_na_mean_weights_igm.loc[(df_na_mean_weights_igm["Mean"]<-0.15)]
                glm_res_more = df_na_mean_weights_igm.loc[(df_na_mean_weights_igm["Mean"]>0.15)]
                df_na_mean_weights_igm = pd.concat([glm_res_less,glm_res_more])
                df_features = df_ha_mean_weights_igg.index.to_list()+df_na_mean_weights_igg.index.to_list()+df_ha_mean_weights_igm.index.to_list()+df_na_mean_weights_igm.index.to_list()+["survival_group","Group","weight"]
                df_HA_NA_igg_igm  =  DataTools.get_df_with_cols(df_HA_NA_igg_igm,df_features)
                '''
                rt.multivaraite_analysis("{}_{}_{}_{}".format(prot[0],prot[1],prot[2],prot[3]),glmoutputpath,stdmethod=rt.standardizemethod[0],group_type=rt.compared_groups[0],standarise=True
                                    ,cv=cv,data=df_HA_NA_igg_igm,l1_l2_weight_df = df_l1_l2_weights)

                

    res_df = pd.DataFrame(columns=["Auc","FPR","TPR"],index = protien_isotype)
    
    res_df.index.set_names('Protien-Iso')
    for protien_name in rt.GLM_res_dict.keys():
        res_df.at[protien_name,"Auc"] = rt.GLM_res_dict[protien_name]['roc_auc_info']['auc']
        res_df.at[protien_name,"FPR"] = rt.GLM_res_dict[protien_name]['roc_auc_info']['rates']['false positive rate']
        res_df.at[protien_name,"TPR"] = rt.GLM_res_dict[protien_name]['roc_auc_info']['rates']['true positive rate']
        
    rt.writeOutputFile(output_path,"RoC_AUC_protiens_iso.xlsx",res_df,write_index = True)
#rapamain(isUnivariate=True)
combine_IgG_IgM(False)

def perform_GLM(protien_isotype,inputfilename,no_features=[],model_name="",outputfilename="",l1=0.01,l2=0.05):
    rootpath = os.path.join(rt.ROOT,rt.PROJECT_PATH,rt.project_name,rt.OUTPUT_PARENT_FOLDER)
    path = os.path.join(rootpath,protien_isotype)
    df_ha_na = rt.readInputFiles(path,inputfilename)
    df_ha_na = rt.get_standardised_df(df_ha_na)

    outcome = "survival_group"
    x_cols_list = rt.DataTools.list_removeItemIfExists(list(df_ha_na.columns),outcome)

    glmoutputpath = os.path.join(path,"NoteBookOutput")
    

    no_of_feature_with_auc_dict = dict.fromkeys(no_of_features)

    final_res_auc = pd.DataFrame(columns=["Auc","FPR","TPR"],index = no_of_features)
    
    final_res_auc.index.set_names('Sigfeatures')
   
    for num_features in no_features:
    
        resdict = rt.glm.glm_LOO(glmoutputpath, df_ha_na, y_col_name=outcome, x_cols_list = x_cols_list, model_name = model_name+"_{}".format(num_features), 
                                heatmap_figsize=(30,20), bar_figsize=(50,25),
                                alpha=l1, L1_wt=l2, heatmap_annotate_text=True, logistic=True,featureSelection=True,
                                num_feature = num_features)

        no_of_feature_with_auc_dict.update({num_features:resdict})

#write the all the auc details to a excel sheet

    for feature_num in no_of_feature_with_auc_dict.keys():
        final_res_auc.at[feature_num,"Auc"] = no_of_feature_with_auc_dict[feature_num]['roc_auc_info']['auc']
        final_res_auc.at[feature_num,"FPR"] = no_of_feature_with_auc_dict[feature_num]['roc_auc_info']['rates']['false positive rate']
        final_res_auc.at[feature_num,"TPR"] = no_of_feature_with_auc_dict[feature_num]['roc_auc_info']['rates']['true positive rate']

    rt.writeOutputFile(path,outputfilename,final_res_auc,write_index = True)

no_of_features = [10,20,30,40,50,60,70,80,90,100,110,120]
#perform_GLM("HA_IgG","GLM_HA_IgG_ShaprioTest_Live_Dead.xlsx",no_of_features,"GLM_shaprio_HA_IgG_0.01_0.2_FS","RoC_AUC_for_different_features_by_importance_loo.xlsx",l2=0.2)
perform_GLM("HA_IgM","GLM_HA_IgM_ShaprioTest_Live_Dead.xlsx",no_of_features,"GLM_shaprio_HA_IgM_0.01_0.2_FS","RoC_AUC_for_different_features_by_importance_loo.xlsx",l2=0.2)
