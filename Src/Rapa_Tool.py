from ast import Gt
from distutils.log import debug
from email.headerregistry import Group
#import logging
import Logger
import os
from datetime import datetime
import sys

import statsmodels.api as sm
from scipy.stats import shapiro
#log = logging.getLogger("my_logger")
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder
import seaborn as sns
from matplotlib_venn import venn2 , venn2_unweighted
from matplotlib.offsetbox import AnchoredText
import numpy as np
import pandas as pd
import LogisticRegression as LR
import matplotlib.gridspec as gridspec
ROOT = r"D:\PhD-Tomer"
LIB_PATH = "Tools_and_Libraries"
PROJECT_PATH = "Projects"
project_name = "Rapamycin"
DATA_SOURCE_FOLDER = "input"
OUTPUT_PARENT_FOLDER = "output"
runId = "Exp_1"

UNIVARIATE = True
#compared_groups=["Rapa_Live_Rapa_Dead","PBS_Live_PBS_Dead","Rapa_Dead_PBS_Dead","Rapa_Live_PBS_Live","Live_Dead"]
if UNIVARIATE:
    compared_groups = ["Live_Dead","PBS","Rap"]
else:
    compared_groups = ["Live_Dead"]
#GROUP = ["All"]
strains = ["X31","Vie1203","Cal09"]

loggerobj =Logger.get_Logger("Rapa_Logger")
import_Lib_path = os.path.join(ROOT,LIB_PATH)
sys.path.append(import_Lib_path)
from LielTools import FileTools 
from LielTools import DataTools
from LielTools import TimeTools
from LielTools import GLM_functions as glm
from SigniCorr import signi_corr
from LielTools import PlotTools as pT
#from LielTools import PlotTools
DEBUG = False
if DEBUG:
    alpha_range = [0.1,0.2,0.4, 0.7,1, 2, 4,5, 10, 20,30,50]
    l1_range =    [0.01, 0.05, 0.1,0.2, 0.4, 0.7, 0.99]
else:
    alpha_range = [0.01, 0.05, 0.1, 0.2, 0.4, 0.7, 1, 2, 5, 10, 50]
    l1_range = [0.01, 0.05, 0.1, 0.2, 0.4, 0.7, 0.99]


#keys for the dictionary of types of experiment to run
CO_I_FW = "co_i_fw"
CO_I_FD = "co_i_fd"
CO_I = "co_i"
CO_FW = "co_fw"
CO_FD = "co_fd"
CO = "co"
NCO_FW ="nco_fw"
NCO_FD = "nco_fd"
NCO = "nco"

#col_ignore = ["group_Rapa_Live","group_PBS_Live","group_Rapa_Dead"]

outcomevar = ["survival_group_Dead"]
#outcomeRap = ["group_Rapa_Dead"]
#outcomePbs = ["group_PBS_Dead"]
covar = ["Group_PBS"]
interactPvalueList = ["intract_Group_PBS_pvalue"]

""" expDict = {NCO:"exp_NoCovariate",NCO_FD:"exp_NoCovariate_FDR",NCO_FW:"exp_NoCovariate_FWER",
            CO:"exp_withCovariate",CO_FD:"exp_withCovariate_FDR",CO_FW:"exp_withCovariate_FWER",
            CO_I:"exp_withcovariate_int",CO_I_FD:"exp_withcovariate_int_FDR",CO_I_FW:"exp_withcovariate_int_FWER"}
 
expDict = {NCO:"exp_NoCovariate",NCO_FD:"exp_NoCovariate_FDR",NCO_FW:"exp_NoCovariate_FWER",
            CO_I:"exp_withcovariate_int",CO_I_FD:"exp_withcovariate_int_FDR",CO_I_FW:"exp_withcovariate_int_FWER"}
"""

expDict = {NCO:"exp_NoCovariate"}


standardizemethod = ["ShaprioTest"]

GLM_res_dict={}
#standardizemethod = ["No_logTranform"]
def get_str_of_date(time=False):
    """
    The function take a boolean paramter to whehter to generate a string current time or date
    The function returns a date string or time string of system time.

    Parameters
    -----------
    time:Boolean
        True or False whether to return system time or date in string format
    
    Return
    --------
    dtstr: string
        Return current date or time in format of string 
        date : DD_MM_YYYY
        time: HH_MM_SS
    
    """
    date = datetime.now()
    if(time):
        dtstr = date.strftime("%d")+"_"+date.strftime("%m")+"_"+date.strftime("%Y")
    else:
       dtstr = date.strftime("%H")+"_"+date.strftime("%M")+"_"+date.strftime("%S")
    
    return dtstr



def createExperimentDirectory(expPath,index):
    """
    """
    #name = childFolderLst[index]
    path = ""
    #get the latest folder
    lstdir = FileTools.get_subfolders(expPath,return_full_paths=True)
    # check if the dir contain any run folder 
    if len(lstdir) >0:
        latestFolderpath = FileTools.get_folder_time(lstdir)
        i = latestFolderpath.find("Exp_")
        lstexpid = latestFolderpath[i+4:]
        newID = int(lstexpid)+1
        newexpFolder = latestFolderpath[:i+4]+str(newID)
        path = FileTools.create_folder(newexpFolder)
    else:
        path=  FileTools.create_folder(os.path.join(expPath,runId))
    
    return path 

def get_item_Index_in_Lst(lst,search):
    idx = -1
    for i in range(len(lst)):
        if lst[i].find(search)>0:
            idx = i
            break
    return idx

def createOutputDirectories(output_path,protien_isotype,foldername,experimentFolderName):
    """
    """
    #before start of the run create folder with the run and # id R###
    output = os.path.join(output_path,protien_isotype)
    path = FileTools.create_folder(output)
    folderPath = FileTools.create_folder(os.path.join(path,foldername))

    #path = createExperimentDirectory(folderPath,-1)
    path = FileTools.create_folder(os.path.join(folderPath,experimentFolderName))
    return path 

def readInputFiles(file_path,file_name,read_index = None,keep_default_na=True):
    """
    """
    try:
        dataFrame = FileTools.read_excel(os.path.join(file_path,file_name),indexCol=read_index,keep_default_na=keep_default_na)
    except FileNotFoundError:
        raise
    return dataFrame

def writeOutputFile(file_path,filename,data,write_index = False):
    """
    """
    FileTools.write2Excel(os.path.join(file_path,filename),data,write_index)

def processDataFrame(df,col_ignore=[]):
    processdata = DataTools.get_df_without_cols(df,col_ignore)
    return processdata

def startuniAnalysis(protien_isotype,inputFilePath,experimentPath,exptype,groupType="",stdmethod="All"):
    """
    """
    #p = experimentPath[:experimentPath.find(expDict[exptype])]
    path = os.path.join(ROOT,PROJECT_PATH,project_name,OUTPUT_PARENT_FOLDER)
    processedFile = path+r"\{}".format(protien_isotype)
    #data = get_processed_data(protien_isotype,inputpath=inputFilePath,outputpath=p,stMethod=stdmethod,Gtype=groupType)
    data = get_processed_data(protien_isotype,inputpath=inputFilePath,stMethod=stdmethod,processedFileoutputpath = processedFile,Gtype=groupType,univariate = True)
    outcomevar[0]= "survival_group"
    #log the ratio of outcome variable 
    #print("The ratio of live:dead outcome is: \n",data[outcomevar[0]].value_counts())
    # below code only for 4 group comparsion
    '''
    if groupType == "Live_Dead":
        #loggerobj.info('dataframe - \n {}'.format(data[["survival_group","survival_group_bin"]].to_string()))
        outcomevar[0]= "survival_group_Dead"
    else:
        #loggerobj.info('dataframe - \n {}'.format(data[["group","group_bin"]].to_string()))
        outcomevar[0]= "Group_PBS"
    '''
    res = performGLM(data,groupType,exptype)
    outputFilename = "GLM_"+protien_isotype+"_"+stdmethod+"_"+groupType+".xlsx"
    writeOutputFile(experimentPath,outputFilename,res,write_index=True)
    
    return res

    
def multivaraite_analysis(protien_iso,outputpath,inputFilePath,stdmethod = "All",processedFilepath="",group_type="",standarise = False,
                          cv = False):
    """
    """
    data = get_processed_data(protien_iso,inputpath=inputFilePath,stMethod=stdmethod,processedFileoutputpath = processedFilepath,Gtype=group_type)
    if standarise:
        data = get_standardised_df(data)
    '''
    fig = signi_corr.plot_signi_corr(data,signif_by_method="FDR",figsize=(20,15),aster_size=5,
                    numbers_size=5,ticks_fontsize=10,title="correlation {}".format(protien_iso))
    fig_filename = "Correlation_{}.jpg".format(protien_iso)
    fig_path = os.path.join(outputpath,fig_filename)
    fig.savefig(fig_path, bbox_inches='tight')
    '''

    
    outcome = "survival_group"
    #remove the outcome column dataframe 
    x_cols_list = DataTools.list_removeItemIfExists(list(data.columns),outcome)
    
    common_peptide_list = get_common_peptides_between_LR_EFT(processedFilepath,protien_iso,stdmethod,group_type)
    #if cross_validation is false read from the already generated l1 and l2 value list 
    if not cv:
        l1_l2_df = readInputFiles(processedFilepath,"L1_L2_weights_{}_{}.xlsx".format(stdmethod,protien_iso),read_index=0)
    resdict= LR.glm_nested_loo(data,l1_l2_df,outcome, x_cols_list,outputpath,protien_iso,stdmethod,alpha_range,l1_range,5,bar_figsize=(50,25),heatmap_figsize=(30,20)
                    , common_peptide_highlight=common_peptide_list)
    GLM_res_dict.update({protien_iso:resdict})
    '''
    elasticnetTune = glm.tune_GLM_elastic_net(outputpath,model_df = data,y_col_name=outcome,x_cols_list = x_cols_list,
    model_name="GLM_{}_{}".format(protien_iso,stdmethod),alphas=alpha_range,l1s=l1_range,cv_folds=5)
    #loggerobj.debug("table for set of alpha:%s and l1:%s\n",alpha_range,l1_range)
    #loggerobj.debug("{}".format(elasticnetTune.to_string()))
    l1,l2 = get_l1_L2_tunned_parameter(protien_iso,stdmethod)
    glm.glm_LOO(outputpath, model_df = data,y_col_name=outcome,x_cols_list = x_cols_list, model_name="GLM_{}_{}".format(stdmethod,protien_iso),
            alpha=l1, L1_wt=l2, heatmap_annotate_text=True, logistic=True,bar_figsize=(20,18),heatmap_figsize=(20,18))
    '''

    
def performGLM(df,gtype,exptype):
    predicators = df.columns.tolist()
    std = True
    #sort the dataframe based on the index
    df = df.sort_index()
    if exptype == "nco":
        result = glm.GLMAnalysis(df,predicators,outcomevar,standardize = std)
        result = get_pvalueAdjustment(result,['pvalue'],"FDR")
        result = get_pvalueAdjustment(result,['pvalue'],"FWER")
    '''
    match exptype:
        case "nco":
            result = glm.GLMAnalysis(df,predicators,outcomevar,standardize = std)
        case "nco_fd":
            
            result = glm.GLMAnalysis(df,predicators,outcomevar,standardize = std)

            #result = get_pvalueAdjustment(result,"FDR")
            result = get_pvalueAdjustment(result,['pvalue'],"FDR")
            #pvalue = signi_corr.multipAdjustPvalsMat(result,['pvalue'],'FDR')
            
        case "nco_fw":
            result = glm.GLMAnalysis(df,predicators,outcomevar,standardize = std)
            result = get_pvalueAdjustment(result,['pvalue'],"FWER")
            #pvalue = signi_corr.multipAdjustPvalsMat(result,['pvalue'],'FWER')
        case "co":
            result = glm.GLMAnalysis(df,predicators,outcomevar,covar,interaction =False,standardize = std)
            #log.debug("----------Only covariate-----------------")
            #log.debug(checkSign_covaraiate_FromGLMModel(result,covar))
            
        case "co_fw":
            result = glm.GLMAnalysis(df,predicators,outcomevar,covar,interaction =False,standardize = std)
            pvalue = signi_corr.multipAdjustPvalsMat(result,['pvalue'],'FWER')
            
            #log.debug(checkSign_covaraiate_FromGLMModel(result,covar))
            
            result.update(pvalue)
        case "co_fd":
            result = glm.GLMAnalysis(df,predicators,outcomevar,covar,interaction =False,standardize = std)
            pvalue = signi_corr.multipAdjustPvalsMat(result,['pvalue'],'FDR')
            
            #log.debug(checkSign_covaraiate_FromGLMModel(result,covar))
           
            result.update(pvalue)
        case "co_i_fd":
            result = glm.GLMAnalysis(df,predicators,outcomevar,covar,interaction =True,standardize = std)
            result = get_pvalueAdjustment(result,interactPvalueList,"FDR")
            loggerobj.debug(checkSign_covaraiate_FromGLMModel(result,interactPvalueList,0.20))
        case "co_i_fw":
            result = glm.GLMAnalysis(df,predicators,outcomevar,covar,interaction =True,standardize = std)
            result = get_pvalueAdjustment(result,interactPvalueList,"FWER")
            loggerobj.debug(checkSign_covaraiate_FromGLMModel(result,interactPvalueList))
        case "co_i":
            result = glm.GLMAnalysis(df,predicators,outcomevar,covar,interaction =True,standardize = std)
            loggerobj.debug(checkSign_covaraiate_FromGLMModel(result,interactPvalueList))
    '''
    return result




def checkSign_covaraiate_FromGLMModel(df,CovariateList=[],significance_value=0.05):
    
    covariateSignPercentageDict = dict.fromkeys(CovariateList)
    for covariate in CovariateList:
        pvalDf = df[df[covariate]<significance_value]
        count = pvalDf.shape[0]
        #log_df_row(pvalDf)
        #log.debug("\n",pvalDf.iloc[:,:-1])
        loggerobj.info('dataframe - \n {}'.format(pvalDf.loc[:, pvalDf.columns != 'res'].to_string()))
        loggerobj.debug("*"*20)
        temp= pvalDf.loc[:, pvalDf.columns != 'res']
        for idx, idxname in enumerate(temp.index):
            for g in groups:
                if g == "PBS":
                    formula = f"PBS: {temp.ParamConst[idxname]} + {temp.ParamBeta[idxname]:.3f} *{idxname} + {temp.Group_PBS[idxname]:.3f} *PBS + {temp.intract_Group_PBS[idxname]:.3f} *{idxname} *PBS"
                else:
                    formula = f"RAP: {temp.ParamConst[idxname]} + {temp.ParamBeta[idxname]:.3f} *{idxname}"
                loggerobj.debug("Model:%s",formula)
        loggerobj.debug("*"*20)
        covariateSignPercentageDict[covariate] = (count/df.shape[0])*100
   
    return covariateSignPercentageDict

def draw_lollipop_chart(df,colList=[],fig_rows=4, fig_cols=5, figsize=(30, 20),
                        title='',
                        title_fontsize=20,
                        title_color='black',ticks_fontsize=14,
                        xlabel = "",ylabel="",yscale="linear", xlabelFont = 14,ylabelFont = 14,    
                        markerFmt="o",lineFormat="-",significance = "FWER",fig_save_path="",isAdjusted = False):

    
    fig, axes = plt.subplots(fig_rows, fig_cols, figsize=figsize,squeeze=False)

    for colname in colList:
        for row in range(fig_rows):
            for col in range(fig_cols):
                #sort the values based on the pvalues 
                y=df[colname].to_numpy()
                x = df.index.to_numpy()
                my_color = np.where(y>=0, 'orange', 'skyblue')
                markerline, stemline, baseline = axes[row,col].stem(x,y,markerfmt = markerFmt,linefmt = lineFormat)
                stemline.set(linewidth=1.5, color=my_color)
                #markerline.set(color=my_color)
                axes[row,col].set_xticklabels(x,rotation=(90),fontsize=ticks_fontsize)
                #axes[row,col].set_xticks(range(0,len(x)))
                axes[row,col].set_ylabel(ylabel)
                axes[row,col].set_xlabel(xlabel)
                #markerline
                for idx,idxname in enumerate(df.index):
                    if isAdjusted == False:
                        if df.pvalue[idxname] < 0.0005:
                            an = "***"
                        elif df.pvalue[idxname] <0.005:
                            an = "**"
                        elif df.pvalue[idxname] <0.05:
                            an = "*"
                        else:
                            an = ""
                    else:
                        if df.pvalue[idxname] < 0.0005 and df.pvalue_FDR[idxname] < 0.2:
                            an = "***"
                        elif df.pvalue[idxname] <0.005 and df.pvalue_FDR[idxname] < 0.2:
                            an = "**"
                        elif df.pvalue[idxname] <0.05 and df.pvalue_FDR[idxname] < 0.2:
                            an = "*"
                        else:
                            an = ""
                    if df.iloc[idx][colname] <0:
                        yidx = df.iloc[idx][colname]-.25
                    else:
                        yidx = df.iloc[idx][colname]+.75
                    text_color = "black"
                    axes[row,col].annotate(an,xy=(idx,yidx), weight='bold', ha='center',
                                     va='center', rotation=90, color=text_color)
                #axes[row,col].set_yticks(y)
                #axes[row,col].set_ylim(min(y)-1,max(y)+1)
    #set the size of the figure
    #split df based on negative and positive value 
    
    
    fig.suptitle(title, fontsize=title_fontsize)
    fig.tight_layout()
    plt.savefig(fig_save_path, bbox_inches='tight')
    plt.close("all")

    
def get_log_transformation(df,standardizeMethod = "NoTransform"):
    
    data = df.copy()
    col_name = df.columns
    if standardizeMethod == "TransformAll":
        for col in col_name:
            if data[col].dtype.name != "object":
                data.loc[data[col]>20,col] = np.log(data[col])
    elif standardizeMethod == "ShaprioTest":
        loggerobj.debug("Log tranformation based on shaprio-wilk test")
        print("shaprio-wilk test")
        for col in col_name:
            if data[col].dtype.name != "object":
                #H0 - that population  is normaly distributed
                #Ha - that population is not normally distributed
                if(shapiro(df[col]).pvalue<=0.05): #we reject H0
                    loggerobj.debug("peptide:{}".format(col))
                    data.loc[data[col]<20,col] = 20 
                    data[col] = np.log(data[col])
    return data
    
def get_shaprio_test(df):
    res ={}
    for col in df.columns:
        res.update({col:shapiro(df[col]).pvalue})
    return res

def get_standardised_df(df):
    standardizeFunc = lambda col: (col - np.nanmean(col)) / np.nanstd(col)
    stdDf = df.copy()
    
    for prec in stdDf.columns:
        if stdDf[prec].dtype.name != "object":
            if len(stdDf[prec].unique()) > 2:
                stdDf[[prec]] = stdDf[[prec]].apply(standardizeFunc)
    return stdDf

    
def check_for_normal_distribution(df):
    #figpath = os.path.join(experimentPath,protien_isotype+"BeforeTranform_qq"+".jpg")
    #tranpath = os.path.join(experimentPath,protien_isotype+"afterTranform_qq"+".jpg")
    plotDf = df.iloc[:,:30]
    plotDf = get_standardised_df(plotDf)
    nrows,ncols= plotDf.shape
    res = get_shaprio_test(plotDf)
    for key,value in res.items():    
        loggerobj.debug("peptide-%s-pvalue:%.4f",key,value)
    #pT.plot_QQ_plot(plotDf,figpath,fig_rows=ncols//6,fig_cols = 6,figsize=(30,20),title=protien_isotype+"Peptide distribution")
    logDf = get_log_transformation(plotDf)
    #logDf = get_standardised_df(logDf)
    #pT.plot_QQ_plot(logDf,tranpath,fig_rows=ncols//6,fig_cols = 6,figsize=(30,20),title=protien_isotype+"Peptide distribution")
    
    #pT.plot_columns_dist(plotDf,figpath,fig_rows=ncols//6,fig_cols = 6,figsize=(30,20),title=protien_isotype+"Peptide distribution")
    logDf = get_log_transformation(plotDf)
    #pT.plot_columns_dist(logDf,tranpath,fig_rows=ncols//6,fig_cols = 6,figsize=(30,20),title=protien_isotype+"Peptide distribution")


def get_pvalueAdjustment(df,column,Method="FWER"):
    data = df.copy()
    if Method == 'FDR':
        methodSM = 'fdr_bh'
    elif Method == 'FWER':
        methodSM = 'holm'

    for pval in column:
        m = pval+"_"+Method
        data[m] = np.nan
        for strain in strains:
            index = data.index.str.contains(strain)
            temp = data[index].copy()
            if temp.shape[0]!=0:
                pvalue = signi_corr.multipAdjustPvalsMat(temp[[pval]],Method)
                pvalue.columns = [m]
    
            #pvalsFlattened = temp[["intract_Group_RAP_pvalue"]].values
            #pval = pvalsFlattened.flatten()
            
            #temp.loc[:,Method] = sm.stats.multipletests(pval, method= methodSM)[1]
            #data.update(temp)
                #data.loc[:,m] = pvalue
                data.update(pvalue)
        data[m]=data[m].fillna(1)
    return data

def get_processed_data(protien_isotype,inputpath="",processedFileoutputpath="",logtransform=True,
                        stMethod="All",convert_categorical=True,Gtype="",univariate = False):

    
    #check if the file exist     
    optFilename = "GLM_"+protien_isotype+"_"+stMethod+"_"+Gtype+".xlsx"

    if os.path.exists(os.path.join(processedFileoutputpath,optFilename)):
        #read the out file
        df = readInputFiles(processedFileoutputpath,optFilename)
    else:
        df = get_df_with_filtered_peptide(protien_isotype,inputpath,Gtype)
        df = get_log_transformation(df,stMethod)

        if univariate:
            df = get_df_based_on_groups(df,Gtype,"Group")

    
        #below code is only needed for 4 group comparsions 
            '''
            if Gtype == "Live_Dead":
                df = DataTools.get_df_without_cols(df,["Group"])
                df =  DataTools.convert_CategoricalToBinary(df)
            else:
                df = DataTools.get_df_without_cols(df,["survival_group"])
                df = DataTools.convert_CategoricalToBinary(df)
            '''
            df = DataTools.get_df_without_cols(df,["group","Group"])
            #df = DataTools.convert_CategoricalToBinary(df)
            print(df['survival_group'].value_counts())
            df["survival_group"] = df["survival_group"].replace(["Live","Dead"],[1,0])
            #df['survival_group']=np.where(df['survival_group']=='Dead',0,1)
            print(df['survival_group'].value_counts())
        else:
            #df = DataTools.convert_CategoricalToBinary(df)
            df["survival_group"] = df["survival_group"].replace(["Live","Dead"],[1,0])
            df["Group"] = df["Group"].replace(["PBS","RAP"],[0,1])
        writeOutputFile(processedFileoutputpath,optFilename,df)
    
    return df


def get_df_with_filtered_peptide(protien_isotype,inputpath="",Gtype =""):
    """
    """
    filteredPeptideFilename = "filtered_peptide_lists_paired_"
    
    path = os.path.join(inputpath,"Array_"+protien_isotype)
    '''
    array_dataFilename = "array_data_RAPA3_{}.xlsx".format(protien_isotype)
    array_data = readInputFiles(path,array_dataFilename)
    '''
    array_data = get_rapa_array_df(inputpath,protien_isotype)
    #demographic columns to extract from the peptide array dataframe 
    demographycolnames = readInputFiles(inputpath,"Demographic.xlsx")
    col_names = demographycolnames.iloc[:,1].values.tolist()

    #get peptides from  the blind filtered list based on the groups

    if Gtype == "Live_Dead":
        filteredPeptideFilename = "filtered_Live_Dead_peptide_lists_{}.xlsx".format(protien_isotype)
    elif Gtype == "Rap":
        filteredPeptideFilename = "filtered_peptide_lists_paired_Rapa_Live_Rapa_Dead_{}.xlsx".format(protien_isotype)
    elif Gtype == "PBS":
        filteredPeptideFilename = "filtered_peptide_lists_paired_PBS_Live_PBS_Dead_{}.xlsx".format(protien_isotype)
    

    #below code is for the 4 group comparision
    ''''
    if Gtype == "Live_Dead":
        filteredPeptideFilename = "filtered_Live_Dead_peptide_lists_{}.xlsx".format(protien_isotype)
    else:
        filteredPeptideFilename = filteredPeptideFilename+Gtype+"_"+protien_isotype+".xlsx"
    '''
    filteredPeptidedf = readInputFiles(path,filteredPeptideFilename)
    filteredpeptideColName = filteredPeptidedf.iloc[:,1].values.tolist()    
    
    #in df only have the col present in the significant peptide file and demographic file 
    col_include = filteredpeptideColName+col_names
    array_data  = DataTools.get_df_with_cols(array_data,col_include)
    return array_data

def get_compared_groups_on_pvalue(expPath,protien_isotype,stMethod,adjustPvalue=False):
    col_name = ["Antigen","P-value","group_names","Adjusted_pvalue"]
    
    data = pd.DataFrame(columns= col_name)
    for g in compared_groups:

        tempdf = pd.DataFrame()
        filename = "GLM_"+protien_isotype+"_"+stMethod+"_"+g+".xlsx"
        df=  readInputFiles(expPath,filename,0)
        df = df.dropna()
        if adjustPvalue:
            df= df[df["pvalue_FDR"]<=0.2]
            tempdf["Adjusted_pvalue"] = df["pvalue_FDR"].to_numpy()
        else:
            df = df[df["pvalue"]<=0.05]
        tempdf["Antigen"] = df.index.to_list()
        tempdf["P-value"] = df["pvalue"].to_numpy()
        match g:
            case "Rapa_Live_Rapa_Dead":
                tempdf["group_names"] = ["Rapa_Live vs. Rapa_Dead"]*df.shape[0]
            case "PBS_Live_PBS_Dead":
                tempdf["group_names"] = ["PBS_Live vs. PBS_Dead"]*df.shape[0]
            case "Rapa_Dead_PBS_Dead":
                tempdf["group_names"] = ["Rapa_Dead vs. PBS_Dead"]*df.shape[0]
            case "Rapa_Live_PBS_Live":
                tempdf["group_names"] = ["Rapa_Live vs. PBS_Live"]*df.shape[0]
            case "Live_Dead":
                tempdf["group_names"] = ["Live vs. Dead"]*df.shape[0]
            case "PBS":
                tempdf["group_names"] = ["PBS_Live vs. PBS_Dead"]*df.shape[0]
            case "Rap":
                tempdf["group_names"] = ["Rapa_Live vs. Rapa_Dead"]*df.shape[0]

        data = pd.concat([data,tempdf])
    return data
            



def draw_venndiagram(df,dfcompare,label=("",""),fig_rows=4, fig_cols=5, figsize=(30, 20),fig_save_path = "",title = "",
                        title_fontsize =14,colormap = ("orange","blue"),item_fontsize = 8,isannotate = False):
    
    sns.set_palette(sns.color_palette("Set2"))

    #fig, axes = plt.subplots(fig_rows, fig_cols, figsize=figsize,squeeze=False)
    fig = plt.figure(figsize= figsize)
    subplot = 220
    for i in range(len(compared_groups)):
        subplot+=1
        group_name = get_comparisionname_from_group(compared_groups[i])
        #create a set for the particular group 
        LRGroup = set(df.loc[df["group_names"]==group_name,"Antigen"].to_numpy())
        EFTGroup = set(dfcompare.loc[dfcompare["group_names"]==group_name,"Antigen"].to_numpy())
        if i ==2:
            ax = plt.subplot2grid((2,1),(i//2,i%2),colspan=4)
        else:
            ax= fig.add_subplot(subplot)
        
        v = venn2_unweighted([LRGroup,EFTGroup],set_labels = label,ax = ax,set_colors=colormap,alpha=0.5)
        if isannotate:
            if(LRGroup-EFTGroup):
                v.get_label_by_id('10').set_text("\n".join(LRGroup-EFTGroup))
                v.get_label_by_id('10').set_fontsize(item_fontsize)
            else:
                v.get_label_by_id('10').set_text("")
            if(EFTGroup-LRGroup):
                v.get_label_by_id('01').set_text("\n".join(EFTGroup-LRGroup))
                v.get_label_by_id('01').set_fontsize(item_fontsize)
            if (LRGroup&EFTGroup):
                v.get_label_by_id('11').set_text("\n".join(LRGroup&EFTGroup))
                v.get_label_by_id('11').set_fontsize(item_fontsize)
        ax.set_title(label = group_name,loc="center",fontsize=20)
    
    '''
    gs = gridspec.GridSpec(2, 2)
    group_size = len(compared_groups)
    i = 0 
    lr_eft_intersection_dict = {}
    for row in range(fig_rows):
        for col in range(fig_cols):
             if (i<group_size):
                 group_name = get_comparisionname_from_group(compared_groups[i])
                #create a set for the particular group 
                 LRGroup = set(df.loc[df["group_names"]==group_name,"Antigen"].to_numpy())
                 EFTGroup = set(dfcompare.loc[dfcompare["group_names"]==group_name,"Antigen"].to_numpy())
                 lr_eft_intersection_dict.update({compared_groups[i]:LRGroup&EFTGroup})
                 ax = plt.subplot(gs[0, 0:2])
                 v = venn2_unweighted([LRGroup,EFTGroup],set_labels = label,ax = ax,set_colors=colormap,alpha=0.5)
                 if isannotate:
                    if(LRGroup-EFTGroup):
                        v.get_label_by_id('10').set_text("\n".join(LRGroup-EFTGroup))
                        v.get_label_by_id('10').set_fontsize(item_fontsize)
                    else:
                        v.get_label_by_id('10').set_text("")
                    if(EFTGroup-LRGroup):
                        v.get_label_by_id('01').set_text("\n".join(EFTGroup-LRGroup))
                        v.get_label_by_id('01').set_fontsize(item_fontsize)
                    if (LRGroup&EFTGroup):
                        v.get_label_by_id('11').set_text("\n".join(LRGroup&EFTGroup))
                        v.get_label_by_id('11').set_fontsize(item_fontsize)
                 axes[row][col].set_title(label = group_name,loc="center")
            
             else:
                fig.delaxes(axes[row,col])
             i+=1
    '''
    fig.suptitle(title, fontsize=title_fontsize)
    fig.tight_layout()
    fig.savefig(fig_save_path, bbox_inches='tight')
    plt.close("all")
    #return lr_eft_intersection_dict




def get_comparisionname_from_group(groupname):
    name =""
    match groupname:
            case "Rapa_Live_Rapa_Dead":
                name = "Rapa_Live vs. Rapa_Dead"
            case "PBS_Live_PBS_Dead":
                name = "PBS_Live vs. PBS_Dead"
            case "Rapa_Dead_PBS_Dead":
                name = "Rapa_Dead vs. PBS_Dead"
            case "Rapa_Live_PBS_Live":
                name = "Rapa_Live vs. PBS_Live"
            case "Live_Dead":
                name = "Live vs. Dead"
            case "PBS":
                name = "PBS_Live vs. PBS_Dead"
            case "Rap":
                name = "Rapa_Live vs. Rapa_Dead"
    
    return name


def get_fisher_test_sig_peptide(expPath,expType,protien_isotype,stMethod):

        eftpath = os.path.join(ROOT,PROJECT_PATH,project_name,DATA_SOURCE_FOLDER)
        if (expType =="nco"):
            df = get_compared_groups_on_pvalue(expPath,protien_isotype,stMethod)
        else:
            df = get_compared_groups_on_pvalue(expPath,protien_isotype,stMethod,True)
        path = os.path.join(eftpath,"Array_"+protien_isotype)
        
        #p = os.path.join(expPath,"All_sig_peptide.xlsx")
        writeOutputFile(expPath,"All_sig_peptide_{}.xlsx".format(stMethod),df)

def plot_venn_diagram(expPath,expType,protien_isotype,stMethod):
    #temp need to change this is future 
    
        eftpath = os.path.join(ROOT,PROJECT_PATH,project_name,DATA_SOURCE_FOLDER)
        filename = os.path.join(expPath,"All_sig_peptide_{}.xlsx".format(stMethod))
        path = os.path.join(eftpath,"Array_"+protien_isotype)
        ''''
        if (expType =="nco"):
            df = get_compared_groups_on_pvalue(expPath,protien_isotype,stMethod)
        else:
            df = get_compared_groups_on_pvalue(expPath,protien_isotype,stMethod,True)
       
        
        #create an excel file 
        #p = os.path.join(expPath,"All_sig_peptide.xlsx")
        #writeOutputFile(expPath,"All_sig_peptide_{}.xlsx".format(stMethod),df)
        '''
        df = readInputFiles(expPath,"All_sig_peptide_{}.xlsx".format(stMethod))
        sigLive_dead = "all_sig_peptides_Live_Dead_peptide_{}.xlsx".format(protien_isotype)
        sig_all_4_comp = "all_sig_4_comp_peptides_{}.xlsx".format(protien_isotype)
        sigLive_deaddf = readInputFiles(path,sigLive_dead,None)
        sig_all_4_compdf = readInputFiles(path,sig_all_4_comp,None)
        df2 = pd.concat([sig_all_4_compdf,sigLive_deaddf])
        title = "Significant peptide {}".format(protien_isotype)
        #create a folder figures
        fig_path = FileTools.create_folder(os.path.join(expPath,"Figures"))
        fig_path = os.path.join(fig_path,"Sig_peptide_venn_{}_{}.jpg".format(stMethod,protien_isotype))
        return draw_venndiagram(df,df2,label=("Logistic Regression","Fisher's exact Test"),fig_rows=3,fig_cols=2,title=title,fig_save_path=fig_path,item_fontsize=6,isannotate=True)

def plot_lollipop(expPath,expType,protien_isotype,stMethod):
    
    for g in compared_groups:
        filename = "GLM_"+protien_isotype+"_"+stMethod+"_"+g+".xlsx"
        comparisonGroup = get_comparisionname_from_group(g)
        df=  readInputFiles(expPath,filename,0)
        fig_filename = "sig_peptides_lollipop_{}_{}_{}.jpg".format(stMethod,protien_isotype,comparisonGroup)
        fig_path = FileTools.create_folder(os.path.join(expPath,"Figures"))
        fig_path = os.path.join(fig_path,fig_filename)
        plottitle = "significant peptide {} {}".format(protien_isotype,comparisonGroup)
        df = df.sort_values(by="ParamBeta").dropna() #drop rows which have beta values as nan
            #sort df based on pvalues 
        if (expType =="nco"):
            draw_lollipop_chart(df,["ParamBeta"],fig_rows=1,fig_cols=1,title=plottitle,yscale="symlog",ticks_fontsize=10,fig_save_path=fig_path)

        else:
            draw_lollipop_chart(df,["ParamBeta"],fig_rows=1,fig_cols=1,title=plottitle,yscale="symlog",ticks_fontsize=10,fig_save_path=fig_path,isAdjusted=True)

        #draw_lollipop_chart(df,["ParamBeta"],fig_rows=1,fig_cols=1,title=plottitle,yscale="symlog",ticks_fontsize=10,fig_save_path=fig_path)


def get_last_createdFolder(path):
     #folderLst = FileTools.get_subfolders(path)
     return  FileTools.get_latest_folder_from_path(path)

def get_labels_forcategorical(df,col_name):

    d = df.copy()
    le = LabelEncoder()

    label = le.fit_transform(d[col_name])

    #loggerobj.debug("converted col_name:%s and {}}",col_name,label)
    d = d.drop(col_name,axis=1)
    col_name = col_name+"_bin"
    d[col_name] = label

    return d

def get_df_based_on_groups(df,group,col_name):
    
    data = df.copy()
    match group:
        case "Rapa_Live_Rapa_Dead":
            data = data[(data[col_name]=="Rapa_Dead") | (data[col_name]=="Rapa_Live")]
        case "PBS_Live_PBS_Dead":
            data = data[(data[col_name]=="PBS_Dead") | (data[col_name]=="PBS_Live")]
        case "Rapa_Dead_PBS_Dead":
            data = data[(df[col_name]=="Rapa_Dead")| (data[col_name]=="PBS_Dead")]
        case "Rapa_Live_PBS_Live":
            data = data[(df[col_name]=="Rapa_Live") | (data[col_name]=="PBS_Live")]
        case "Rap":
            data = data[(df[col_name]=="RAP") | (data[col_name]=="RAP")]
        case "PBS":
            data = data[(df[col_name]=="PBS") | (data[col_name]=="PBS")]
    
    return data

    

def get_rapa_array_df(array_data_path,protien_isotype):

    path = os.path.join(array_data_path,"Array_"+protien_isotype)
    array_dataFilename = "array_data_RAPA3_{}.xlsx".format(protien_isotype)
    array_data = readInputFiles(path,array_dataFilename)
    return array_data

def rapa_draw_df_box_plot(x_col,y_col_names,modeldf,df = None,output_file_path=None, fig_rows=4, fig_cols=5, figsize=(30, 20),
                          title='', title_fontsize=18, title_y=1.03, font_scale=1, sns_style='ticks',hue_col=None,plot_with_boder = [],
                          **boxplot_kwargs):
    """
    """
    plt.close('all')
    sns.set(font_scale=font_scale)
    sns.set_style(sns_style)

    num_figs = len(y_col_names)
    if fig_cols * fig_rows < num_figs:
        print('plot_boxplot_subplots: number of columns', num_figs, 'is smaller than fig_cols*fig_rows')

    i = 0
    fig, axes = plt.subplots(fig_rows, fig_cols, figsize=figsize)
    for row in range(fig_rows):
        for col in range(fig_cols):
            if (i < num_figs):
                pT.plot_boxplot(df[x_col], df[y_col_names[i]],
                             seriesHue=None if hue_col is None else df[hue_col],
                             ax=axes[row, col],
                             xy_title_fontsize=10,xRotation=0,figsize=(15,12))
                #get the pvalues (FDR and FWER and not corrected)
                pvalue = modeldf.at[y_col_names[i],"pvalue"]
                fdr = modeldf.at[y_col_names[i],"pvalue_FDR"]
                fwer = modeldf.at[y_col_names[i],"pvalue_FWER"]

                #annotate this values of the fig
                annText = "P: {:.4f}\nP-FDR: {:.4f}\nP-FWER: {:.4f}".format(pvalue,fdr,fwer)

                at = AnchoredText(annText, prop=dict(size=7), frameon=False, loc='upper right')
                #at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
                axes[row,col].add_artist(at)
                axes[row,col].tick_params(axis='x', labelsize=8)
                axes[row,col].tick_params(axis='y', labelsize=8)
                if y_col_names[i] in plot_with_boder:
                    axes[row,col].spines['bottom'].set_color('red')
                    axes[row,col].spines['top'].set_color('red')
                    axes[row,col].spines['left'].set_color('red')
                    axes[row,col].spines['right'].set_color('red')
                    '''
                    ax = axes[row,col].axis()
                    rec = plt.Rectangle((ax[0] - 0.7, ax[2] - 0.5), (ax[1] - ax[0]) + 1, (ax[3] - ax[2]) + 0.4, fill=False, lw=2, edgecolor="red")
                    rec = axes[row,col].add_patch(rec)
                    rec.set_clip_on(False)
                    '''
                '''
                x = df[y_col_names[i]].max()
                text_color = "black"
                axes[row,col].annotate(annText,xy=(x-200,100), weight='bold', ha='center',
                                     va='center', rotation=90, color=text_color)
                '''
                i = i + 1
            else:
                fig.delaxes(axes[row,col])
    fig.suptitle(title, fontsize=title_fontsize, y=title_y)
    fig.tight_layout()

    if output_file_path is not None:
        fig.savefig(output_file_path, bbox_inches='tight', dpi=500)
    plt.close("all")
    return axes



def get_l1_L2_tunned_parameter(protien,logTranform):
    l1 = 0.0
    l2 = 0.0
    if protien == "HA_IgM":
        if logTranform == "ShaprioTest":
            l1 = 10.0
            l2 = 0.005
        elif logTranform == "TransformAll":
            l1 = 5.0
            l2 = 0.005
        else:
            l1 = 0.01
            l2 = 0.05
    elif protien == "HA_IgG":
        if logTranform == "ShaprioTest":
            l1 = 0.01
            l2 = 0.05
            
        elif logTranform == "TransformAll":
            l1 = 0.01
            l2 = 0.99
        else:
            l1 = 0.001
            l2 = 0.4
    
    elif protien == "NA_IgM":
        if logTranform == "ShaprioTest":
            l1 = 0.1
            l2 = 0.05
        elif logTranform == "TransformAll":
            l1 = 0.005
            l2 = 0.7
        else:
            l1 = 0.001
            l2 = 0.7
    else:
         if logTranform == "ShaprioTest":
            l1 = 0.001
            l2 = 0.01
         elif logTranform == "TransformAll":
            l1 = 0.005
            l2 = 0.99
         else:
            l1 = 0.2
            l2 = 0.01

    return l1,l2

def get_common_peptides_between_LR_EFT(df1,df2,group):
    #temp need to change this is future 
    
        '''
        eftpath = os.path.join(ROOT,PROJECT_PATH,project_name,DATA_SOURCE_FOLDER)
        filename = os.path.join(expPath,"All_sig_peptide_{}.xlsx".format(stMethod))
        path = os.path.join(eftpath,"Array_"+protien_isotype)
        
        if (expType =="nco"):
            df = get_compared_groups_on_pvalue(expPath,protien_isotype,stMethod)
        else:
            df = get_compared_groups_on_pvalue(expPath,protien_isotype,stMethod,True)
       
        
        #create an excel file 
        #p = os.path.join(expPath,"All_sig_peptide.xlsx")
        #writeOutputFile(expPath,"All_sig_peptide_{}.xlsx".format(stMethod),df)
        
        df = readInputFiles(expPath,"All_sig_peptide_{}.xlsx".format(stMethod))
        sigLive_dead = "all_sig_peptides_Live_Dead_peptide_{}.xlsx".format(protien_isotype)
        #sig_all_4_comp = "all_sig_4_comp_peptides_{}.xlsx".format(protien_isotype)
        sigLive_deaddf = readInputFiles(path,sigLive_dead,0)
        #sig_all_4_compdf = readInputFiles(path,sig_all_4_comp,None) 
        '''
        group_name = get_comparisionname_from_group(group)
        #create a set for the particular group 
        LRGroup = set(df1.loc[df1["group_names"]==group_name,"Antigen"].to_numpy())
        EFTGroup = set(df2.loc[df2["group_names"]==group_name,"Antigen"].to_numpy())
        return (LRGroup,EFTGroup)
