import pandas as pd
import numpy as np
import os 
import sys
from sklearn.linear_model import LogisticRegression
import matplotlib.pyplot as plt
ROOT = r"D:\PhD-Tomer"
LIB_PATH = "Tools_and_Libraries"
import_Lib_path = os.path.join(ROOT,LIB_PATH)
sys.path.append(import_Lib_path)

from LielTools import GLM_functions as glm
from LielTools import FileTools 
def log_regression(df,y_col_name, x_cols_list):
    """
    """


def glm_nested_loo(df,l1_l2weights_df,y_col_name, x_cols_list,outpath,protien_isotype,lgtransmethod,alpharange,l1range,cvfolds,
                    heatmap_figsize=(12, 8), bar_figsize=(10,8),heatmap_annotate_text=True,cross_valid= False,common_peptide_highlight = []):
    """
    """
    alphaL1dict = {}
    GLM_res_list, GLM_test_res_list = [], []
    l1_wt,l2_wt=[],[]
    for i_index_test in range(len(df.index)):
        print("Test point index :",i_index_test)
        ind_train = list(range(len(df.index)))
        ind_train.remove(i_index_test)
        ind_train = np.array(ind_train)
        ind_test = np.array([i_index_test])

        if cross_valid:
            elasticnetTune = glm.tune_GLM_elastic_net(outpath,model_df = df.iloc[ind_train],y_col_name=y_col_name,x_cols_list = x_cols_list,
                            model_name="GLM_{}_{}".format(protien_isotype,lgtransmethod),alphas=alpharange,l1s=l1range,cv_folds=cvfolds)
        

            elasticnetTune = elasticnetTune.astype(float)
            maxelasticnettune = elasticnetTune.idxmax(axis =1)
            alpha = maxelasticnettune.index[0]
            l1 = maxelasticnettune.iloc[0]

            l1_wt.append(alpha)
            l2_wt.append(l1)
            
        else:
            alpha = l1_l2weights_df.at[i_index_test,"Alpha"]
            l1 =  l1_l2weights_df.at[i_index_test,"L1"]
            trainres,testres = glm.glm_LOO_LR(df, y_col_name, x_cols_list,ind_train,ind_test,
                                        alpha=alpha, L1_wt=l1)

        
        GLM_res_list.append(trainres)
        GLM_test_res_list.append(testres)

    res = glm.glm_loo_plot(fig_path = outpath,model_test_res_list = GLM_test_res_list,model_name = "GLM_{}_{}".format(lgtransmethod,protien_isotype),y_col_name = y_col_name,
                            no_data_pts =df.shape[0],model_train_res_list = GLM_res_list,
                            heatmap_figsize=heatmap_figsize, bar_figsize=bar_figsize,heatmap_annotate_text=heatmap_annotate_text, plot_auc =True,
                            color_specific_yticklabels = common_peptide_highlight)
    if cross_valid:
        alphaL1dict.update({"Alpha":l1_wt,"L1":l2_wt})
        l1l2df = pd.DataFrame(alphaL1dict)
        FileTools.write2Excel(outpath + "/L1_L2_weights_{}_{}.xlsx".format(lgtransmethod,protien_isotype),l1l2df)
    return res




'''
def plot_roc(GLM_test_auc_list,save_path=None, save_auc_txt=None, show_if_none=False,
             title='ROC curve',fig_rows = 8, fig_cols=5,fig_size=(20,15)):
    no_rows = len(GLM_test_auc_list)
    if fig_rows*fig_cols<no_rows:
        print("Subplots value is less few plots might be removed")

    fig, axes = plt.subplots(fig_rows, fig_cols, figsize=fig_size)
    for row in fig_rows:
        for col in fig_cols:



    fig = plt.figure()
    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(fpr, tpr, label='AUC = {:.3f}'.format(auc))
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title(title)
    plt.legend(loc='best')

    PlotTools.savePlt(savePath=save_path, showIfNone=show_if_none)
    if save_auc_txt is not None:
        FileTools.write_list_to_txt([auc], save_auc_txt)

    return fig
'''