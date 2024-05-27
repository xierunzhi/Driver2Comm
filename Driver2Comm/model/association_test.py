"""
association test part of Driver2Comm
association test : fisher exact test
"""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import statsmodels.stats.multitest
class AssociationTest(object):
    """
    Association part of Driver2Comm
    """
    def __init__(self, internal_factor : pd.DataFrame,external_factor : pd.DataFrame , threshold = 0.05):
        """

        :param internal_factor:     internal matrix to be analysis
        :param external_factor:
        :param threshold:
        """
        self.internal = internal_factor
        self.external = external_factor
        self.threshold = threshold
        self.external_combination = None


    def my_fisher_exact_test(self,external:pd.Series,internal:pd.Series):
        """

        :param external:
        :param internal:
        :return:
        """
        a = sum((external+internal)==0)
        d = sum((external + internal) == 2)
        b = sum((internal-external)==1)
        c = sum((internal-external)==-1)
        _, pvalue = fisher_exact([[a, b], [c, d]],alternative='greater')
        return pvalue

    def association_test(self):
        """
        Association test version 2 ,which give a whole p values list and perform a multi testing correction (Benjamini/Hochberg)

        :return:
        """
        nfp = self.external.shape[0]
        ndriver = self.internal.shape[0]
        ret = list()
        for i in range(ndriver):
            pvalue_list = []
            internal_vec = self.internal.iloc[i,:]
            for idx in range(nfp):
                external_vec = self.external.iloc[idx,:]
                pvalue_list.append(self.my_fisher_exact_test(external_vec,internal_vec))
            # perform a multiTest correction
            pass_list,adjust_pvalue,_,_ = statsmodels.stats.multitest.multipletests(pvalue_list,alpha= self.threshold,method = 'fdr_bh')
            pass_idx = np.where(pass_list)
            self.external_combination = self.external.take(pass_idx[0])
            pvalue_list = np.take(adjust_pvalue,pass_idx[0])
            pvalue_sorted_idx = np.argsort(pvalue_list)
            pvalue_list = list(np.take(pvalue_list,pvalue_sorted_idx))
            tested_fp_list = list(np.take(pass_idx[0],pvalue_sorted_idx))
            print("Number of {}-associated CCC signature: {}".format(self.internal.index[i],len(pvalue_list)))
            ret.append({'internal':self.internal.index[i],'external matrix':self.external_combination.T
                           ,'idx of passed FP':tested_fp_list,'internal vec':internal_vec
                        ,'pvalue of passed FP':pvalue_list})
        return ret



