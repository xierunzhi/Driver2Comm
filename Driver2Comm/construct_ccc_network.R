library(CytoTalk)
library(reticulate)
use_condaenv('driver2comm',required = TRUE)
celltype_list = c('Cd8+Tcells','Macrophages','Tumor')
ct_combn = combn(celltype_list,2)
input_dir = 'F:/pancancer/21_NG_breast_cancer/data/resample_input/'
patient_list <- list.dirs(input_dir,recursive = F,full.names = F)
dir_out<- 'D:/test/cytotalk_output'
for(i in 1:length(patient_list)){
  patient_id = patient_list[i]
  input_path = paste(input_dir,patient_id,sep = '/')
  print(paste(patient_id, 'is processing ......'))
  outputPATH = paste(dir_out,'/',patient_id,sep = '')
  for(j in 1:(length(ct_combn)/2)){
    type_a = ct_combn[j]
    type_b = ct_combn[j+1]
    path_AB = paste(outputPATH,'/',type_a,'-',type_b,sep = '')
    if(!dir.exists(path_AB)){
      dir.create(path_AB,recursive = T)
    }
    lst_scrna <- read_matrix_folder(input_path)
    print(paste(type_a,'and ',type_b,'are processing ...'))
    run_cytotalk(lst_scrna,
                             type_a,
                             type_b,
                             cutoff_a = 0.05,
                             cutoff_b = 0.05,
                             dir_out = path_AB)
  }
}

#need to check cutoff
