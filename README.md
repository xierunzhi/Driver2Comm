# Driver2Comm

A general computationalframework for associating cancer cell-intrinsic driver genes with -extrinsic cell-cell communication

## Installation

### Required Packages

The following software must be installed on your machine:

Python: tested with version 3.10

R :  tested with R 4.0.3

#### Python dependencies

* numpy (tested v1.26.2)
* pandas (tested v2.1.4)
* statsmodels (tested v0.14.1)
* scipy (tested v1.11.4)
* networkx (tested v3.2.1)
* matplotlib (tested v)
* pcst-fast (tested v1.0.7)

#### R dependencies

* CytoTalk(tested v0.99.0)

 #### Preparing the virtual environment

1. Create a new conda environment using the environment.yml file or the requirements.txt file with one of the following commands:

   ```bash
   conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
   conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
   conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda
   conda env create -f environment.yml
   # or
   conda create --name Driver2Comm --file requirements.txt
   
   ```

2. Install R and the Cytotalk package

   ``` bash
   conda install R
   R
   ```

   ```r
   > install.packages("devtools")
   > devtools::install_github("tanlabcode/CytoTalk")
   ```

##### Notation:

* For detail information of Cytotalk installation, please refer to [tanlabcode/Cytotalk](https://github.com/tanlabcode/CytoTalk)  

## Tutorial

### Preparation

According to Cytotalk's requirement, the input data folder should be agreed with following struture.

```txt
── Cytotalk_input
   ├─ patient1
   		├─ scRNAseq_Tumorcells.csv
       	├─ scRNAseq_BCells.csv
       	├─ scRNAseq_EndothelialCells.csv
       	├─ scRNAseq_Fibroblasts.csv
       	├─ scRNAseq_LuminalEpithelialCells.csv
       	├─ scRNAseq_Macrophages.csv
       	└─ scRNAseq_CD8+TCells.csv
   ├─ patient2
   		├─ scRNAseq_Tumorcells.csv
       	├─ scRNAseq_BCells.csv
       	├─ scRNAseq_EndothelialCells.csv
       	├─ scRNAseq_Fibroblasts.csv
       	├─ scRNAseq_Macrophages.csv
       	└─ scRNAseq_CD8+TCells.csv
   ├─ patient3
   		├─ ...
   └─ patient4
   		├─ ...
```

### Usage

You can apply Driver2Comm in the following step.

- Step0: Construction of cell-cell communication network
- Step1: Identification of cancer-driver assoication communicaiton patterns
- Step2: Downstream analysis of Driver2Com

#### Step0: construction of cell-cell communication network

```bash
Rscript construction_ccc_network.R
```

#### Step1: Identification of cancer-driver assoicatied communicaiton patterns

After construction of CCC network by Cytotalk, we started to identify the cancer-driver associated communicaiton patterns. 
The rest analysis of Driver2Comm is fully based on python implement. Based on the excellent flexibility and interactivity of jupyter notebook, we recommend you to use notebook for all of Driver2Comm's remaining analyses.

```python
d2c_instance = Driver2Comm(
					minsup = ,
					num_celltype=3,
					mode = 'driver',
					outputPATH ='./result',
                 	subtype_annotation = ,
                 	visualize = False,
                 	cancer_type = 'Lung'
    """
    parameters:
    __________________

        :param minsup: Hyperparameter, the minimal suppport of gSpan algorithm
        :param inputPATH: the input file path of Cytotalk input
        :param patient_metadata: the metadata of patients, must contain: driver gene of each patient
        :param num_celltype: the number of cell type in Multi-Cell-Type-Communication MCTC networks
        :param mode: cancer subtype or driver
        :param outputPATH: the output file path
        :param subtype_annotation: the subtype_annotation file path
        :param visualize: if user visualize the cancer driver associated CCC signature or not
        :param cancer_type: type of cancer
    __________________
        
    """
)
# run Driver2Comm to identify cancer driver associated CCC signatures.
d2c_instance.run()
```

#### Step2: Downstream analysis of Driver2Comm

We offer several ways to visualize and evaluate the cancer driver associated CCC signatures identified by Driver2Comm.

```python
d2c_instance.display_associated_FP()
# display all cancer driver associated CCC signatures on screen.

```

## Maintainer

Runzhi Xie (rzxie@stu.xidian.edu.cn)

## Acknowledgement

The gSpan implementation in Driver2Comm was based in part on source code of [gspan-mining](https://github.com/betterenvi/gSpan) package.

## Citation