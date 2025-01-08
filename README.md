# Driver2Comm

A general computational framework for associating cancer cell-intrinsic driver genes with -extrinsic cell-cell communication

## New

1. Driver2Comm is now available for **any non-tumor cell type of interest** in the TME. (May 27th, 2024)

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

#### R dependencies

* CytoTalk(tested v0.99.0)

 #### Preparing the virtual environment

1. Create a new conda environment using the environment.yml file or the requirements.txt file with one of the following commands:

   ```bash
   conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
   conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
   conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda
   conda env create -f environment.yml
   
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
3. Install Driver2Comm
   install Driver2Comm using pip
   ``` bash
      pip install Driver2Comm -U
   ```
   
##### Notation:

* For detail information of Cytotalk installation, please refer to [tanlabcode/Cytotalk](https://github.com/tanlabcode/CytoTalk)  


## Tutorial

### Preparation

the input of Driver2Comm include two part:

1. intrinsic factors: the driver information of patients in this study 
2. Extrinsic factors: the MCTC networks of patients in this study 

#### Intrinsic factor
Driver2Comm needs patient driver information as input，where column1 is patients' id, the remaining columns are driver involving in your study. 1 means the driver variated in this patient. 0 means not.

|      | Patient_id | Driver1 | Driver2 |
| ---- | ---------- | ------ | ------ |
| 0    |      patient1     | 1       |    0    |
| 1    |       patient2       |   0    |   0    |
| 2    |       patient3       |   1   |   0    |
| 3    |      patient41      | 1       |    1    |


#### Extrinsic factor
If you use CytoTalk as your MCTC network construction algorthm, the input data folder should be agreed with following structure.

```txt
── CytoTalk_output
   ├─ patient1
   		├─ Celltype1-Celltype2
   			├─ FinalNetwork.txt
       	├─ Celltype1-Celltype3
       		├─ FinalNetwork.txt
       	└─ Celltype2-Celltype3
       		├─ FinalNetwork.txt
   ├─ patient2
   		├─ Celltype1-Celltype2
   			├─ FinalNetwork.txt
       	├─ Celltype1-Celltype3
       		├─ FinalNetwork.txt
       	└─ Celltype2-Celltype3
       		├─ FinalNetwork.txt
   ├─ patient3
   		├─ ...
   └─ patient4
   		├─ ...
   		
```

```python
from Driver2Comm import formulate_c2c_network
inputPATH = './data/CytoTalk_output' 
outputPATH = './data/Driver2Comm_input'
patient_list = ['patient1','patient2',...,'patientn']
formulate_c2c_network(inputPATH,outputPATH,patient_list)
```

**Warning !!!**
**1. the order of `patient_list` should exactly coincide with  patient driver information matrix**

Then in your output directory, you will find a `c2c_network.data` , and the  `c2c_network.data` is in following format.

Driver2Comm can extend to other CCC inferences method, in this case, you can also directly formulate you cell-cell communication as the following format.

```txt
t # 1 			# t # 1 means the MCTC network of patient 1
v 0 FERMT3__Macrophages # vertex data was formualed as 'v vertex_id vertex_label'
v 1 MITF__Macrophages
e 0 1 0					# edge was formualed as 'e vertex_id_a vertex_id_b edge label'
...
t # n
v 1 SOCS1__Cd8+Tcells
v 1 CD52__Cd8+Tcells
e 0 1 0
t # -1			# t # -1 represent as the end of this file
```
**2. the label of celltype should not contain any space**
for example:
"CD8+ T" **×** 
"CD8+_T" **✓**
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
import Driver2Comm as dc
c2c_network = dc.read_c2c_network(os.path.join('./data/c2c_network.data'))
model = dc.Driver2Comm(
					c2c_network = c2c_network,
                    patient_driver_info=patient_driver_info,
                    minsup= 20,
                    association_test_threshold=0.05,
                    outputPATH= './data/Driver2Comm_output',
                    cancer_type='Brca'
    """
    parameters:
    __________________
        :param c2c_network: the formulated cell-cell communication network data
        :param patient_driver_info: the patient metadata formulated as the form shown in previous section "intrinsic factor"
        :param minsup: Hyperparameter, the minimal suppport of frequent subnetwork extraction algorithm
        :param association_test_threshold: Hyperparameter, the threshold of association testing
        :param outputPATH: the path where output files place
        :param cancer_type: type of cancer
    __________________
        
    """
)
# run Driver2Comm to identify cancer driver associated CCC signatures.
model.run()
```

#### Step2: Downstream analysis of Driver2Comm

We offer several ways to visualize and evaluate the cancer driver associated CCC signatures identified by Driver2Comm.

```python
# display all cancer driver associated CCC signatures on screen.
model.display_associated_FP()
# output association testing result of Driver2Comm as a Pandas Dataframe
associated_mat = model.output_assocaited_mat()

```
### Downstream evaluation

#### Intrinsic-Extrinsic (IE) pathway identification
we provide two function to get IE pathway
##### to get the IE pathway of top k driver-associated CCC signature
```python
import Driver2Comm as dc
interested_driver = 'ESR1__Tumor' # your interested driver
patient_path = './xxx' # the path containing cytotalk result of interested patient
k = 10 # number of interested signatures
top_k_ie_pathway_mat = dc.generate_top_k_ie_pathway(model.association_test_ret,model.frequent_subgraphs,interested_driver, k=k,patient_path=patient_path)
```
##### to get the IE pathway of specific interested gene in CCC signature
```python
import Driver2Comm as dc
interested_driver = 'ESR1__Tumor' # your interested driver
patient_path = './xxx' # the path containing cytotalk result of interested patient
interested_gene = 'APOE_Macrophages' # your interested gene
pcsf_network, original_network = dc.read_graph_from_cytotalk_output(patient_path)
ie_ret = dc.get_ie_pathway(interested_driver,interested_gene,pcsf_network,original_network)
```

## Maintainer

Runzhi Xie (rzxie@stu.xidian.edu.cn)

## Acknowledgement

The gSpan implementation in Driver2Comm was based in part on source code of [gspan-mining](https://github.com/betterenvi/gSpan) package.

## Citation