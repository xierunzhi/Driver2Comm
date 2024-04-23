import Driver2Comm as dc
import pandas as pd
patient_metadata = pd.read_csv('F:/Driver2Comm/Experiment/data/21_brca/brca_annotation.csv')
a = dc.Driver2Comm(
    minsup=3,
    inputPATH='D:/Driver2Comm/data/cytotalk_input/',
    patient_metadata=patient_metadata,
    outputPATH='D:/Driver2Comm/data',
    cancer_type='Brca',
    visualize=True,
)
a.run()