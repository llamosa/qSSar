from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors, MACCSkeys
from rdkit.Chem import Draw
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
#from structure import molecule
import glob as glob
import numpy as np
#from scipy.spatial.distance import cdist
#%matplotlib inline
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import os
#from matplotlib.pyplot import imshow
#from PIL import Image
from IPython.display import Image, display

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


smiles = []
names = []

supp = Chem.SmilesMolSupplier('/data/datasets/mllamosa/SMILES_price_SAR/ZincCompounds_WaitOK.smi')
mlist = [m for m in supp]
finger_dimension = 167
step = int(len(supp) / 10)

#for struct_file in glob.glob('C:/Users/mllamosa.PROSTATECENTRE/Dropbox/2017/UBC_Prostate/Tox/deepchem/xyz/home/mllamosa/data/tox21_deepchem/xyz/*xyz'):
#for struct_file in glob.glob('/home/fer19x/data/data61/data/10k_diastereomers_B3LYP_G4MP2/x*'):
        #print(struct_file)
for j in range(100):

    i=0
    molecules = []
    fps = []
    id = []
    print("Split",j)
    for m in mlist[j*step:(j+1)*step]:
        #print(smi)
        try:
            #m = Chem.MolFromSmiles(smi)
            Chem.Kekulize(m)
            molecules.append(m)
            fps.append(np.array([fp for fp in MACCSkeys.GenMACCSKeys(m).ToBitString()]))
            #fps.append(np.array([fp for fp in AllChem.GetMorganFingerprintAsBitVect(m, 2).ToBitString()]))
            id.append(i)
        except:
            fps.append(np.repeat(0,finger_dimension))
            id.append(i)
            print(i,m)
        i = i + 1

    header = ["mol"]

    for i in range(finger_dimension):
        header.append("fps"+ str(i))


    fps = np.array(fps).reshape(len(fps),finger_dimension)
    id = np.array(id)
    id =id.reshape(len(fps),1)
    data = np.hstack((id,fps))
    header = np.array(header).reshape(1,len(header))
    data_header = np.vstack((header,data))
    np.savetxt("/data/datasets/mllamosa/SMILES_price_SAR/ZincCompounds_WaitOK_maccs.tab_" + str(step), data_header, delimiter=",", fmt="%s")
