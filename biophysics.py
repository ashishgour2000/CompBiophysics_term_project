#Sl No. 23
#Souvik Das 22AT91R09
#Yash Jain 19CH10078
#Ashish Gour 18CS30008

#CompBiophysics term project code

#Readme
# 1. This code takes .pdb file as input and generate a 1D binary array as output
# 2. Have python3 and numpy installed in your machine
# 3. Save pdb file and the code in the same folder/directory
# 4. Update the 'pdbfile' variable to the name of the pdb file

#import libraries
import numpy as np

#pdb file to be read
pdbfile = '1a8g.pdb'

#reading pdb file
list2 = []
try:
  with open(pdbfile) as pdbfile:
    for line in pdbfile:
        if line[:4] == 'ATOM':
            list1 = []
            list1 = line.split(' ')
            while("" in list1):
              list1.remove("")
            list2.append(list1)
except:
  print("No such pdb file exists in the same folder.")
            
#Remove H atoms from both protein and ligand structure
#Storing ligand structure in lig_list and protein structure in protein_list
lig_list = []
protein_list = []
h_count = 0
for every in list2:
  if every[2][0] =='H':
    h_count = h_count+1
  elif every[3] == "LIG":
    lig_list.append(every)
  else:
    protein_list.append(every)

#calculating euclidean distance between each protein-ligand atom pairs
distance_array = np.zeros((len(lig_list), len(protein_list)))
count=0
for i in range(0,len(lig_list)):
  for j in range(0,len(protein_list)):
    p1 = np.array([float(lig_list[i][6]), float(lig_list[i][7]), float(lig_list[i][8])])
    p2 = np.array([float(protein_list[j][6]), float(protein_list[j][7]), float(protein_list[j][8])])
    squared_dist = np.sum((p1-p2)**2, axis=0)
    dist = np.sqrt(squared_dist)
    distance_array[i][j] = dist

#take threshold as input
threshold = input("What's the threshold distance: ")
threshold = float(threshold)

#generating protein fingerprint
final_array = np.zeros((1, len(protein_list)))
for j in range(0,len(protein_list)):
  flag=0
  for i in range(0,len(lig_list)):
    if distance_array[i][j] <= threshold:
      flag=1
  final_array[0][j]=flag

#output
print("No of binding sites in protein: ",int(np.sum(final_array)))
print("Protein fingerprint")
#print(final_array[0])
for i in range(0,len(protein_list)):
    print(final_array[0][i])