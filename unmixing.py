import numpy as np
import tifffile as tiff

filename = 'sox10BAC_AviRangap_18s3_subset_Subset1-16.lsm'
image = tiff.imread(filename)
ind = image.shape
# channel = ind[1]

# find the individual spectra: autofluoresence, mCherry, CFP
mat_g = image[0,0,:,363:379,959:965].mean(axis=(1,2))
mat_r1 = image[0,0,:,38:45,503:508].mean(axis=(1,2))
mat_r2 = image[0,0,:,57:68,528:535].mean(axis=(1,2))
mat_r = (mat_r1+mat_r2)/2
mat_r[0:15]=0
mat_b = image[0,0,:,127:177,773:785].mean(axis=(1,2))


indvars = np.column_stack((mat_g,mat_r,mat_b))
jacobian = np.zeros(shape=(indvars.shape[1],indvars.shape[1]))
jvector = np.zeros(shape=(indvars.shape[1],1))
coefcube = np.zeros(shape=(indvars.shape[1], ind[3], ind[4] ))

# get the jacobian from indvars
for i in range(0, indvars.shape[1]):
    for j in range(0, i + 1 ):
        jacobian[i,j] = np.multiply(indvars.T[i], indvars.T[j]).sum()
        if i != j:
            jacobian[j,i] = jacobian[i,j]

for zp in range(0,ind[1]):
    for i in range(0,ind[3]):
        for j in range(0,ind[4]):
            data = image[0,zp,:,i,j]
            jvector = np.multiply(data, indvars.T).sum(axis=1)
            coefcube[:,i,j] = np.linalg.solve(jacobian, jvector)
    coefcube[coefcube<0]=0
    coefcube[coefcube>1]=1

    temp = image[0,zp,:,:,:].sum(axis = (0))
    test = np.multiply(temp,coefcube)
    # save tiff files
    for ch in range(0,indvars.shape[1]):
        newfilename = filename[0:-4] + '-z' + (str(zp)) + '-ch'+(str(ch)) + '.tif'
        tiff.imsave(newfilename, np.uint16(test[ch]))
