import numpy as np
import tifffile as tiff
import time

filename = '/Users/rose/scripts/unmixing_tests/sox10BAC_AviRangap_18s3_subset_Subset1.lsm'

image = tiff.imread(filename)
#image = image.asarray()
ind = image.shape

mat_g = np.array([23.3529411764706,66.0084033613445,119.193277310924,3.17647058823529,99.7226890756303,22.4201680672269,1464.19327731092,1639.48739495798,1655.83193277311,1964.35294117647,3094.43697478992,3203.25210084034,3423.16806722689,3034.48739495798,2524.19327731092,1756.52941176471,64.1596638655462,111.478991596639,1737.19327731092,1568.30252100840,2378.53781512605,2414.48739495798,2757.32773109244,2078.40336134454,1908.86554621849,1282.04201680672,1101.15966386555,757.571428571429,751.193277310924	,506.915966386555	,529.033613445378	,479.445378151261])
mat_r = np.array([0	,0	,0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	398.322916666667,	19.6197916666667,	40.4218750000000,	531.817708333333,	1057.75520833333,	1839.18229166667,	2562.25520833333,	2810.59375000000,	2336.18229166667,	1859.68229166667,	1349.33854166667,	1072.64062500000,	1064.69270833333,	832.760416666667,	737.255208333333,	581.734375000000,	327.791666666667])
mat_b = np.array([40.7692307692308,	75.3499245852187,	128.188536953243,	3.98793363499246,	97.5671191553545,	31.6983408748115,	1193.49321266968,	1646.87631975867,	1439.75263951735,	1399.86123680241,	1574.42081447964,	1274.33484162896,	1098.04374057315,	810.850678733032,	648.865761689291,	408.399698340875,	21.1628959276018,	36.6530920060332,	321.734539969834,	238.428355957768,	193.601809954751,	246.591251885370,	247.185520361991,	161.921568627451,	109.594268476621,	91.5158371040724,	85.2126696832579,	105.184012066365,	78.5837104072398,	75.8174962292609,	55.3619909502262,	54.2941176470588])


# find the individual spectra: autofluoresence, mCherry, CFP the pixel positions are from the 7 z-stacks subset
# mat_g =np.mean(numpy.float64(image[0,0,:,362:378,958:964]) ,axis=(1,2)) #image[0,0,:,363:379,959:965].mean(axis=(1,2))
# mat_r1 = np.mean(image[0,0,:,38:45,503:508],axis=(1,2))
# mat_r2 = np.mean(image[0,0,:,57:68,528:535],axis=(1,2))
# mat_r = (mat_r1+mat_r2)/2
# mat_r[0:15]=0
# mat_b = image[0,0,:,127:177,773:785].mean(axis=(1,2))



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
    t1 = time.perf_counter()
    for i in range(0,ind[3]):
        for j in range(0,ind[4]):
            data = np.copy(image[0,zp,:,i,j])
            jvector = np.multiply(data, indvars.T).sum(axis=1)
            coefcube[:,i,j] = np.linalg.solve(jacobian, jvector)
    coefcube[coefcube<0]=0
    coefcube[coefcube>1]=1
    t2 = time.perf_counter()
    print("Elapsed Time: " + str(t2 - t1))

    temp = image[0,zp,:,:,:].sum(axis = (0))
    test = np.multiply(temp,coefcube)
    unmix = (test/test.max())*(2**16)
    # save tiff files
    for ch in range(0,indvars.shape[1]):
        newfilename = filename[0:-4] + '-z' + (str(zp)) + '-ch'+(str(ch)) + '.tif'
        tiff.imsave(newfilename, np.uint16(unmix[ch]))
