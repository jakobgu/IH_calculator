
# cwd=%pwd

cwd='/Users/gutzmannjj/IH_analysis'

## importing the relevant packages
import neo
import numpy as np
import matplotlib.pyplot as plt
import openpyxl as ox
import pandas as pd

## defining the function to calculate IH for 5 different current commands from a given .abf file
def IH_calculation(filename): 
    ## loading the data and metadata from the file-to-be-analyzed
    r = neo.AxonIO(filename)
    bl = r.read_block(lazy=False, cascade=True)

    ## defining the relevant measures for partitioning of the signal (e.g. number of runs and number of sweeps per run)
    sweeps=len(bl.segments)
    runs=int(sweeps/9)
    len_of_set=int(sweeps/runs)
    samplingrate=int(bl.segments[0].analogsignals[0].sampling_rate)/1000
    len_of_signal=int(len(bl.segments[0].analogsignals[0]))
    x_axis=np.linspace(1,int(len_of_signal/samplingrate),num=len_of_signal)

    ## setting up new arrays of appropriate sizes to store the extracted and sorted data in
    ave_signal_1=np.zeros((len_of_signal,runs))
    ave_signal_2=np.zeros((len_of_signal,runs))
    ave_signal_3=np.zeros((len_of_signal,runs))
    ave_signal_4=np.zeros((len_of_signal,runs))
    ave_signal_5=np.zeros((len_of_signal,runs))
    ave_signal_6=np.zeros((len_of_signal,runs))
    ave_signal_7=np.zeros((len_of_signal,runs))
    ave_signal_8=np.zeros((len_of_signal,runs))
    ave_signal_9=np.zeros((len_of_signal,runs))
    s1_data=np.zeros((len_of_signal,sweeps))

    ## extracting the data from the neo.core object
    for i in range(sweeps):
        ans = bl.segments[i].analogsignals # every segment (i.e. sweep) contains two signals
        s1 = ans[0] # Signal 1 is membrane voltage
        s2 = ans[1] # Signal 2 is the command current
        s1_data[:,i]=s1.magnitude[:,0] # separating the data from its unit

    ## plotting each intividual signal
        plt.subplot(2,1,1)
        plt.plot(x_axis, s1_data[:,i],'k',linewidth=0.5)
        plt.ylabel("Membrane potential [mV]")

        plt.subplot(2,1,2)
        plt.plot(x_axis, s2.magnitude,'k',linewidth=0.5)
        plt.xlabel("Time [msec]")
        plt.ylabel("command current [pA]")

    ## averaging all signals that were elicited by the same command current (every 10th signal)
    for i in range(runs):
        ave_signal_1[:,i]=s1_data[:,(i*len_of_set)]
        ave_signal_2[:,i]=s1_data[:,(i*len_of_set+1)]
        ave_signal_3[:,i]=s1_data[:,(i*len_of_set+2)]
        ave_signal_4[:,i]=s1_data[:,(i*len_of_set+3)]
        ave_signal_5[:,i]=s1_data[:,(i*len_of_set+4)]
        ave_signal_6[:,i]=s1_data[:,(i*len_of_set+5)]
        ave_signal_7[:,i]=s1_data[:,(i*len_of_set+6)]
        ave_signal_8[:,i]=s1_data[:,(i*len_of_set+7)]
        ave_signal_9[:,i]=s1_data[:,(i*len_of_set+8)]

    ave_signal_1=np.mean(ave_signal_1,1) 
    ave_signal_2=np.mean(ave_signal_2,1) 
    ave_signal_3=np.mean(ave_signal_3,1) 
    ave_signal_4=np.mean(ave_signal_4,1) 
    ave_signal_5=np.mean(ave_signal_5,1) 
    ave_signal_6=np.mean(ave_signal_6,1) 
    ave_signal_7=np.mean(ave_signal_7,1) 
    ave_signal_8=np.mean(ave_signal_8,1) 
    ave_signal_9=np.mean(ave_signal_9,1) 

    ## finding H-current value and location for the first 5 current steps
    H1=np.amin(ave_signal_1[int(280*samplingrate):int(400*samplingrate)])
    H1_loc=((np.argmin(ave_signal_1[int(280*samplingrate):int(400*samplingrate)]))+(280*samplingrate))/50
    real_H1=H1-np.mean((ave_signal_1[int(1250*samplingrate):int(1275*samplingrate)]))
    H2=np.amin(ave_signal_2[int(280*samplingrate):int(400*samplingrate)])
    H2_loc=((np.argmin(ave_signal_2[int(280*samplingrate):int(400*samplingrate)]))+(280*samplingrate))/50
    real_H2=H2-np.mean((ave_signal_2[int(1250*samplingrate):int(1275*samplingrate)]))
    H3=np.amin(ave_signal_3[int(280*samplingrate):int(400*samplingrate)])
    H3_loc=((np.argmin(ave_signal_3[int(280*samplingrate):int(400*samplingrate)]))+(280*samplingrate))/50
    real_H3=H3-np.mean((ave_signal_3[int(1250*samplingrate):int(1275*samplingrate)]))
    H4=np.amin(ave_signal_4[int(280*samplingrate):int(400*samplingrate)])
    H4_loc=((np.argmin(ave_signal_4[int(280*samplingrate):int(400*samplingrate)]))+(280*samplingrate))/50
    real_H4=H4-np.mean((ave_signal_4[int(1250*samplingrate):int(1275*samplingrate)]))
    H5=np.amin(ave_signal_5[int(280*samplingrate):int(400*samplingrate)])
    H5_loc=((np.argmin(ave_signal_5[int(280*samplingrate):int(400*samplingrate)]))+(280*samplingrate))/50
    real_H5=H5-np.mean((ave_signal_5[int(1250*samplingrate):int(1275*samplingrate)]))

    ## initiating new H-current array and populating it with H-current values
    IHval=np.zeros((5,1))
    IHval[0,:]=real_H1
    IHval[1,:]=real_H2
    IHval[2,:]=real_H3
    IHval[3,:]=real_H4
    IHval[4,:]=real_H5


    # plotting the averages on top of all individual signals
    plt.subplot(2,1,1)
    plt.plot(x_axis, ave_signal_1,'r', label=s1.name,linewidth=1)
    plt.plot(H1_loc, H1, marker='*', color='b')
    plt.plot(x_axis, ave_signal_2,'r', label=s1.name,linewidth=1)
    plt.plot(H2_loc, H2, marker='*', color='b')
    plt.plot(x_axis, ave_signal_3,'r', label=s1.name,linewidth=1)
    plt.plot(H3_loc, H3, marker='*', color='b')
    plt.plot(x_axis, ave_signal_4,'r', label=s1.name,linewidth=1)
    plt.plot(H4_loc, H4, marker='*', color='b')
    plt.plot(x_axis, ave_signal_5,'r', label=s1.name,linewidth=1)
    plt.plot(H5_loc, H5, marker='*', color='b')
    plt.plot(x_axis, ave_signal_6,'r', label=s1.name,linewidth=1)
    plt.plot(x_axis, ave_signal_7,'r', label=s1.name,linewidth=1)
    plt.plot(x_axis, ave_signal_8,'r', label=s1.name,linewidth=1)
    plt.plot(x_axis, ave_signal_9,'r', label=s1.name,linewidth=1)

    ## displaying the plot
    plt.show()

    plt.scatter([-200, -150, -100, -50, 0], IHval)
    plt.show()

    return(IHval)


## opening the xcel-file that contains the names of all abf files (row 5) and the whether or not this file contains analyzable IH currents (row 22) 
excel_filename=cwd+'/cell_names.xlsx'
wb = ox.load_workbook(excel_filename, data_only=True)
sheet = wb.get_sheet_by_name('Sheet1')
df = pd.DataFrame(sheet.values)
files=np.array(df.iloc[4:5,1:])
cell_decision=np.array(df.iloc[18:19,1:])
IH_check=np.array(df.iloc[22:23,1:])

## creating a pandas DataFrame that will contain the clalculations made in the IH_calculation function, seperated in columns for each analyzed cell
IH_df=pd.DataFrame([['-200pA:'],['-150pA:'],['-100pA:'],['-50pA:'],['0pA:']])
cell_list=[]
for i in range(len(IH_check[0])):
    if ((IH_check[0][i]==1) & (cell_decision[0][i]==1)):
        cell_list.append(files[0][i])

for ii in range(len(cell_list)):
    IHval=IH_calculation(cwd+'/'+cell_list[ii]+'.abf')
    IH_df.insert(loc=ii+1, column=cell_list[ii], value=IHval)

print(IH_df)