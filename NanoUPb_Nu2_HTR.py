import pandas as pd
import csv
import glob
import os

datalist = glob.glob(os.getcwd()+'/*.csv')
print(datalist)
x = input("timeslice:")
timeslice = int(x)

integration_time = 20*timeslice

for item in datalist:
#    df = pd.read_csv("output_026_20us_Plesov1.csv")
    df = pd.read_csv(item)
    csvname = item.split("/")[-1]
    g = open("integrated/"+str(integration_time)+"us_"+csvname, 'w')
    writer = csv.writer(g, lineterminator='\n')
    csvlist = []
    csvlist.append("Timestamp")
    csvlist.append("U(238)")
    csvlist.append("Pb(206)")
    writer.writerow(csvlist)
    for i in range(int(len(df["Step"])/timeslice)):
        csvlist = []
        sum238 = 0
        sum206 = 0
        for j in range(timeslice):
            sum238 += df["U(238)"][i*timeslice+j]
            sum206 += df["Pb(206)"][i*timeslice+j]
        csvlist.append(df["Timestamp"][i*timeslice+j])
        csvlist.append(sum238/timeslice)
        csvlist.append(sum206/timeslice)
        writer.writerow(csvlist)
    g.close()

    #    print(sum)
