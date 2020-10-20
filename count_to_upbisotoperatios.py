import csv
import numpy as np

mode1 = input("Pb/Pb correction mode: 1.NIST glass or 2.zircon: ")
mode2 = input("Pb/U correction mode: 1.NIST glass or 2.zircon: ")
mode3 = input("Pb/Pb correction mode: 1.NISTSRM612 or 2.NISTSRM610: ")
mode4 = input("Pb/U correction mode: 1.91500 zircon or 2.GJ-1: ")

x = input("number of unknowns: ")

x = int(x)+6


list = ["206Pb/238U", "207Pb/235U", "208Pb/232Th","207Pb/206Pb", "208Pb/206Pb", "232Th/238U"]
if int(mode3)==1:
    ref_nist = [0.254525667358606, 96.4094995340435,  0.548399498523207, 0.90726, 2.164, 3.00]
elif int(mode3)==2:
    ref_nist = [0.257, 97.94378,  0.565373, 0.9096, 2.168612, 3.00]
else:
    print("false")

if int(mode4) == 1:
    ref_zircon = [0.17928,1.86777719,0.054,0.07556,0.855]
elif int(mode4) == 2:
    ref_zircon = [0.09761,0.8093,0.03011,0.06014,0.087354714]
else:
    print("false")

ratio = []
ratio_err = []

len = len(list)

ls = []
num = 0

f = open("output_0909_2.csv",'r').readlines()
for line in f:
    a = line.split(',')
    ls.append(a)
    num = num+1

if int(mode2) == 1:
    for j in range(3):
        seq = int(num//(6+x))
        f_pb_u = []
        err_f_pb_u = []

        for i in range(seq):
            pbu = []
            pbu.append(float(ls[1+x*i][14+2*j]))
            pbu.append(float(ls[2+x*i][14+2*j]))
            pbu.append(float(ls[3+x*i][14+2*j]))
            pbu.append(float(ls[1+x*(i+1)][14+2*j]))
            pbu.append(float(ls[2+x*(i+1)][14+2*j]))
            pbu.append(float(ls[3+x*(i+1)][14+2*j]))

            sum = 0
            sum_sq = 0
            for item in pbu:
                sum = sum+item
                b = sum/6
            for item in pbu:
                sum_sq = sum_sq + (item-b)**2
                #print(sum_sq)
            f_pb_u.append(ref_nist[j]/b)
            err_f_pb_u.append(np.sqrt(sum_sq/5))

        c = []
        d = []
        for i in range(x):
            for k in range(seq):
                c.append(float(ls[1+i+k*x][14+2*j])*f_pb_u[k])
                d.append(c[i+k*x]*np.sqrt((float(ls[1+i+k*x][15+2*j])/float(ls[1+i+k*x][14+2*j]))**2+(err_f_pb_u[k]/f_pb_u[k])**2))
        ratio.append(c)
        ratio_err.append(d)


    if int(mode1) == 1:

        for j in range(3):

            seq = int(num//(6+x))
            f_pb_pb = []
            err_f_pb_pb = []

            for i in range(seq):
                pb = []
                pb.append(float(ls[1+x*i][20+2*j]))
                pb.append(float(ls[2+x*i][20+2*j]))
                pb.append(float(ls[3+x*i][20+2*j]))
                pb.append(float(ls[1+x*(i+1)][20+2*j]))
                pb.append(float(ls[2+x*(i+1)][20+2*j]))
                pb.append(float(ls[3+x*(i+1)][20+2*j]))

                sum = 0
                sum_sq = 0
                for item in pb:
                    sum = sum+item
                    b = sum/6
                for item in pb:
                    sum_sq = sum_sq + (item-b)**2
                    #print(sum_sq)
                f_pb_pb.append(ref_nist[j+3]/b)
                err_f_pb_pb.append(np.sqrt(sum_sq/5))

            c = []
            d = []
            for i in range(x):
                for k in range(seq):
                    c.append(float(ls[1+i+k*x][20+2*j])*f_pb_pb[k])
                    d.append(c[i+k*x]*np.sqrt((float(ls[1+i+k*x][21+2*j])/float(ls[1+i+k*x][20+2*j]))**2+(err_f_pb_pb[k]/f_pb_pb[k])**2))
            ratio.append(c)
            ratio_err.append(d)

    elif int(mode1) == 2:

        for j in range(2):
            seq = int(num//(6+x))
            f_pb_pb = []
            err_f_pb_pb = []

            for i in range(seq):
                pb = []
                pb.append(float(ls[4+x*i][20+2*j]))
                pb.append(float(ls[5+x*i][20+2*j]))
                pb.append(float(ls[6+x*i][20+2*j]))
                pb.append(float(ls[4+x*(i+1)][20+2*j]))
                pb.append(float(ls[5+x*(i+1)][20+2*j]))
                pb.append(float(ls[6+x*(i+1)][20+2*j]))

                sum = 0
                sum_sq = 0
                for item in pb:
                    sum = sum+item
                    b = sum/6
                for item in pb:
                    sum_sq = sum_sq + (item-b)**2
                    #print(sum_sq)
                f_pb_pb.append(ref_zircon[j+3]/b)
                err_f_pb_pb.append(np.sqrt(sum_sq/5))

            c = []
            d = []
            for i in range(x):
                for k in range(seq):
                    c.append(float(ls[1+i+k*x][18+2*j])*f_pb_pb[k])
                    d.append(c[i+k*x]*np.sqrt((float(ls[1+i+k*x][19+2*j])/float(ls[1+i+k*x][18+2*j]))**2+(err_f_pb_pb[k]/f_pb_pb[k])**2))
            ratio.append(c)
            ratio_err.append(d)




    g = open('isotopic_ratio_0909_2.csv', 'w')
    writer = csv.writer(g, lineterminator='\n')



    csvlist = []

    csvlist.append("206Pb/238U")
    csvlist.append("1SD")
    csvlist.append("207Pb/235U")
    csvlist.append("1SD")
    csvlist.append("208Pb/232Th")
    csvlist.append("1SD")
    csvlist.append("207Pb/206U")
    csvlist.append("1SD")
    csvlist.append("208Pb/206Pb")
    csvlist.append("1SD")
    writer.writerow(csvlist)



    for i in range(x):
        csvlist = []
        for j in range(len):
            csvlist.append(ratio[j][i])
            csvlist.append(ratio_err[j][i])
        csvlist.append(ratio[0][i]*ratio[3][i]*137.88)
        csvlist.append(np.sqrt(ratio_err[0][i]*ratio_err[0][i]/ratio[0][i]/ratio[0][i]+ratio_err[3][i]*ratio_err[3][i]/ratio[3][i]/ratio[3][i]))
        writer.writerow(csvlist)


    g.close()


elif int(mode2) == 2:
    for j in range(3):
        seq = int(num//(6+x))
        f_pb_u = []
        err_f_pb_u = []
        ave_pb_u = []

        for i in range(seq):
            pbu = []
            pbu.append(float(ls[4+x*i][14+2*j]))
            pbu.append(float(ls[5+x*i][14+2*j]))
            pbu.append(float(ls[6+x*i][14+2*j]))
            pbu.append(float(ls[4+x*(i+1)][14+2*j]))
            pbu.append(float(ls[5+x*(i+1)][14+2*j]))
            pbu.append(float(ls[6+x*(i+1)][14+2*j]))

            sum = 0
            sum_sq = 0
            for item in pbu:
                sum = sum+item
                b = sum/6
            for item in pbu:
                sum_sq = sum_sq + (item-b)**2

                #print(sum_sq)
            ave_pb_u.append(b)
            f_pb_u.append(ref_zircon[j]/b)
            err_f_pb_u.append(np.sqrt(sum_sq/5))

        c = []
        d = []
#        print(ave_pb_u)
        print(f_pb_u[0], err_f_pb_u[0]/ave_pb_u[0])
#        print(seq)

        for k in range(1):
            for l in range(x+6):
                c.append(float(ls[1+l][14+2*j])*f_pb_u[k])
                d.append(c[l]*np.sqrt((float(ls[1+l][15+2*j])/float(ls[1+l][14+2*j]))**2+(err_f_pb_u[k]/ave_pb_u[k])**2))
            #    c.append(float(ls[1+i+k*x][14+2*j])*f_pb_u[k])
            #    d.append(c[i+k*x]*np.sqrt((float(ls[1+i+k*x][15+2*j])/float(ls[1+i+k*x][14+2*j]))**2+(err_f_pb_u[k]/ave_pb_u[k])**2))

                #print((float(ls[1+i+k*x][15+2*j])/float(ls[1+i+k*x][14+2*j])))
#                print((err_f_pb_u[k]), ave_pb_u[k])
#                print((err_f_pb_u[k]/ave_pb_u[k]))
        ratio.append(c)
        ratio_err.append(d)
#    print(ratio)
#    print(c)
    if int(mode1) == 1:

        for j in range(3):

            seq = int(num//(6+x))
            f_pb_pb = []
            err_f_pb_pb = []
            ave_pb_pb = []

            for i in range(seq):
                pb = []
                pb.append(float(ls[1+x*i][20+2*j]))
                pb.append(float(ls[2+x*i][20+2*j]))
                pb.append(float(ls[3+x*i][20+2*j]))
                pb.append(float(ls[1+x*(i+1)][20+2*j]))
                pb.append(float(ls[2+x*(i+1)][20+2*j]))
                pb.append(float(ls[3+x*(i+1)][20+2*j]))

                sum = 0
                sum_sq = 0
                for item in pb:
                    sum = sum+item
                    b = sum/6
                for item in pb:
                    sum_sq = sum_sq + (item-b)**2
                    #print(sum_sq)
                f_pb_pb.append(ref_nist[j+3]/b)
                err_f_pb_pb.append(np.sqrt(sum_sq/5))
                ave_pb_pb.append(b)
#                print(ave_pb_pb)
                print(f_pb_pb[0], err_f_pb_pb[0]/ave_pb_pb[0])


            c = []
            d = []
#            print(ls)
            for i in range(x):
#                for k in range(1):
#                    c.append(float(ls[1+i+k*x][20+2*j])*f_pb_pb[k])
#                    d.append(c[i+k*x]*np.sqrt((float(ls[1+i+k*x][21+2*j])/float(ls[1+i+k*x][20+2*j]))**2+(err_f_pb_pb[k]/ave_pb_pb[k])**2))
                for k in range(1):
                    for l in range(x+6):
                        #print(x)
                        c.append(float(ls[1+l][20+2*j])*f_pb_pb[k])
                        d.append(c[l]*np.sqrt((float(ls[1+l][21+2*j])/float(ls[1+l][20+2*j]))**2+(err_f_pb_pb[k]/ave_pb_pb[k])**2))

            #for l in range(6):
                #c.append(float(ls[x+l][20+2*j])*f_pb_pb[0])
                #d.append(c[x+1+l]*np.sqrt((float(ls[x+l][21+2*j])/float(ls[x+l][20+2*j]))**2+(err_f_pb_pb[0]/ave_pb_pb[0])**2))
            ratio.append(c)
            ratio_err.append(d)


    elif int(mode1) == 2:

        for j in range(2):
            seq = int(num//(6+x))
            f_pb_pb = []
            err_f_pb_pb = []

            for i in range(seq):
                pb = []
                pb.append(float(ls[4+x*i][20+2*j]))
                pb.append(float(ls[5+x*i][20+2*j]))
                pb.append(float(ls[6+x*i][20+2*j]))
                pb.append(float(ls[4+x*(i+1)][20+2*j]))
                pb.append(float(ls[5+x*(i+1)][20+2*j]))
                pb.append(float(ls[6+x*(i+1)][20+2*j]))

                sum = 0
                sum_sq = 0
                for item in pb:
                    sum = sum+item
                    b = sum/6
                for item in pb:
                    sum_sq = sum_sq + (item-b)**2
                    #print(sum_sq)
                f_pb_pb.append(ref_zircon[j+3]/b)
                err_f_pb_pb.append(np.sqrt(sum_sq/5))

            c = []
            d = []
            for i in range(x):
                for k in range(seq):
                    c.append(float(ls[1+i+k*x][18+2*j])*f_pb_pb[k])
                    d.append(c[i+k*x]*np.sqrt((float(ls[1+i+k*x][19+2*j])/float(ls[1+i+k*x][18+2*j]))**2+(err_f_pb_pb[k]/f_pb_pb[k])**2))

            ratio.append(c)
            ratio_err.append(d)



    g = open('isotopic_ratio_0909_2.csv', 'w')
    writer = csv.writer(g, lineterminator='\n')


    csvlist = []

    csvlist.append("206Pb/238U")
    csvlist.append("1SD")
    csvlist.append("207Pb/235U")
    csvlist.append("1SD")
    csvlist.append("208Pb/232Th")
    csvlist.append("1SD")
    csvlist.append("207Pb/206Pb")
    csvlist.append("1SD")
    csvlist.append("208Pb/206Pb")
    csvlist.append("1SD")
    csvlist.append("Th/U")
    csvlist.append("1SD")
    csvlist.append("207Pb/235U(206Pb/238U*207Pb/206Pb*137.88)")
    csvlist.append("1SD")
    csvlist.append("206Pb/238U age / Ma")
    csvlist.append("1SD")
    csvlist.append("207Pb/235U age / Ma")
    csvlist.append("1SD")
    csvlist.append("208Pb/232Th age / Ma")
    csvlist.append("1SD")

    writer.writerow(csvlist)

    for i in range(x+6):
        csvlist = []
        for j in range(len):
            csvlist.append(ratio[j][i])
            csvlist.append(ratio_err[j][i])
        csvlist.append(ratio[0][i]*ratio[3][i]*137.88)
        csvlist.append(ratio[0][i]*ratio[3][i]*137.88*np.sqrt(ratio_err[0][i]*ratio_err[0][i]/ratio[0][i]/ratio[0][i]+ratio_err[3][i]*ratio_err[3][i]/ratio[3][i]/ratio[3][i]))
        csvlist.append(np.log(1+ratio[0][i])/1.55125e-10/1e+6)
        csvlist.append(np.log(1+ratio[0][i]+ratio_err[0][i])/1.55125e-10/1e+6 - np.log(1+ratio[0][i])/1.55125e-10/1e+6)
        csvlist.append(np.log(1+ratio[0][i]*ratio[3][i]*137.88)/0.00000000098485/1e+6)
        csvlist.append(np.log(1+ratio[0][i]*ratio[3][i]*137.88+ratio[0][i]*ratio[3][i]*137.88*np.sqrt(ratio_err[0][i]*ratio_err[0][i]/ratio[0][i]/ratio[0][i]+ratio_err[3][i]*ratio_err[3][i]/ratio[3][i]/ratio[3][i]))/0.00000000098485/1e+6 - np.log(1+ratio[0][i]*ratio[3][i]*137.88)/0.00000000098485/1e+6)
        csvlist.append(np.log(1+ratio[2][i])/0.000000000049475/1e+6)
        csvlist.append(np.log(1+ratio[2][i]+ratio_err[2][i])/0.000000000049475/1e+6 - np.log(1+ratio[2][i])/0.000000000049475/1e+6)


        writer.writerow(csvlist)

    g.close()



else:
    print("false")
