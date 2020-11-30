import numpy as np
import pexpect
#Run in folder mono of the OPCD data


path1="/Users/lennertthormaehlen/Documents/Uni/Heidelberg/Master/Astro-approach/OP/OPCD_3.3/mono/"
path2="/Users/lennertthormaehlen/Documents/Uni/Promotion/axionflux/SolarAxionFlux/data/opacity_tables/OP/"
ite=150
jne=66
Z=[1,2,6,7,8,10,11,12,13,14,16,18,20,24,25,26,28]
element_name = ["H","He","C","N","O","Ne","Na","Mg","Al","Si","S","Ar","Ca","Cr","Mn","Fe","Ni"]
for iz in range(0,17):
    label='temp'+'m'+str(int(Z[iz])).zfill(2)+'.'+str(ite).zfill(3)+'.'+str(jne)
    try:
        header=np.genfromtxt(path1+label,max_rows=1,skip_header=2)
        headerlength=int(header[1]-header[0])+4
        s=np.genfromtxt(path1+label,skip_header=headerlength)
        s[:,0]=np.log(s[:,0])
        filename="opacity_table_"+element_name[iz]+"_"+str(int(ite))+"_"+str(int(jne))+".dat"
        np.savetxt(path2+filename,s,fmt=["%9.5e","%9.2e"],header="OPACITY Project tables",delimiter="")
    except OSError:
        proc=pexpect.spawn('./monop.out')
        proc.sendline('m'+str(int(Z[iz])).zfill(2))
        proc.sendline(str(ite).zfill(3))
        proc.sendline(str(jne))
        label='temp'+'m'+str(int(Z[iz])).zfill(2)+'.'+str(ite).zfill(3)+'.'+str(jne)
        proc.sendline(label)
        proc.expect('WRITTEN')
        print(label)
        header=np.genfromtxt(path1+label,max_rows=1,skip_header=2)
        headerlength=int(header[1]-header[0])+4
        s=np.genfromtxt(path1+label,skip_header=headerlength)
        s[:,0]=np.log(s[:,0])
        filename="opacity_table_"+element_name[iz]+"_"+str(int(ite))+"_"+str(int(jne))+".dat"
        np.savetxt(path2+filename,s,fmt=["%9.5e","%9.2e"],header="OPACITY Project tables",delimiter="")
