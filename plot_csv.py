/******************************** plot_csv.py **************************************/
#!/usr/bin/env python3
import sys, os
import csv
import matplotlib.pyplot as plt


if len(sys.argv)<2:
print("Usage: python3 plot_csv.py output.csv"); sys.exit(1)
path=sys.argv[1]
os.makedirs('plots', exist_ok=True)
cols=[]
with open(path,'r',newline='') as f:
rdr=csv.reader(f)
header=next(rdr)
N=len(header)
for _ in range(N): cols.append([])
for row in rdr:
if len(row)!=N: continue
for i in range(N): cols[i].append(float(row[i]))


names=header
T=cols[0]
for i in range(1,len(names)):
plt.figure()
plt.plot(T, cols[i])
plt.xlabel('t [s]')
plt.ylabel(names[i])
plt.title(f"{names[i]} vs time")
plt.grid(True)
out=f"plots/{i:02d}_{names[i].replace(' ','_').replace('/','-')}.png"
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print("Saved plots in ./plots/")
