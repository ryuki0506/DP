import numpy as np
import matplotlib.pyplot as plt

data=np.loadtxt("result.txt",delimiter=",")
data_row=data.shape[0]
data_column=data.shape[1]

ps=[data[i+1][0] for i in range(data_row-1)]
lengths=data[0][1:]
conGRs=[data[i+1][data_column-1] for i in range(data_row-1)]

pmax=np.max(ps)
pmin=np.min(ps)
steps=len(ps)

def outputGR_vs_len(plistsize):
  plistsize=plistsize-1
  fig, ax=plt.subplots()
  for i in range(plistsize):
    pi=(i*steps)//plistsize
    ax.plot(lengths,data[1+pi][1:],label=f"$p=${ps[pi]}")
  
  ax.plot(lengths,data[data_row-1][1:],label=f"$p=${ps[data_row-2]}")
  ax.set_title('length vs GR')
  ax.set_xlabel('length')
  ax.set_ylabel('$¥alpha_p$')
  ax.axhline(np.log(2),ls='--')
  ax.legend(loc='upper left')
  
  plt.savefig(f'GR_vs_length_max={len(lengths)}')
  plt.show()
  
def outputGR_vs_p():
  fig,ax = plt.subplots()
  ax.plot(ps,conGRs)
  ax.set_title('p vs GR')
  ax.set_xlabel('p')
  ax.set_ylabel('$¥alpha_p$')
  
  plt.savefig(f'GR_vs_p_max={len(lengths)}')
  plt.show()
  
  
outputGR_vs_len(5)
outputGR_vs_p()