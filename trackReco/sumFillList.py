import pandas as pd

ds = 'Run-1a'
df = df = pd.read_csv('fills_'+ds+'.csv')

print(df.sum())

