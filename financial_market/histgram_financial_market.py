import numpy as np

from matplotlib.image import NonUniformImage
import matplotlib.pyplot as plt
import pandas as pd

#usecols_list=['F','G','H','I']
#for my_usecols in usecols_list:
df=pd.read_excel('data_financial_market.xlsx',\
                      sheet_name='Sheet5',usecols='A:D')
#n, bins_edges, =plt.hist(df['Sharpe'],color='skyblue',edgecolor='red')
plt.hist(df['Sharpe'],color='skyblue',edgecolor='red')
plt.xlabel('Sharpe')
plt.savefig('./Sharpe.jpg')
plt.close()

plt.hist(df['income'],color='skyblue',edgecolor='red')
plt.xlabel('Annualized rate of return (%)')
plt.savefig('./income.jpg')
plt.close()

plt.hist(df['fluctuation'],color='skyblue',edgecolor='red')
plt.xlabel('Annualized rate of volatility (%)')
plt.savefig('./fluctuation.jpg')
plt.close()

plt.hist(df['income_extra'],color='skyblue',edgecolor='red')
plt.xlabel('Annual abnormal return')
plt.savefig('./income_extra.jpg')
plt.close()