import pandas as pd
import numpy as np
column_data = ["Dev positions", "cX", "cY", "frequency"]
numDevs = 5
df = pd.DataFrame(data=np.empty((numDevs,len(column_data))), columns=column_data)

df.loc[2, column_data] = [(0,2),1,1,1]
print(df)