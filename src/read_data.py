import pandas as pd

df = pd.read_csv("data/triazine_data.csv")

print(df)
print("\nNumber of compounds:", len(df))
print("\nColumn names:")
print(df.columns.tolist())