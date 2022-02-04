import pandas as pd

df = pd.read_csv(r"C:\Users\joann\PycharmProjects\magisterka\out\viruses_types.csv")
df.columns = ["HOST", "REGION", "number", "YEAR", "GROUP"]

df.drop('number',axis='columns', inplace=True)
print (df)
df = df.to_html()
with open(r"C:\Users\joann\PycharmProjects\magisterka\out\viruses_types.html", "w") as file:
    file.write(df)
