##import the necessary libraries and the GDSC dataset
import numpy as np
import pandas as pd
##import the data
gdsc = pd.read_excel("GDSC.xlsx")
## identify the key variables
gdsc.head()
##description of the data
gdsc.shape
gdsc.columns
gdsc.info()
gdsc.describe()
gdsc["DRUG_NAME"].nunique()
gdsc["TCGA_DESC"].nunique()
gdsc["CELL_LINE_NAME"].nunique()
## Identify the missing values
missing = gdsc.isnull().sum().sort_values(ascending=False)
missing.head(20)
## Identify the duplicates
gdsc = gdsc.drop_duplicates()
gdsc
# Check minimum and maximum AUC values
print("Minimum AUC:", gdsc["AUC"].min())
print("Maximum AUC:", gdsc["AUC"].max())
# Basic statistics
min_val = gdsc["LN_IC50"].min()
max_val = gdsc["LN_IC50"].max()
median_val = gdsc["LN_IC50"].median()

print("Minimum:", min_val)
print("Maximum:", max_val)
print("Median:", median_val)
# IC50_uM back-transformed from LN_IC50
gdsc["IC50_uM"] = np.exp(gdsc["LN_IC50"])
gdsc["CNA_bool"] = gdsc["CNA"].map({
    "Y": True,
    "N": False
})
gdsc["Gene Expression_bool"] = gdsc["Gene Expression"].map({
    "Y": True,
    "N": False
})
gdsc["Methylation_bool"] = gdsc["Methylation"].map({
    "Y": True,
    "N": False
})
gdsc
# drug sensitivity plot
plt.hist(gdsc['LN_IC50'], bins=50)
plt.xlabel('LN_IC50')
plt.ylabel('Frequency')
plt.title('Distribution of Drug Sensitivity')
plt.show()
# Calculate average sensitivity and variance for each drug
drug_stats = gdsc.groupby('DRUG_NAME').agg(
    Median_LN_IC50=('LN_IC50', 'median'),
    Median_AUC=('AUC', 'median'),
    Var_LN_IC50=('LN_IC50', 'var'),
    Cell_Line_Count=('CELL_LINE_NAME', 'count')
).reset_index()

# Top 5 most effective drugs (Lowest mean LN_IC50)
most_effective = drug_stats.sort_values(by='Median_LN_IC50').head(5)
print("Top 5 Most Effective Drugs (Lowest LN_IC50):")
print(most_effective[['DRUG_NAME', 'Median_LN_IC50', 'Median_AUC']])

# Top 5 least effective drugs (Highest mean LN_IC50)
least_effective = drug_stats.sort_values(by='Median_LN_IC50', ascending=False).head(5)
print("\nTop 5 Least Effective Drugs (Highest LN_IC50):")
print(least_effective[['DRUG_NAME', 'Median_LN_IC50', 'Median_AUC']])

# Top 5 most variable drugs (Context-dependent/Targeted treatments)
most_variable = drug_stats.sort_values(by='Var_LN_IC50', ascending=False).head(5)
print("\nTop 5 Most Variable Drugs (Potential Personalized Medicine Candidates):")
print(most_variable[['DRUG_NAME', 'Var_LN_IC50']])
# cancer cell line analysis
# Group by cancer type
cancer_summary = (
    gdsc.groupby("TCGA_DESC")
      .agg(
          Median_LN_IC50=("LN_IC50", "median"),
          n_Cell_Lines=("CELL_LINE_NAME", "nunique")
      )
      .sort_values("Median_LN_IC50")
)

cancer_summary.head(10)
# Identify Robust vs Small Samples
def sample_quality(n):
    if n < 5:
        return "⚠ Very small sample"
    elif n <= 20:
        return "Moderate"
    else:
        return "Robust"

cancer_summary["Sample_Quality"] = (
    cancer_summary["n_Cell_Lines"]
    .apply(sample_quality)
)

cancer_summary.head(10)
# Most Sensitive Cancer Types
cancer_summary.head(10)
# Most Resistant Cancer Types
cancer_summary.tail(10)
# Visualisation of drug-sensitive cancer types
import seaborn as sns
import matplotlib.pyplot as plt

top_sensitive = cancer_summary.head(10).reset_index()

plt.figure(figsize=(10,5))

sns.barplot(
    data=top_sensitive,
    x="TCGA_DESC",
    y="Median_LN_IC50"
)

plt.title("Most Drug-Sensitive Cancer Types")
plt.ylabel("Median LN_IC50")

plt.show()
# Cancer-Type-Selective Drug Responses
# group by cancer type and drug
drug_cancer = (
    gdsc.groupby(["DRUG_NAME", "TCGA_DESC"])
      .agg(
          Median_LN_IC50=("LN_IC50", "median")
      )
      .reset_index()
)
# Sensitivity range for each drug
# Large range = highly selective drug
drug_selectivity = (
    drug_cancer.groupby("DRUG_NAME")
    .agg(
        Min_LN_IC50=("Median_LN_IC50", "min"),
        Max_LN_IC50=("Median_LN_IC50", "max")
    )
)

drug_selectivity["Range"] = (
    drug_selectivity["Max_LN_IC50"]
    - drug_selectivity["Min_LN_IC50"]
)

drug_selectivity = (
    drug_selectivity
    .sort_values("Range", ascending=False)
)

drug_selectivity.head(10)
# Most and Least Sensitive Cancer Types
# Most sensitive cancer type
most_sensitive = (
    drug_cancer.loc[
        drug_cancer.groupby("DRUG_NAME")["Median_LN_IC50"].idxmin()
    ]
    .rename(columns={
        "TCGA_DESC": "Most_Sensitive_Cancer",
        "Median_LN_IC50": "Min_LN_IC50"
    })
)

# Least sensitive cancer type
least_sensitive = (
    drug_cancer.loc[
        drug_cancer.groupby("DRUG_NAME")["Median_LN_IC50"].idxmax()
    ]
    .rename(columns={
        "TCGA_DESC": "Least_Sensitive_Cancer",
        "Median_LN_IC50": "Max_LN_IC50"
    })
)
top_selective = (
    drug_selectivity
    .merge(
        most_sensitive[[
            "DRUG_NAME",
            "Most_Sensitive_Cancer",
            "Min_LN_IC50"
        ]],
        on="DRUG_NAME"
    )
    .merge(
        least_sensitive[[
            "DRUG_NAME",
            "Least_Sensitive_Cancer",
            "Max_LN_IC50"
        ]],
        on="DRUG_NAME"
    )
)

top_selective.head(10)
# add pathway information
pathway_info = (
    gdsc[["DRUG_NAME", "TARGET_PATHWAY"]]
    .drop_duplicates()
)

top_selective = top_selective.merge(
    pathway_info,
    on="DRUG_NAME",
    how="left"
)
top6 = top_selective.head(6)

top6
plt.figure(figsize=(10,5))

sns.barplot(
    data=top6.reset_index(),
    x="DRUG_NAME",
    y="Range"
)

plt.xticks(rotation=45)

plt.title("Top Cancer-Type-Selective Drugs")

plt.ylabel("LN_IC50 Range")

plt.show()
# Cancer Type Sensitivity Boxplot
plt.figure(figsize=(14,6))

sns.boxplot(
    data=gdsc,
    x="TCGA_DESC",
    y="LN_IC50"
)

plt.xticks(rotation=90)

plt.title("Cancer Type Sensitivity Distribution")

plt.show()
# drug-cancer type heatmap
heatmap_data = (
    drug_cancer
    .pivot(
        index="DRUG_NAME",
        columns="TCGA_DESC",
        values="Median_LN_IC50"
    )
)

plt.figure(figsize=(14,10))

sns.heatmap(
    heatmap_data,
    cmap="coolwarm",
    center=0
)

plt.title("Drug × Cancer Type Sensitivity Heatmap")

plt.show()
# Genomic & Molecular Influences on Drug Response
# LN_IC50 vs AUC Correlation
corr_gdsc = gdsc[["LN_IC50", "AUC"]].dropna()

# Sample dataset of 5000
sample_df = corr_gdsc.sample(5000, random_state=42)

from scipy.stats import pearsonr, spearmanr

# Pearson correlation
pearson_r, pearson_p = pearsonr(
    sample_df["LN_IC50"],
    sample_df["AUC"]
)

# Spearman correlation
spearman_rho, spearman_p = spearmanr(
    sample_df["LN_IC50"],
    sample_df["AUC"]
)

print("Pearson r =", round(pearson_r, 3))
print("Spearman rho =", round(spearman_rho, 3))
# Visualisation of the relationship
gdsc["MSI_STATUS"].value_counts()
# MSI Status
gdsc["Microsatellite instability Status (MSI)"].value_counts()
# Create two groups
msi_h = gdsc[gdsc["Microsatellite instability Status (MSI)"] == "MSI-H"]["LN_IC50"]

mss = gdsc[
    gdsc["Microsatellite instability Status (MSI)"]=="MSS/MSI-L"]["LN_IC50"]
# calculate medians
print("MSI-H median:", msi_h.median())
print("MSS/MSI-L median:", mss.median())

#Mann–Whitney U Test
from scipy.stats import mannwhitneyu

stat, p = mannwhitneyu(
    msi_h,
    mss,
    alternative="two-sided"
)

print("P-value:", p)
# MSI sensitivity plot
import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(6,5))

sns.boxplot(
    data= gdsc,
    x= "Microsatellite instability Status (MSI)",
    y="LN_IC50"
)

plt.title("Drug Sensitivity by MSI Status")

plt.show()

# Identify MSI-H-Selective Drugs
from scipy.stats import mannwhitneyu
import pandas as pd

results = []

for drug in gdsc["DRUG_NAME"].unique():

    sub = gdsc[gdsc["DRUG_NAME"] == drug]

    group1 = sub[
        sub["Microsatellite instability Status (MSI)"] == "MSI-H"
    ]["LN_IC50"]

    group2 = sub[
        sub["Microsatellite instability Status (MSI)"] == "MSS/MSI-L"]["LN_IC50"]

    # Skip small groups
    if len(group1) < 3 or len(group2) < 3:
        continue

    stat, p = mannwhitneyu(
        group1,
        group2,
        alternative="two-sided"
    )

    delta = (
        group1.median()
        - group2.median()
    )

    results.append([
        drug,
        delta,
        p
    ])

results_df = pd.DataFrame(
    results,
    columns=[
        "DRUG_NAME",
        "Delta_LN_IC50 (MSI-H - MSS median)",
        "p_value"
    ]
)

# Bonferroni Correction
n_tests = len(results_df)

results_df["bonferroni_p"] = (
    results_df["p_value"] * n_tests
)

# identify the Significant MSI-H-Selective Drugs
significant = results_df[
    (results_df["bonferroni_p"] < 0.05)
].sort_values("Delta_LN_IC50 (MSI-H - MSS median)")

# add pathways
pathways = (
    gdsc[["DRUG_NAME", "TARGET_PATHWAY"]]
    .drop_duplicates()
)

significant = significant.merge(
    pathways,
    on="DRUG_NAME",
    how="left"
)
significant.head(10)
# visualisation
top = significant.head(10)

plt.figure(figsize=(8,5))

sns.barplot(
    data=top,
    x="Delta_LN_IC50 (MSI-H - MSS median)",
    y="DRUG_NAME"
)

plt.title("MSI-H Selective Drugs")

plt.show()
#Growth Properties
gdsc["Growth Properties"].value_counts()

#create 2 groups
suspension = gdsc[
    gdsc["Growth Properties"] == "Suspension"
]["LN_IC50"]

adherent = gdsc[
    gdsc["Growth Properties"] == "Adherent"
]["LN_IC50"]

#calculate median
print("Suspension median:", suspension.median())
print("Adherent median:", adherent.median())

from scipy.stats import mannwhitneyu

stat, p = mannwhitneyu(
    suspension,
    adherent,
    alternative="two-sided"
)

print("P-value:", p)
# Identify Drugs Preferentially Effective in Suspension Cells
results = []

for drug in gdsc["DRUG_NAME"].unique():

    sub = gdsc[gdsc["DRUG_NAME"] == drug]

    group1 = sub[
        sub["Growth Properties"] == "Suspension"
    ]["LN_IC50"]

    group2 = sub[
        sub["Growth Properties"] == "Adherent"
    ]["LN_IC50"]

    # Skip small groups
    if len(group1) < 3 or len(group2) < 3:
        continue

    stat, p = mannwhitneyu(
        group1,
        group2,
        alternative="two-sided"
    )

    delta = (
        group1.median()
        - group2.median()
    )

    results.append([
        drug,
        delta,
        p
    ])

#convert into dataframe
import pandas as pd

results_df = pd.DataFrame(
    results,
    columns=[
        "DRUG_NAME",
        "Delta_LN_IC50",
        "p_value"
    ]
)
#multiple correction
results_df["bonferroni_p"] = (
    results_df["p_value"] * len(results_df)
)

#significant suspension selective drugs
significant = results_df[
    results_df["bonferroni_p"] < 0.05
]

significant = significant.sort_values(
    "Delta_LN_IC50"
)


#add pathwy information
pathways = (
    gdsc[["DRUG_NAME", "TARGET_PATHWAY"]]
    .drop_duplicates()
)

significant = significant.merge(
    pathways,
    on="DRUG_NAME",
    how="left"
)

significant.head(10)
# visualize
top = significant.head(10)

plt.figure(figsize=(8,5))

sns.barplot(
    data=top,
    x="Delta_LN_IC50",
    y="DRUG_NAME"
)

plt.title("Suspension-Selective Drugs")

plt.show()
# Genomic, transcriptomic and epigenomic influence on drug response
# inspect the columns
gdsc[["CNA", "Gene Expression", "Methylation"]].head()

# CNA availability
cna_yes = gdsc[gdsc["CNA_bool"] == True]["LN_IC50"]
cna_no = gdsc[gdsc["CNA_bool"] == False]["LN_IC50"]

print("CNA available median:", cna_yes.median())
print("No CNA median:", cna_no.median())

#statistical test
from scipy.stats import mannwhitneyu

stat, p = mannwhitneyu(
    cna_yes,
    cna_no,
    alternative="two-sided"
)

print("P-value:", p)

# effect size
delta = cna_yes.median() - cna_no.median()

print("Delta median:", delta)
#visualisation
import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(6,5))

sns.boxplot(
    data=gdsc,
    x="CNA_bool",
    y="LN_IC50"
)

plt.title("Effect of CNA Availability on Drug Sensitivity")

plt.show()
# Gene expression
expr_yes = gdsc[
    gdsc["Gene Expression_bool"] == True
]["LN_IC50"]

expr_no = gdsc[
    gdsc["Gene Expression_bool"] == False
]["LN_IC50"]

print("expr_yes median:", expr_yes.median())
print("expr_no median:", expr_no.median())
#statistical test
stat, p = mannwhitneyu(
    expr_yes,
    expr_no
)

# effect size
delta = (
    expr_yes.median()
    - expr_no.median()
)

print("P-value:", p)
print("Delta:", delta)
# visualisation
plt.figure(figsize=(6,5))

sns.boxplot(
    data=gdsc,
    x="Gene Expression_bool",
    y="LN_IC50"
)

plt.title("Effect of Gene Expression Availability")

plt.show()
# Methylation
Methylation_yes = gdsc[
    gdsc["Methylation_bool"] == True
]["LN_IC50"]

Methylation_no = gdsc[
    gdsc["Methylation_bool"] == False
]["LN_IC50"]

print("methylation_yes median:", Methylation_yes.median())
print("methylation_no median:", Methylation_no.median())

#statistical test
stat, p = mannwhitneyu(
    Methylation_yes,
    Methylation_no
)

# effect size
delta = (
   Methylation_yes.median()
    - Methylation_no.median()
)

print("P-value:", p)
print("Delta:", delta)
# visualisation
plt.figure(figsize=(6,5))

sns.boxplot(
    data=gdsc,
    x="Methylation_bool",
    y="LN_IC50"
)

plt.title("Effect of Methylation Availability on drug sensitivity")

plt.show()
# Mean Z-Score Heatmap
heatmap_data = (
    gdsc.groupby([
        "TCGA_DESC",
        "TARGET_PATHWAY"
    ])["Z_SCORE"]
    .mean()
    .reset_index()
)
heatmap_pivot = heatmap_data.pivot(
    index="TCGA_DESC",
    columns="TARGET_PATHWAY",
    values="Z_SCORE"
)
plt.figure(figsize=(14,8))

sns.heatmap(
    heatmap_pivot,
    cmap="coolwarm",
    center=0
)

plt.title("Mean Z-SCORE by Cancer Type and Pathway")

plt.show()
# Analysis of Drug response across the pathways
gdsc['TARGET_PATHWAY'].value_counts()
pathway_response = gdsc.groupby(
    'TARGET_PATHWAY'
)['LN_IC50'].mean().sort_values()

pathway_response
plt.figure(figsize=(10,6))

sns.barplot(
    x=pathway_response.values,
    y=pathway_response.index
)

plt.xlabel("Mean LN_IC50")
plt.ylabel("Target Pathway")
plt.title("Drug Sensitivity Across Target Pathways")

plt.show()
heatmap_data = gdsc.pivot_table(
    values='LN_IC50',
    index='TCGA_DESC',
    columns='TARGET_PATHWAY',
    aggfunc='mean'
)
plt.figure(figsize=(14,8))

sns.heatmap(
    heatmap_data,
    cmap='coolwarm',
    linewidths=0.5
)

plt.title("Cancer Type vs Pathway Drug Sensitivity")

plt.show()
#pathway level ranking
pathway_summary = gdsc.groupby(
    'TARGET_PATHWAY'
).agg({
    'LN_IC50': ['mean', 'median', 'std', 'count']
})

pathway_summary
mit = gdsc[gdsc["TARGET_PATHWAY"] == "Mitosis"]

mit_sorted = mit.groupby("DRUG_NAME")[["LN_IC50", "AUC"]].median().reset_index().sort_values(by="LN_IC50")

print(mit_sorted)