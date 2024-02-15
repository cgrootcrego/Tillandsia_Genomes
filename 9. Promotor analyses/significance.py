import pandas as pd
import scipy.stats as stats
import sys

# Function to check normality and apply the appropriate test
def compare_distributions(background, genes_of_interest, column):
    # Shapiro-Wilk normality test
    bg_pvalue = stats.shapiro(background[column])[1]
    goi_pvalue = stats.shapiro(genes_of_interest[column])[1]

    # Normality threshold (alpha level)
    alpha = 0.05

    # If both groups pass the normality test, perform t-test; otherwise, perform Wilcoxon rank-sum test
    if bg_pvalue > alpha and goi_pvalue > alpha:
        stat, pvalue = stats.ttest_ind(background[column], genes_of_interest[column], equal_var=False)
        test_used = "t-test"
    else:
        stat, pvalue = stats.ranksums(background[column], genes_of_interest[column])
        test_used = "Wilcoxon rank-sum test"

    return test_used, pvalue

# Function to perform Chi-square test on a 2x2 contingency table
def chi_square_test(background, genes_of_interest, column):
    a = ((background[column] > 0).sum())
    b = ((background[column] == 0).sum())
    c = ((genes_of_interest[column] > 0).sum())
    d = ((genes_of_interest[column] == 0).sum())

    contingency_table = [[a, b], [c, d]]
    chi2, p, dof, expected = stats.chi2_contingency(contingency_table)
    return chi2, p, contingency_table

# Main function to compare background set and genes of interest
def main(background_file, genes_of_interest_file):
    # Read the files
    background = pd.read_csv(background_file, sep="\t", index_col=0)
    genes_of_interest = pd.read_csv(genes_of_interest_file, sep="\t", index_col=0)

    # Iterate through each column
    for column in background.columns:
        bg_median = background[column].median()
        bg_mean = background[column].mean()
        goi_median = genes_of_interest[column].median()
        goi_mean = genes_of_interest[column].mean()

        test_used, pvalue = compare_distributions(background, genes_of_interest, column)
        chi2, chi_p, contingency_table = chi_square_test(background, genes_of_interest, column)

        # Output results
        print(f"For column '{column}':")
        print(f"Background median: {bg_median}, mean: {bg_mean}")
        print(f"Genes of Interest median: {goi_median}, mean: {goi_mean}")
        print(f"Test used for comparing distributions: {test_used}")
        print(f"P-value from {test_used}: {pvalue}")
        print("Contingency Table (Counts > 0 vs Counts = 0):")
        print(pd.DataFrame(contingency_table, columns=['Count > 0', 'Count = 0'],
                           index=['Background', 'Genes of Interest']))
        print(f"Chi-square statistic: {chi2}, P-value from Chi-square test: {chi_p}\n")

# Retrieve file paths from command-line arguments
if len(sys.argv) < 3:
    print("Usage: python script.py <background_file_path> <genes_of_interest_file_path>")
    sys.exit(1)

background_file_path = sys.argv[1]
genes_of_interest_file_path = sys.argv[2]

# Call the main function
main(background_file_path, genes_of_interest_file_path)
