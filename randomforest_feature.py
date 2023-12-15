import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc, accuracy_score, precision_score, recall_score, f1_score
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.metrics import accuracy_score
import numpy as np
import seaborn as sns

# read the data
xls_file_path = 'Feature_genera.xls'
csv_file_path = 'LEfSeLM_108_0.01_cutoff.csv'
output_file_path = './micro_gene_sample.csv'

# read the important features name
important_genes_df = pd.read_excel(xls_file_path)
important_genes = important_genes_df.iloc[:, 0].tolist()

# read the csv file
csv_df = pd.read_csv(csv_file_path)

# assume the first column is the micro name
filtered_csv_df = csv_df[csv_df[csv_df.columns[0]].isin(important_genes)]

# keep the first column and the last column
filtered_csv_df.to_csv(output_file_path, index=False)

# read the data
gene_features_df = pd.read_csv('micro_gene_sample.csv')  # 假设基因特征表是CSV文件
other_features_df = pd.read_csv('Final_patient_info.csv', encoding='ISO-8859-1', error_bad_lines=False)

# transpose the gene_features_df
gene_features_df = gene_features_df.set_index('Unnamed: 0').T

# extract the time and group from other_features_df
time_group_df = other_features_df[['Specimen_number', 'Time', 'group', 'Batch']]  # 假设样本名列名为 'sample_name'

# merge the gene_features_df and time_group_df
merged_df = pd.merge(gene_features_df, time_group_df, left_index=True, right_on='Specimen_number')

# remove the Specimen_number column
columns = [col for col in merged_df.columns if col != 'Specimen_number']
# combined the columns
df = merged_df[['Specimen_number'] + columns].copy()

# save the merged_df
df.to_csv('merged_features.csv', index=False)
df.loc[:, 'AI_Occurred'] = df['group'].apply(lambda x: 0 if 'noAI' in x else 1)
df.loc[:, 'LM_Occurred'] = df['group'].apply(lambda x: 0 if 'noLM' in x else 1)
df.iloc[:, 1:22] = df.iloc[:, 1:22].astype('float64')
df_cox = df.drop(columns=['group', 'Specimen_number'])

# save df
df.to_csv('df_r.csv', index=False)
X = df.drop(columns=['Specimen_number', 'Time', 'group', 'LM_Occurred'])  # 自变量，删除样本名和时间列
y = df['LM_Occurred']  # 因变量

def randomforest_kfold(X, y, n_splits, n_estimators, random_state):
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    accuracies = []
    precisions = []
    recalls = []
    f1s = []
    fprs_train = []
    tprs_train = []
    aucs_train = []
    fprs_test = []
    tprs_test = []
    aucs_test = []
    feature_importance_list = []

    # iterate over the folds
    for train_idx, test_idx in cv.split(X, y):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        # create the model
        model = RandomForestClassifier(n_estimators=n_estimators, max_depth=5, random_state=random_state,
                                       class_weight='balanced', min_samples_leaf=2)
        model.fit(X_train, y_train)

        # calculate accuracy, precision, recall and f1
        y_pred = model.predict(X_test)
        accuracies.append(accuracy_score(y_test, y_pred))
        precisions.append(precision_score(y_test, y_pred))
        recalls.append(recall_score(y_test, y_pred))
        f1s.append(f1_score(y_test, y_pred))

        # calculate auc and roc train
        y_scores_train = model.predict_proba(X_train)[:, 1]
        fpr_train, tpr_train, _ = roc_curve(y_train, y_scores_train)
        roc_auc_train = auc(fpr_train, tpr_train)
        fprs_train.append(fpr_train)
        tprs_train.append(tpr_train)
        aucs_train.append(roc_auc_train)

        # calculate auc and roc test
        y_scores_test = model.predict_proba(X_test)[:, 1]
        fpr_test, tpr_test, _ = roc_curve(y_test, y_scores_test)
        roc_auc_test = auc(fpr_test, tpr_test)
        fprs_test.append(fpr_test)
        tprs_test.append(tpr_test)
        aucs_test.append(roc_auc_test)
        # calculate feature importance
        feature_importance_list.append(model.feature_importances_)

    # calculate mean feature importance
    mean_feature_importance = np.mean(feature_importance_list, axis=0)
    features = pd.DataFrame({'Feature': X.columns, 'Importance': mean_feature_importance})


    return {
        'features': features,
        'mean_accuracy': np.mean(accuracies),
        'mean_precision': np.mean(precisions),
        'mean_recall': np.mean(recalls),
        'mean_f1': np.mean(f1s),
        'fprs_train': fprs_train,
        'tprs_train': tprs_train,
        'aucs_train': aucs_train,
        'fprs_test': fprs_test,
        'tprs_test': tprs_test,
        'aucs_test': aucs_test
    }


def plot_kfold(result, flag, color, legend, style):

    mean_fpr = np.linspace(0, 1, 100)
    tprs = []

    # Interpolation of TPRs for each fold
    for i in range(len(result['fprs_' + flag])):
        tprs.append(np.interp(mean_fpr, result['fprs_' + flag][i], result['tprs_' + flag][i]))

    # Calculation of average TPR and standard deviation
    mean_tpr = np.mean(tprs, axis=0)
    std_tpr = np.std(tprs, axis=0)

    # Compute confidence intervals
    tpr_upper = np.minimum(mean_tpr + std_tpr, 1)
    tpr_lower = np.maximum(mean_tpr - std_tpr, 0)

    # Add [0,0] at the start for all arrays
    mean_fpr = np.insert(mean_fpr, 0, 0)
    mean_tpr = np.insert(mean_tpr, 0, 0)
    tpr_lower = np.insert(tpr_lower, 0, 0)
    tpr_upper = np.insert(tpr_upper, 0, 0)

    # print acuracy recall precision f1
    print(f"Accuracy: {np.mean(result['mean_accuracy']):.2f} ± {np.std(result['mean_accuracy']):.2f}")
    print(f"Precision: {np.mean(result['mean_precision']):.2f} ± {np.std(result['mean_precision']):.2f}")
    print(f"Recall: {np.mean(result['mean_recall']):.2f} ± {np.std(result['mean_recall']):.2f}")
    print(f"F1: {np.mean(result['mean_f1']):.2f} ± {np.std(result['mean_f1']):.2f}")

    # Plotting the mean ROC curve
    plt.plot(mean_fpr, mean_tpr, color=color, linestyle=style, lw=2,
             label=f'{legend} (Auc = {np.mean(result["aucs_" + flag]):.2f})')

    # Plotting confidence intervals
    plt.fill_between(mean_fpr, tpr_lower, tpr_upper, color=color, alpha=0.2)


def randomforest_kfold_split(df):

    dfs = {}

    # 拆分数据集
    dfs['df_0'] = df
    dfs['df_1'] = df[df['Batch'] == 'Union74']
    dfs['df_2'] = df[df['Batch'] == 'Union34']
    dfs['df_3'] = df[df['Time'] <= 24]
    dfs['df_4'] = df[(df['Time'] > 24) & (df['Time'] <= 36)]
    dfs['df_5'] = df[(df['Time'] > 36) & (df['Time'] <= 48)]
    dfs['df_6'] = df[df['Time'] > 48]
    row_name = df.drop(columns=['Specimen_number', 'Time', 'group', 'LM_Occurred', 'Batch']).columns.tolist()

    feature_importances = {}
    # read hyperparameters from txt file
    with open('hyperparameters2.txt', 'r') as f:
        lines = f.readlines()

    hyperparameters = []
    for line in lines:
        # Split the line by comma and remove '\n', then convert each element to int
        processed_line = tuple(int(x) for x in line.strip().split(','))
        hyperparameters.append(processed_line)

    colors = [['blue', 'green', 'red'], ['yellow', 'purple', 'orange']]  # 不同数据集的颜色
    legend = [['Overall Training Dataset ROC Curve', 'Cohort 1 - Training Dataset ROC Curve',
               'Cohort 2 - Training Dataset ROC Curve'],
              ['Overall Testing Dataset ROC Curve', 'Cohort 1 - Testing Dataset ROC Curve',
               'Cohort 2 - Testing Dataset ROC Curve']]  # 不同数据集的图例

    for i in range(3):
        df_i = dfs[f'df_{i}']
        X_sub = df_i.drop(columns=['Specimen_number', 'Time', 'group', 'LM_Occurred', 'Batch'])
        y_sub = df_i['LM_Occurred']
        result = randomforest_kfold(X_sub, y_sub, n_splits=5, n_estimators=80, random_state=46)
        flag1 = 'train'
        flag2 = 'test'
        plot_kfold(result, flag1, colors[0][i], legend[0][i], style='solid')
        plot_kfold(result, flag2, colors[1][i], legend[1][i], style='dashed')
        feature_importances[f'DataSet{i}'] = result['features']['Importance']

    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curves for RandomForest in LM Classification with Batch effect removed Data', fontsize=10)
    plt.legend(loc="lower right", fontsize=6.7)
    plt.savefig('ROC_Curves_for_RandomForest_in_LM_Classification_with_Batch_Effect_Removed_Data.pdf', dpi=1000)
    plt.show()

    # combine all the features importance
    all_features_df = pd.concat(feature_importances, axis=1)
    # rename the columns and index
    all_features_df.columns = ['Overall', 'Cohort 1', 'Cohort 2']
    all_features_df.index = row_name

    # plot the heatmap
    plt.figure(figsize=(20, 8))
    sns.heatmap(all_features_df, annot=True, cmap='viridis')
    plt.title('Feature Importance Across Different data subgroups')
    plt.ylabel('Feature')
    plt.xlabel('Groups')
    plt.show()
    # save the data of feature importance
    all_features_df.to_csv('feature_importance_2(batch effect removed).csv')


def generate(df):

    # read the data
    column_names = df.columns.tolist()
    # remove the AI_Occurred column
    column_names.remove('AI_Occurred')
    # insert the AI_Occurred column to the last
    new_position = 24
    column_names.insert(new_position, 'AI_Occurred')

    df = df[column_names]
    total_non_zero_counts = (df.iloc[:, 1:25] != 0).sum()

    total_proportions = (total_non_zero_counts / len(df)).round(4)


    subsets = [
        df[df['Batch'] == 'Union74'],
        df[df['Batch'] == 'Union34'],
        df[df['Time'] <= 24],
        df[(df['Time'] > 24) & (df['Time'] <= 36)],
        df[(df['Time'] > 36) & (df['Time'] <= 48)],
        df[df['Time'] > 48]
    ]

    proportions = {'Overall': total_proportions}


    subset_labels = ['Batch1 (Union74)', 'Batch2(Union34)', '2-Year', '3-Year', '4-Year', '5.5-Year']


    for i, subset in enumerate(subsets):
        filtered_non_zero_counts = (subset.iloc[:, 1:25] != 0).sum()
        filtered_proportions = (filtered_non_zero_counts / len(subset)).round(4)
        proportions[subset_labels[i]] = filtered_proportions


    proportions_df = pd.DataFrame(proportions)


    plt.figure(figsize=(15, 8))
    sns.heatmap(proportions_df, cmap='viridis', annot=True)
    plt.title('Gene Prevelence in Different Subsets')
    plt.ylabel('Gene')
    plt.xlabel('Sample Type')
    plt.show()

    proportions_df.to_csv('proportions_2.csv')



def drug_resist(df):
    # create a empty dictionary to store different dataframes
    dfs = {}

    # split the dataset
    dfs['df_0'] = df[df['AI_Occurred'] == 1]
    dfs['df_1'] = df[df['AI_Occurred'] == 0]

    row_name = df.drop(
        columns=['Specimen_number', 'Time', 'group', 'LM_Occurred', 'Batch', 'AI_Occurred']).columns.tolist()

    feature_importances = {}
    # read hyperparameters from txt file
    with open('hyperparameters2.txt', 'r') as f:
        lines = f.readlines()

    hyperparameters = []
    for line in lines:
        # Split the line by comma and remove '\n', then convert each element to int
        processed_line = tuple(int(x) for x in line.strip().split(','))
        hyperparameters.append(processed_line)


    colors = [['blue', 'green', 'red'], ['yellow', 'purple', 'orange']]  # 不同数据集的颜色
    legend = [['AI group Training Dataset ROC Curve', 'Non-AI group Training Dataset ROC Curve'],
              ['AI group Testing Dataset ROC Curve', 'Non-AI group Testing Dataset ROC Curve']]  # 不同数据集的图例

    for i in range(2):
        df_i = dfs[f'df_{i}']
        X_sub = df_i.drop(columns=['Specimen_number', 'Time', 'group', 'LM_Occurred', 'Batch', 'AI_Occurred'])
        y_sub = df_i['LM_Occurred']
        result = randomforest_kfold(X_sub, y_sub, n_splits=5, n_estimators=100, random_state=46)
        flag1 = 'train'
        flag2 = 'test'
        plot_kfold(result, flag1, colors[0][i], legend[0][i], style='solid')
        plot_kfold(result, flag2, colors[1][i], legend[1][i], style='dashed')
        feature_importances[f'DataSet{i}'] = result['features']['Importance']

    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curves for RandomForest in LM Classification with Batch Effect Removed Data', fontsize=10)
    plt.legend(loc="lower right", fontsize=6.7)
    plt.savefig('ROC_Curves_for_RandomForest_in_LM_Classification_with_batch_effected_Data.pdf', dpi=1000)
    plt.show()



if __name__ == '__main__':

    #randomforest_kfold_split(df)
    drug_resist(df)



