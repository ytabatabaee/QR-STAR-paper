import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def plot_experiment_coef(df):
    df['cd'] = df['cd'] / (2 * df['n'] - 4)
    df['qrcoef'] = df['qrcoef'].apply(lambda x: str("{:.0E}".format(x)))
    sns.set_theme()
    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    df_acc = df.loc[df['k'] == 1000].groupby(['length', 'qrcoef'], as_index=False, sort=False)[
        'acc'].mean()
    df_nclade = df.loc[df['k'] == 1000].groupby(['length', 'qrcoef'], as_index=False, sort=False)[
        'cd'].mean()
    sns.lineplot(ax=axes[0], data=df_acc, x="qrcoef", y="acc",
                 hue='length', palette='bright')
    sns.lineplot(ax=axes[0], data=df_acc, x="qrcoef", y="acc",
                 color='black', err_style='bars', ci=None, linestyle='-.')
    sns.lineplot(ax=axes[1], data=df_nclade, x="qrcoef", y="cd",
                 hue='length', palette='bright')
    sns.lineplot(ax=axes[1], data=df_nclade, x="qrcoef", y="cd",
                 color='black', err_style='bars', ci=None, linestyle='-.')

    for i in range(2):
        axes[i].set_xlabel('Shape Coefficient ($C$)', fontsize=14)
    axes[1].set_ylabel('Rooting Error (nCD)', fontsize=14)
    axes[1].set_ylim(ymin=0, ymax=0.1)
    axes[0].set_ylabel('Proportion Correctly Rooted', fontsize=14)  # str(genes[j]) + ' Genes')
    axes[0].set_ylim(ymin=0, ymax=1)
    axes[0].text(-1.7, 1, 'a)', fontweight='bold', fontsize=18)
    axes[1].text(-1.9, 0.1, 'b)', fontweight='bold', fontsize=18)


    plt.legend(title='Seq Length')
    axes[0].get_legend().remove()
    # fig.suptitle(r'ASTRAL-III S100 (n=101, AD=0.46), QR* (LE) - Experiment on changing the Shape Coefficient ($C$) - 1000 Genes')
    plt.savefig("s100_qrstar_le_exp_coef.pdf", bbox_inches='tight')


def plot_experiment_abratio(df):
    df['cd'] = df['cd'] / (2 * df['n'] - 4)
    df['abratio'] = df['abratio'].apply(lambda x: str("{:.0E}".format(x)))
    sns.set_theme()
    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    df_acc = df.loc[df['k'] == 1000].groupby(['length', 'abratio'], as_index=False, sort=False)[
        'acc'].mean()
    df_nclade = df.loc[df['k'] == 1000].groupby(['length', 'abratio'], as_index=False, sort=False)[
        'cd'].mean()
    sns.lineplot(ax=axes[0], data=df_acc, x="abratio", y="acc",
                 hue='length', palette='bright')
    sns.lineplot(ax=axes[0], data=df_acc, x="abratio", y="acc",
                 color='black', err_style='bars', ci=None, linestyle='-.')
    sns.lineplot(ax=axes[1], data=df_nclade, x="abratio", y="cd",
                 hue='length', palette='bright')
    sns.lineplot(ax=axes[1], data=df_nclade, x="abratio", y="cd",
                 color='black', err_style='bars', ci=None, linestyle='-.')

    for i in range(2):
        axes[i].set_xlabel(r'$\frac{\alpha_{max}}{\beta_{min}}$', fontsize=20)
    axes[1].set_ylabel('Rooting Error (nCD)', fontsize=14)
    axes[1].set_ylim(ymin=0, ymax=0.1)
    axes[0].set_ylabel('Proportion Correctly Rooted', fontsize=14)  # str(genes[j]) + ' Genes')
    axes[0].set_ylim(ymin=0, ymax=1)
    axes[0].text(-1.7, 1, 'a)', fontweight='bold', fontsize=18)
    axes[1].text(-1.9, 0.1, 'b)', fontweight='bold', fontsize=18)

    plt.legend(title='Seq Length')
    axes[0].get_legend().remove()
    # fig.suptitle(r'ASTRAL-III S100 (n=101, AD=0.46), QR* (LE) - Experiment on changing $\frac{\alpha_{max}}{\beta_{min}}$ - 1000 Genes')
    plt.savefig("s100_qrstar_le_exp_abratio.pdf", bbox_inches='tight')


def plot_experiment_200taxa_trues_trueg(df):
    df = df.replace('QR*', 'QR-STAR')
    df['cd'] = df['cd'] / (2 * df['n'] - 4)
    fig = plt.figure(figsize=(20, 20), dpi=80)
    sns.set_theme()
    g = sns.FacetGrid(df.loc[(df['n'] == 201)], col='height',
                      row='speciation', margin_titles=True, aspect=0.8)
    g.map(sns.boxplot, "k", "cd", "method", palette='Blues', meanprops={"marker": "o",
                                                                        "markerfacecolor": "white",
                                                                        "markeredgecolor": "black",
                                                                        "markersize": "4"}, showmeans=True)
    g.set_titles(row_template='{row_name}', col_template='{col_name}')
    # g.add_legend()
    g.add_legend(loc='lower center', ncol=3, fontsize=14)
    g.set(ylim=(0, 0.15))
    axes = g.axes.flatten()
    axes[0].set_title("500K")
    axes[1].set_title("2M")
    axes[2].set_title("10M")
    g.set_xlabels("")  # number of genes
    g.set_ylabels("")  # ("Species Tree Error (nCD)")
    g.fig.subplots_adjust(top=0.8, bottom=0.18, left=0.14)
    g.fig.text(0.34, 0.09, s='Number of Genes', fontsize=14)
    g.fig.text(0.07, 0.35, s='Rooting Error (nCD)', rotation=90, fontsize=14)
    g.fig.suptitle('True gene trees', fontsize=14, fontweight="bold", x=0.41, y=0.88)
    # plt.title('True gene trees', fontsize=14, fontweight="bold")
    g.fig.text(0.59, 0.75, s='AD = 0.09', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.59, 0.43, s='AD = 0.21', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.405, 0.75, s='AD = 0.34', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.405, 0.43, s='AD = 0.34', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.22, 0.75, s='AD = 0.68', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.22, 0.43, s='AD = 0.69', bbox=dict(facecolor='white', alpha=0.7))
    plt.savefig("qrvsqrstar_trueg_trues.pdf", bbox_inches='tight')


def plot_experiment_200taxa_trues_estg(df):
    df = df.replace('QR*', 'QR-STAR')
    df['cd'] = df['cd'] / (2 * df['n'] - 4)
    fig = plt.figure(figsize=(20, 20), dpi=80)
    sns.set_theme()
    g = sns.FacetGrid(df.loc[(df['n'] == 201)], col='height',
                      row='speciation', margin_titles=True, aspect=0.8)
    g.map(sns.boxplot, "k", "cd", "method", palette='Blues', meanprops={"marker": "o",
                                                                        "markerfacecolor": "white",
                                                                        "markeredgecolor": "black",
                                                                        "markersize": "4"}, showmeans=True)
    g.set_titles(row_template='{row_name}', col_template='{col_name}')
    g.add_legend(loc='lower center', ncol=3, fontsize=14)
    g.set(ylim=(0, 0.15))
    axes = g.axes.flatten()
    # axes[0].set_title("500K (High ILS)")
    # axes[1].set_title("2M (Medium ILS)")
    # axes[2].set_title("10M (Low ILS)")
    axes[0].set_title("500K")
    axes[1].set_title("2M")
    axes[2].set_title("10M")
    g.set_xlabels("")  # number of genes
    g.set_ylabels("")  # ("Species Tree Error (nCD)")
    g.fig.subplots_adjust(top=0.8, bottom=0.18, left=0.14)
    g.fig.text(0.34, 0.09, s='Number of Genes', fontsize=14)
    g.fig.text(0.07, 0.35, s='Rooting Error (nCD)', rotation=90, fontsize=14)
    g.fig.suptitle('Estimated gene trees', fontsize=14, fontweight="bold", x=0.41, y=0.88)
    g.fig.text(0.59, 0.75, s='AD = 0.09', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.59, 0.43, s='AD = 0.21', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.405, 0.75, s='AD = 0.34', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.405, 0.43, s='AD = 0.34', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.22, 0.75, s='AD = 0.68', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.22, 0.43, s='AD = 0.69', bbox=dict(facecolor='white', alpha=0.7))
    plt.savefig("qrvsqrstar_estg_trues_journal.pdf", bbox_inches='tight')


def lineplot_ad_200taxa_trues_estg(df):
    df = df.replace('QR*', 'QR-STAR')
    df['cd'] = df['cd'] / (2 * df['n'] - 4)
    fig = plt.figure(figsize=(5, 5), dpi=80)
    sns.set_theme()
    sns.lineplot(data=df.loc[(df['n'] == 201)], x='ad', y='cd', hue='method', palette='Dark2', ci=None, marker='o',
                 linestyle='', linewidth=0, err_style='bars', edgecolor=None)
    plt.savefig("qrvsqrstar_trueg_trues_ad_vs_ncd.pdf", bbox_inches='tight')


def lineplot_gtee_200taxa_trues_estg(df):
    df = df.replace('QR*', 'QR-STAR')
    df['cd'] = df['cd'] / (2 * df['n'] - 4)
    fig = plt.figure(figsize=(20, 20), dpi=80)
    sns.set_theme()
    g = sns.FacetGrid(df.loc[(df['n'] == 201)], col='height', margin_titles=True, height=3, aspect=2)
    g.map(sns.lineplot, "gtee", "cd", "method", palette='Dark2',
          ci=None)  # , marker='o', linestyle='', linewidth = 0, err_style='bars')
    g.set_titles(col_template='{col_name}')
    g.add_legend()
    g.set(ylim=(0, 0.075))
    axes = g.axes.flatten()
    axes[0].set_title("500K (High ILS)")
    axes[1].set_title("2M (Medium ILS)")
    axes[2].set_title("10M (Low ILS)")
    g.set_xlabels("")  # number of genes
    g.set_ylabels("")  # ("Species Tree Error (nCD)")
    g.fig.subplots_adjust(top=0.8, bottom=0.18, left=0.14)
    g.fig.text(0.45, 0.02, s='GTEE (nRF)')
    g.fig.text(0.08, 0.25, s='Rooting Error (nCD)', rotation=90)
    # sns.relplot(
    #    data=df.loc[(df['n'] == 201) & (df['k'] == 1000)],
    #    x="gtee", y="cd",
    #    hue="method", size="k", col="height",
    #    kind="line", size_order=[50, 200, 1000], palette=sns.color_palette("rocket_r", n_colors=2),
    #    height=5, aspect=.75, facet_kws=dict(sharex=True),
    # )
    plt.savefig("qrvsqrstar_estg_trues_gtee_vs_ncd.pdf", bbox_inches='tight')


def plot_experiment_avian_astral(df):
    df = df.replace('QR*', 'QR-STAR')
    df = df.replace('optimal', 'Optimal')
    df['cd'] = df['cd'] / (2 * 48 - 4)
    fig = plt.figure(figsize=(20, 20), dpi=80)
    df['seqlen'] = df['seqlen'].replace([True], 'True gene trees')
    df['seqlen'] = df['seqlen'].replace([500], 'Est. gene trees (500 bp)')
    sns.set_theme()
    g = sns.FacetGrid(df, col='ils', row='seqlen', margin_titles=True,
                      col_order=['0.5X', '1X', '2X'], row_order=['True gene trees', 'Est. gene trees (500 bp)'])
    g.map(sns.boxplot, "k", "cd", "method", palette='PuRd', meanprops={"marker": "o",
                                                                       "markerfacecolor": "white",
                                                                       "markeredgecolor": "black",
                                                                       "markersize": "4"}, showmeans=True)
    g.set_titles(row_template='{row_name}', col_template='{col_name}')
    g.add_legend(loc='lower center', ncol=3, fontsize=12)
    g.set(ylim=(0, 0.6))
    # g.set(ylim=(0, 0.4))
    axes = g.axes.flatten()
    axes[0].set_title("0.5X ILS (0.59 AD)")
    axes[1].set_title("1X ILS (0.47 AD)")
    axes[2].set_title("2X ILS (0.35 AD)")
    g.set_xlabels("")  # number of genes
    g.set_ylabels("")  # ("Species Tree Error (nCD)")
    g.fig.subplots_adjust(top=0.8, bottom=0.18, left=0.14)
    g.fig.text(0.32, 0.09, s='Number of Genes', fontsize=14)
    g.fig.suptitle('Avian', fontsize=16, fontweight="bold", x=0.41, y=0.88)
    g.fig.text(0.09, 0.25, s='Rooted Species Tree Error (nCD)', rotation=90, fontsize=14)
    g.fig.text(0.22, 0.44, s='GTEE = 0.59', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.4, 0.44, s='GTEE = 0.54', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.58, 0.44, s='GTEE = 0.57', bbox=dict(facecolor='white', alpha=0.7))
    plt.savefig("avian_astral_cut.pdf", bbox_inches='tight')


def plot_experiment_avian_seqlen(df):
    df = df.replace('QR*', 'QR-STAR')
    df = df.replace('optimal', 'Optimal')
    df['cd'] = df['cd'] / (2 * 48 - 4)
    fig = plt.figure(figsize=(20, 20), dpi=80)
    sns.set_theme()
    g = sns.FacetGrid(df.loc[(df['ils'] == '1X')], col='seqlen', margin_titles=True,
                      col_order=[True, 1500, 1000, 500, 250])
    g.map(sns.boxplot, "k", "cd", "method", palette='PuRd', meanprops={"marker": "o",
                                                                       "markerfacecolor": "white",
                                                                       "markeredgecolor": "black",
                                                                       "markersize": "4"}, showmeans=True)
    g.set_titles(row_template='{row_name}', col_template='{col_name}')
    g.add_legend(loc='upper center', ncol=3, fontsize=14)
    g.set(ylim=(0, 0.6))
    # g.set(ylim=(0, 0.4))
    axes = g.axes.flatten()
    axes[0].set_title("True gene trees")
    axes[1].set_title("1500 bp")
    axes[2].set_title("1000 bp")
    axes[3].set_title("500 bp")
    axes[4].set_title("250 bp")
    g.set_xlabels("")  # number of genes
    g.set_ylabels("")  # ("Species Tree Error (nCD)")
    g.fig.subplots_adjust(top=0.8, bottom=0.18, left=0.14)
    g.fig.text(0.4, 0.01, s='Number of Genes', fontsize=14)
    g.fig.text(0.1, 0.05, s='Rooted Species Tree Error (nCD)', rotation=90, fontsize=14)
    g.fig.text(0.195, 0.7, s='GTEE = 0.0', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.315, 0.7, s='GTEE = 0.31', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.435, 0.7, s='GTEE = 0.39', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.555, 0.7, s='GTEE = 0.54', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.675, 0.7, s='GTEE = 0.67', bbox=dict(facecolor='white', alpha=0.7))
    plt.savefig("avian_astral_seqlen.pdf", bbox_inches='tight')


def plot_experiment_mammalian(df):
    df = df.replace('QR*', 'QR-STAR')
    df = df.replace('optimal', 'Optimal')
    df['cd'] = df['cd'] / (2 * 37 - 4)
    df['seqlen'] = df['seqlen'].replace([True], 'True gene trees')
    df['seqlen'] = df['seqlen'].replace([500], 'Est. gene trees (500 bp)')
    fig = plt.figure(figsize=(20, 20), dpi=80)
    sns.set_theme()
    g = sns.FacetGrid(df.loc[(df['seqlen'] != 1000)], col='ils', row='seqlen', margin_titles=True,
                      col_order=['0.5X', '1X', '2X'], row_order=['True gene trees', 'Est. gene trees (500 bp)'], aspect=0.8)
    g.map(sns.boxplot, "k", "cd", "method", palette='YlGn', meanprops={"marker": "o",
                                                                       "markerfacecolor": "white",
                                                                       "markeredgecolor": "black",
                                                                       "markersize": "4"}, showmeans=True)
    g.set_titles(row_template='{row_name}', col_template='{col_name}')
    g.add_legend(loc='lower center', ncol=3, fontsize=12)
    g.set(ylim=(0, 0.6))
    # g.set(ylim=(0, 0.4))
    axes = g.axes.flatten()
    axes[0].set_title("0.5X ILS (0.54 AD)")
    axes[1].set_title("1X ILS (0.33 AD)")
    axes[2].set_title("2X ILS (0.18 AD)")
    g.set_xlabels("")  # number of genes
    g.set_ylabels("")  # ("Species Tree Error (nCD)")
    g.fig.subplots_adjust(top=0.8, bottom=0.18, left=0.14)
    g.fig.suptitle('Mammalian', fontsize=16, fontweight="bold", x=0.39, y=0.88)
    g.fig.text(0.3, 0.09, s='Number of Genes', fontsize=14)
    g.fig.text(0.08, 0.25, s='Rooted Species Tree Error (nCD)', rotation=90, fontsize=14)
    g.fig.text(0.19, 0.44, s='GTEE = 0.27', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.355, 0.44, s='GTEE = 0.27', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.518, 0.44, s='GTEE = 0.27', bbox=dict(facecolor='white', alpha=0.7))
    plt.savefig("mammalian_astral.pdf", bbox_inches='tight')


def plot_experiment_mammalian_seqlen(df):
    df = df.replace('QR*', 'QR-STAR')
    df = df.replace('optimal', 'Optimal')
    df['cd'] = df['cd'] / (2 * 37 - 4)
    fig = plt.figure(figsize=(20, 20), dpi=80)
    sns.set_theme()
    g = sns.FacetGrid(df.loc[(df['ils'] == '1X')], col='seqlen', margin_titles=True,
                      col_order=[True, 1000, 500])
    g.map(sns.boxplot, "k", "cd", "method", palette='YlGn', meanprops={"marker": "o",
                                                                       "markerfacecolor": "white",
                                                                       "markeredgecolor": "black",
                                                                       "markersize": "4"}, showmeans=True)
    g.set_titles(row_template='{row_name}', col_template='{col_name}')
    g.add_legend(loc='upper center', ncol=3, fontsize=14)
    g.set(ylim=(0, 0.6))
    # g.set(ylim=(0, 0.4))
    axes = g.axes.flatten()
    axes[0].set_title("True gene trees")
    axes[1].set_title("1000 bp")
    axes[2].set_title("500 bp")
    g.set_xlabels("")  # number of genes
    g.set_ylabels("")  # ("Species Tree Error (nCD)")
    g.fig.subplots_adjust(top=0.8, bottom=0.18, left=0.14)
    g.fig.text(0.33, 0.01, s='Number of Genes', fontsize=14)
    g.fig.text(0.09, 0.05, s='Rooted Species Tree Error (nCD)', rotation=90, fontsize=14)
    g.fig.text(0.22, 0.7, s='GTEE = 0.0', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.38, 0.7, s='GTEE = 0.16', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.55, 0.7, s='GTEE = 0.27', bbox=dict(facecolor='white', alpha=0.7))
    plt.savefig("mammalian_astral_seqlen.pdf", bbox_inches='tight')


def plot_experiment_200taxa_astral_estg(df):
    df = df.replace('QR*', 'QR-STAR')
    df = df.replace('optimal', 'Optimal')
    df['cd'] = df['cd'] / (2 * 201 - 4)
    fig = plt.figure(figsize=(20, 20), dpi=80)
    sns.set_theme()
    g = sns.FacetGrid(df.loc[(df['n'] == 201)], col='height',
                      row='speciation', margin_titles=True, aspect=1.3)
    g.map(sns.boxplot, "k", "cd", "method", palette='Reds', meanprops={"marker": "o",
                                                                       "markerfacecolor": "white",
                                                                       "markeredgecolor": "black",
                                                                       "markersize": "4"}, showmeans=True)
    g.set_titles(row_template='{row_name}', col_template='{col_name}')
    g.add_legend(fontsize=14)
    # g.set(ylim=(0, 0.4))
    axes = g.axes.flatten()
    axes[0].set_title("500K")
    axes[1].set_title("2M")
    axes[2].set_title("10M")
    g.set_xlabels("")  # number of genes
    g.set_ylabels("")  # ("Species Tree Error (nCD)")
    g.fig.subplots_adjust(top=0.8, bottom=0.18, left=0.14)
    g.fig.text(0.45, 0.09, s='Number of Genes', fontsize=14)
    g.fig.text(0.09, 0.28, s='Rooted Species Tree Error (nCD)', rotation=90, fontsize=14)
    g.fig.text(0.865, 0.72, s='Estimated gene trees', fontsize=14, fontweight='bold')
    g.fig.text(0.76, 0.75, s='AD = 0.09', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.76, 0.43, s='AD = 0.21', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.52, 0.75, s='AD = 0.34', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.52, 0.43, s='AD = 0.34', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.282, 0.75, s='AD = 0.68', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.282, 0.43, s='AD = 0.69', bbox=dict(facecolor='white', alpha=0.7))
    plt.savefig("qrvsqrstar_estg_astral_full.pdf", bbox_inches='tight')


def plot_experiment_200taxa_astral_trueg(df):
    df = df.replace('QR*', 'QR-STAR')
    df = df.replace('optimal', 'Optimal')
    df['cd'] = df['cd'] / (2 * df['n'] - 4)
    fig = plt.figure(figsize=(20, 20), dpi=80)
    sns.set_theme()
    g = sns.FacetGrid(df.loc[(df['n'] == 201)], col='height',
                      row='speciation', margin_titles=True, aspect=1.3)
    g.map(sns.boxplot, "k", "cd", "method", palette='Reds', meanprops={"marker": "o",
                                                                       "markerfacecolor": "white",
                                                                       "markeredgecolor": "black",
                                                                       "markersize": "4"}, showmeans=True)
    g.set_titles(row_template='{row_name}', col_template='{col_name}')
    g.add_legend(fontsize=14)
    # g.set(ylim=(0, 0.4))
    g.set(ylim=(0, 0.4))
    axes = g.axes.flatten()
    axes[0].set_title("500K")
    axes[1].set_title("2M")
    axes[2].set_title("10M")
    g.set_xlabels("")  # number of genes
    g.set_ylabels("")  # ("Species Tree Error (nCD)")
    g.fig.subplots_adjust(top=0.8, bottom=0.18, left=0.14)
    g.fig.text(0.45, 0.09, s='Number of Genes', fontsize=14)
    g.fig.text(0.09, 0.28, s='Rooted Species Tree Error (nCD)', rotation=90, fontsize=14)
    g.fig.text(0.87, 0.72, s='True gene trees', fontsize=14, fontweight='bold')
    g.fig.text(0.76, 0.75, s='AD = 0.09', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.76, 0.43, s='AD = 0.21', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.52, 0.75, s='AD = 0.34', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.52, 0.43, s='AD = 0.34', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.282, 0.75, s='AD = 0.68', bbox=dict(facecolor='white', alpha=0.7))
    g.fig.text(0.282, 0.43, s='AD = 0.69', bbox=dict(facecolor='white', alpha=0.7))
    plt.savefig("qrvsqrstar_trueg_astral.pdf", bbox_inches='tight')


def plot_experiment_varytaxa_trues_estg(df):
    df = df.replace('QR*', 'QR-STAR')
    # df['cd'] = df['cd'] / (2 * df['n'] - 4)
    fig = plt.figure(figsize=(20, 20), dpi=80)
    sns.set_theme()
    g = sns.FacetGrid(df.loc[(df['height'] == 2000000) & (df['speciation'] == 0.000001)], col='n', col_wrap=3)
    g.map(sns.boxplot, "k", "cd", "method", palette='Purples', meanprops={"marker": "o",
                                                                          "markerfacecolor": "white",
                                                                          "markeredgecolor": "black",
                                                                          "markersize": "4"}, showmeans=True)
    g.set_titles(row_template='{row_name}', col_template='{col_name}')
    g.add_legend()
    axes = g.axes.flatten()
    axes[0].set_title("10 taxa")
    axes[1].set_title("50 taxa")
    axes[2].set_title("100 taxa")
    axes[3].set_title("200 taxa")
    # axes[4].set_title("500 taxa")
    g.set_xlabels("")  # number of genes
    g.set_ylabels("")  # ("Species Tree Error (nCD)")
    g.fig.subplots_adjust(top=0.8, bottom=0.18, left=0.14)
    g.fig.text(0.45, 0.09, s='Number of Genes')
    g.fig.text(0.08, 0.35, s='Rooting Error (nCD)', rotation=90)
    # g.set(ylim=(0, 0.4))
    plt.savefig("qrvsqrstar_varytaxa_estg_trues.pdf", bbox_inches='tight')


def plot_experiment_coef2(df):
    df['cd'] = df['cd'] / (2 * df['n'] - 4)
    df['qrcoef'] = df['qrcoef'].apply(lambda x: str("{:.0E}".format(x)))
    sns.set_theme()
    fig, axes = plt.subplots(2, 4, figsize=(24, 10))
    genes = [1000, 500, 200, 50]
    for j in range(4):
        g = genes[j]
        df_acc = df.loc[df['k'] == g].groupby(['length', 'qrcoef'], as_index=False, sort=False)[
            'acc'].mean()
        df_nclade = df.loc[df['k'] == g].groupby(['length', 'qrcoef'], as_index=False, sort=False)[
            'cd'].mean()
        sns.lineplot(ax=axes[0, j], data=df_acc, x="qrcoef", y="acc",
                     hue='length', palette='bright')
        sns.lineplot(ax=axes[0, j], data=df_acc, x="qrcoef", y="acc",
                     color='black', err_style='bars', ci=None, linestyle='-.')
        sns.lineplot(ax=axes[1, j], data=df_nclade, x="qrcoef", y="cd",
                     hue='length', palette='bright')
        sns.lineplot(ax=axes[1, j], data=df_nclade, x="qrcoef", y="cd",
                     color='black', err_style='bars', ci=None, linestyle='-.')

        for i in range(2):
            axes[i, j].set_xlabel('Shape Coefficient')
        axes[1, j].set_ylabel('nCD')
        axes[1, j].set_ylim(ymin=0, ymax=0.2)
        axes[0, j].set_ylabel('Proportion Correctly Rooted')  # str(genes[j]) + ' Genes')
        axes[0, j].set_ylim(ymin=0, ymax=1)
        axes[0, j].set_title(str(g) + ' Genes')
        axes[0, 0].set_xticklabels(axes[0, 0].get_xticklabels(), rotation=90)

    fig.suptitle('ASTRAL-III S100 (n=101, AD=0.46), QR* (LE) - Experiment on Shape Coefficient')
    plt.savefig("s100_qrstar_le_exp_coef.pdf", bbox_inches='tight')


def plot_experiment_astral3(df):
    df = df.replace('QR*', 'QR-STAR')
    df = df.replace('optimal', 'Optimal')
    df['cd'] = df['cd'] / (2 * 101 - 4)
    fig = plt.figure(figsize=(20, 20), dpi=80)
    sns.set_theme()
    g = sns.FacetGrid(df, col='seqlen', margin_titles=True,
                      col_order=[1600, 800, 400, 200], aspect=1.3)
    g.map(sns.boxplot, "k", "cd", "method", palette='YlGn', meanprops={"marker": "o",
                                                                       "markerfacecolor": "white",
                                                                       "markeredgecolor": "black",
                                                                       "markersize": "4"}, showmeans=True)
    g.set_titles(row_template='{row_name}', col_template='{col_name}')
    g.add_legend(fontsize=12)
    # g.set(ylim=(0, 0.4))
    # g.set(ylim=(0, 0.4))
    axes = g.axes.flatten()
    # axes[0].set_title("500K")
    # axes[1].set_title("2M")
    # axes[2].set_title("10M")
    g.set_xlabels("")  # number of genes
    g.set_ylabels("")  # ("Species Tree Error (nCD)")
    g.fig.subplots_adjust(top=0.8, bottom=0.18, left=0.14)
    g.fig.text(0.45, 0.09, s='Number of Genes', fontsize=14)
    g.fig.text(0.09, 0.28, s='Rooted Species Tree Error (nCD)', rotation=90, fontsize=14)
    # g.fig.text(0.76, 0.75, s='AD = 0.09', bbox=dict(facecolor='white', alpha=0.7))
    # g.fig.text(0.76, 0.43, s='AD = 0.21', bbox=dict(facecolor='white', alpha=0.7))
    # g.fig.text(0.52, 0.75, s='AD = 0.34', bbox=dict(facecolor='white', alpha=0.7))
    # g.fig.text(0.52, 0.43, s='AD = 0.34', bbox=dict(facecolor='white', alpha=0.7))
    # g.fig.text(0.282, 0.75, s='AD = 0.68', bbox=dict(facecolor='white', alpha=0.7))
    # g.fig.text(0.282, 0.43, s='AD = 0.69', bbox=dict(facecolor='white', alpha=0.7))
    plt.savefig("astarl3_estg_astral.pdf", bbox_inches='tight')


if __name__ == '__main__':
    # Experiment 0
    df1 = pd.read_csv('qr-v1.2.4_astral3_coef_exp.csv')
    plot_experiment_coef(df1)
    df2 = pd.read_csv('qr-v1.2.4_astral3_abratio_exp.csv')
    plot_experiment_abratio(df2)

    # Experiment 1
    '''df1 = pd.read_csv('qr-v1.2.4_astral2_trueg_trues.csv')
    plot_experiment_200taxa_trues_trueg(df1)
    df2 = pd.read_csv('qr-v1.2.4_astral2_estg_trues.csv')
    plot_experiment_200taxa_trues_estg(df2)'''

    # Experiment 2
    '''df1 = pd.read_csv('qr-optimal_astral2_trueg_ests.csv')
    plot_experiment_200taxa_astral_trueg(df1)
    df2 = pd.read_csv('qr-optimal_astral2_estg_ests.csv')
    plot_experiment_200taxa_astral_estg(df2)'''

    # Experiment 3
    '''df_avian1 = pd.read_csv('qr-v1.2.4_avian_trueg_astral.csv')
    df_avian2 = pd.read_csv('qr-v1.2.4_avian_estg_astral.csv')
    df_avian = pd.concat([df_avian1, df_avian2])
    plot_experiment_avian_astral(df_avian)
    plot_experiment_avian_seqlen(df_avian)

    df_m1 = pd.read_csv('qr-v1.2.4_mammalian_trueg_astral.csv')
    df_m2 = pd.read_csv('qr-v1.2.4_mammalian_estg_astral.csv')
    df_mammalian = pd.concat([df_m1, df_m2])
    plot_experiment_mammalian(df_mammalian)
    plot_experiment_mammalian_seqlen(df_mammalian)'''

    # df_astral3 = pd.read_csv('qr-v1.2.4_ASTRAL-III_estg_astral.csv')
    # plot_experiment_astral3(df_astral3)
    # plot_experiment_coef(df1[['qrcoef', "replicate", "cd", "acc", "k", "length", "n"]])
