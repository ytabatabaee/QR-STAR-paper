import sys
import dendropy
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

plt.rcParams["figure.figsize"] = (8, 6)
#plt.rcParams["figure.figsize"] = (42, 35)
#plt.rcParams["figure.figsize"] = (15, 12)
#plt.rcParams["figure.figsize"] = (40, 35)
np.set_printoptions(threshold=sys.maxsize)
tns = dendropy.TaxonNamespace()
caterpillars = dendropy.TreeList.get(path='topologies/caterpillar.tre', schema='newick', taxon_namespace=tns)
pseudo_caterpillars = dendropy.TreeList.get(path='topologies/pseudo_caterpillar.tre', schema='newick',
                                            taxon_namespace=tns)
balanced = dendropy.TreeList.get(path='topologies/balanced.tre', schema='newick', taxon_namespace=tns)


def clade_distance(t1, t2):
    t1.encode_bipartitions()
    t2.encode_bipartitions()
    return dendropy.calculate.treecompare.symmetric_difference(t1, t2)


def main():
    rooted_quintet_indices = np.load('rooted_quintet_indices.npy')
    cd = np.zeros((105, 105), dtype=int)

    u_distribution = np.array([5, 3, 3, 4, 2, 2, 1, 1, 2, 1, 1, 2, 4, 1, 1])

    ineqs = []
    invs = []

    for i in range(len(caterpillars)):
        inv, ineq = score(u_distribution, rooted_quintet_indices[i], "c")
        ineqs.append(ineq)
        invs.append(inv)

    for i in range(len(pseudo_caterpillars)):
        inv, ineq = score(u_distribution, rooted_quintet_indices[i + 60], "p")
        ineqs.append(ineq)
        invs.append(inv)

    for i in range(len(balanced)):
        inv, ineq = score(u_distribution, rooted_quintet_indices[i + 75], "b")
        ineqs.append(ineq)
        invs.append(inv)

    tree_pair_df = pd.DataFrame(columns=['t1', 't2', 'rd', 'V(t1,t2)'])
    count_diff = np.zeros((105, 105), dtype=int)
    count_diff_7 = np.zeros((7, 7), dtype=int)

    for i in range(105):
        t1 = get_tree(i)
        for j in range(105):
            t2 = get_tree(j)
            tree_pair_df = tree_pair_df.append({'t1': i, 't2': j, 'rd': clade_distance(t1, t2) / 2,
                                                'V(t1,t2)': len(violates(ineqs[i], ineqs[j], invs[j]))},
                                               ignore_index=True)
            cd[i][j] = clade_distance(t1, t2)  # won't need this and the line after this with the dataframe
            count_diff[i][j] = len(violates(ineqs[i], ineqs[j], invs[j]))

    idx_unrooted = [1 - 1, 2 - 1, 59 - 1, 60 - 1, 67 - 1, 76 - 1, 105 - 1]
    for i in range(7):
        for j in range(7):
            count_diff_7[i][j] = len(violates(ineqs[idx_unrooted[i]], ineqs[idx_unrooted[j]], invs[idx_unrooted[j]]))

    count_diff_c = np.zeros((60, 60), dtype=int)
    for i in range(60):
        for j in range(60):
            count_diff_c[i][j] = len(violates(ineqs[i], ineqs[j], invs[j]))

    count_diff_p = np.zeros((15, 15), dtype=int)
    for i in range(60, 75):
        for j in range(60, 75):
            count_diff_p[i - 60][j - 60] = len(violates(ineqs[i], ineqs[j], invs[j]))

    count_diff_b = np.zeros((30, 30), dtype=int)
    for i in range(75, 105):
        for j in range(75, 105):
            count_diff_b[i - 75][j - 75] = len(violates(ineqs[i], ineqs[j], invs[j]))

    #print(violates(ineqs[80], ineqs[78], invs[78]))
    #plot_conflict_heatmap(count_diff_7, labels=[1, 2, 59, 60, 67, 76, 105])
    plot_conflict_heatmap(count_diff_p, labels=list(range(61, 76, 1)))
    #plot_conflict_heatmap(count_diff, labels=list(range(1, 106, 1)))#labels=[1, 2, 59, 60, 67, 76, 105])
    #plot_clade_vs_conflicts(tree_pair_df.groupby(['rd', 'V_t1(t2)']).size().reset_index().rename(columns={0:'count'}))
    #plot_dist_conflicts(count_diff)
    #plot_path_length_parameter()


def get_tree(idx):
    if idx < 60:
        return caterpillars[idx]
    elif idx >= 60 and idx < 75:
        return pseudo_caterpillars[idx - 60]
    elif idx >= 75:
        return balanced[idx - 75]


def plot_path_length_parameter():
    plt.cla()
    x = np.linspace(0, 3, 100)
    y = np.minimum(np.exp(-x), 1-np.exp(-x))
    plt.axhline(y=0.5, color="black", linestyle=":")
    plt.axvline(color="grey")
    plt.plot(x, y, linewidth=2, label=r"$min(e^{-x}, 1-e^{-x})$")
    plt.xlim(0, 3)
    plt.ylim(0, 1.1)
    plt.grid(linestyle='--', linewidth=0.5)
    plt.title('$f(R)$')
    plt.legend()
    plt.xlabel('$x$')
    plt.savefig('f_r.pdf', bbox_inches='tight')


# def plot_A_k():

def plot_dist_conflicts(count_diff):
    plt.cla()
    ax = sns.displot(data=count_diff.flatten(), linewidth=0, color='black', kde=True, bins= 30, edgecolor='black')
    plt.grid(linestyle='--', linewidth=0.5)
    plt.title('Distribution of $V(R,R^\prime)$')
    plt.xlabel('$V(R,R^\prime)$')
    plt.ylabel('count')
    plt.savefig('dist.pdf', bbox_inches = 'tight')


def plot_clade_vs_conflicts(df):
    plt.cla()
    df = df.iloc[1:, :]
    ax = sns.scatterplot(data=df, y="rd", x="V(t1,t2)", linewidth=0, color='black', size='count', sizes=(10, 200))
    plt.grid(linestyle='--', linewidth=0.5)
    ax.set_ylim(0, 4)
    # sns.scatterplot(data=degree_dist_in, y="count", x="degree", linewidth=0, color='turquoise', label='in-degree')
    # sns.scatterplot(data=degree_dist_out, y="count", x="degree", linewidth=0, color='pink', label='out-degree')
    plt.title('Distribution of Root Distance based on $V(R,R^\prime)$')
    plt.xlabel('$V(R,R^\prime)$')
    plt.ylabel('root distance')
    plt.savefig('cd_dist2.pdf', bbox_inches = 'tight')


def plot_conflict_heatmap(count_diff, labels=[1, 2, 59, 60, 67, 76, 105]):
    plt.cla()
    ax = sns.heatmap(count_diff, annot=True, fmt="d", cmap="YlGnBu", cbar=True, cbar_kws={"shrink": .2})
    ax.set(xlabel="Other Tree (R')", ylabel='Model Tree (R)')
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.tick_params(length=0, bottom=True, top=True, right=True, left=True, labelright=True, labeltop=True,
                   labelbottom=True, labelleft=True)
    ax.set_xticks([x+0.5 for x in range(0, len(labels), 1)])
    #ax.set_xticks([x for x in range(0, len(labels), 1)])
    ax.set_yticks([x+0.5 for x in range(0, len(labels), 1)])
    #ax.set_xticks([x for x in range(0, len(labels), 1)])
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)
    #plt.legend()
    plt.savefig('heatmap_p_noinv.pdf', bbox_inches='tight')


def violates(ineqs0, ineqs1, invs1):
    # 0 is the model tree
    # 1 is the other tree
    violations = []
    for i0 in ineqs0:
        for i1 in ineqs1:
            if i0 == i1[::-1]:
                violations.append(i0)
        #for j1 in invs1:
        #    if i0 == j1:
        #        violations.append(i0)
    return violations


def invariants(u, indices, type):
    equivalence_classes = []
    inequality_classes = []
    if type == 'c':
        equivalence_classes = [[0], [1], [2], [3, 12], [4, 11], [5, 8], [6, 7, 9, 10, 13, 14]]
        inequality_classes = [[0, 1], [0, 3], [1, 4], [3, 4], [4, 6], [2, 1], [2, 5],
                              [5, 4]]  # [a, b] -> a > b, a and b are cluster indices
    elif type == 'b':
        equivalence_classes = [[0], [1, 2], [3, 12], [4, 5, 8, 11], [6, 7, 9, 10, 13, 14]]
        inequality_classes = [[0, 1], [0, 2], [1, 3], [2, 3], [3, 4]]
    elif type == 'p':
        equivalence_classes = [[0], [1, 2], [3, 12], [7, 10], [4, 5, 6, 8, 9, 11, 13, 14]]
        inequality_classes = [[0, 1], [0, 2], [0, 3], [1, 4], [2, 4], [3, 4]]

    invs = []
    for c in equivalence_classes:
        # print([indices[i]+1 for i in c]) #** important
        for i in range(len(c)):
            for j in range(len(c)):
                if i != j:
                    invs.append((indices[c[i]] + 1, indices[c[j]] + 1))

    ineqs = []
    for ineq in inequality_classes:  # distance between clusters
        for i in equivalence_classes[ineq[0]]:
            for j in equivalence_classes[ineq[1]]:
                ineqs.append((indices[i] + 1, indices[j] + 1))

    return invs, ineqs


def score(u_distribution, indices, type):
    return invariants(u_distribution, indices, type)


if __name__ == "__main__":
    main()
