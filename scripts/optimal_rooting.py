"""
Find the optimal rooting of an unrooted tree that minimizes normalized clade distance with a given rooted tree
All rights reserved.
License: 3-Clause BSD,
see https://opensource.org/licenses/BSD-3-Clause
"""
import os

import dendropy
import argparse
import numpy as np


def optimal_rooted_tree(utree, rtree):
    rooted_candidates = []
    cd_values = []

    tree = dendropy.Tree(utree)
    for edge in tree.preorder_edge_iter():
        try:
            tree.reroot_at_edge(edge, update_bipartitions=True)
            rooted_candidates.append(dendropy.Tree(tree))
            cd_values.append(dendropy.calculate.treecompare.symmetric_difference(tree, rtree))
        except:
            continue

    return rooted_candidates, cd_values


def main(args):
    utree_path = args.utree
    rtree_path = args.rtree
    otree_path = args.otree

    tns = dendropy.TaxonNamespace()
    utree = dendropy.Tree.get(path=utree_path, schema='newick', taxon_namespace=tns, rooting="force-unrooted")
    rtree = dendropy.Tree.get(path=rtree_path, schema='newick', taxon_namespace=tns, rooting="force-rooted")
    rooted_utrees, cd_values = optimal_rooted_tree(utree, rtree)
    idx = np.argmin(cd_values)
    rooted_utrees[0].resolve_polytomies(update_bipartitions=True)
    print(idx, cd_values[idx])
    # print(dendropy.calculate.treecompare.symmetric_difference(rtree, rooted_utrees[idx]))
    with open(otree_path, 'w') as fp:
        fp.write(str(rooted_utrees[idx]) + ';\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--utree", type=str, help="unrooted tree",
                        required=True, default=None)
    parser.add_argument("-r", "--rtree", type=str, help="rooted tree",
                        required=True, default=None)
    parser.add_argument("-o", "--otree", type=str, help="rooted version of utree that minimizes nCD with rtree",
                        required=True, default=None)
    main(parser.parse_args())
